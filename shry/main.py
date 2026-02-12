# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
Main task abstraction
"""

# python modules
import ast
import copy
import datetime
import io
import itertools
import logging
import math
import multiprocessing
import operator
import os
import re
import calendar
import shutil
import signal
import sys
from fnmatch import fnmatch
from dataclasses import dataclass

# python modules
import numpy as np
import spglib
import tqdm
from pymatgen.core import Composition, PeriodicSite, Structure
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupOperations
from pymatgen.util.coord import in_coord_list_pbc, lattice_points_in_supercell

# shry modules
from . import const
from .cif_io import (
    _get_cif_block,
    _get_cif_dict,
    _get_magnetic_symops_from_cif,
    _get_symops_from_cif,
    parse_cif_to_structure,
    _str2float,
    _structure_to_cif_string,
    _structure_to_mcif_string,
)
from .core import Substitutor, _get_symmetry_operations_from_spglib
from .patches import _apply_pymatgen_patches

# shry version control
try:
    from ._version import version as shry_version
except (ModuleNotFoundError, ImportError):
    shry_version = "unknown (not installed)"

_apply_pymatgen_patches()


def _chunk_generator(gen, chunk_size):
    """
    Yield chunks from a generator without loading all items into memory.

    Args:
        gen: Generator to chunk
        chunk_size: Number of items per chunk

    Yields:
        Lists of items from the generator
    """
    chunk = []
    for item in gen:
        chunk.append(item)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def _consume_patterns_chunk(patterns_chunk):
    """
    Worker function to consume patterns (for --no-write benchmark).

    Args:
        patterns_chunk: List of (aut, pattern) tuples

    Returns:
        Number of patterns processed
    """
    # Simply consume the patterns (measures pattern generation overhead)
    count = 0
    for _ in patterns_chunk:
        count += 1
    return count


def _write_cif_files(file_data_list):
    """
    Worker function to write CIF files from strings.

    Args:
        file_data_list: List of (filename, cif_string) tuples

    Returns:
        Number of files written
    """
    count = 0
    for filename, cif_string in file_data_list:
        try:
            with open(filename, "w", encoding="utf-8") as f:
                f.write(cif_string)
            count += 1
        except Exception as e:
            logging.error(f"Failed to write {filename}: {e}")
    return count


@dataclass(frozen=True)
class _RunConfig:
    structure_file: str
    from_species: list
    to_species: list
    scaling_matrix: np.ndarray
    symmetrize: bool
    sample: int | str | None
    seed: int
    symprec: float
    atol: float
    angle_tolerance: float
    dir_size: int
    write_symm: bool
    write_ewald: bool
    max_ewald: float | None
    no_write: bool
    no_dmat: bool
    t_kind: str
    no_cache: bool = False


def _math_eval(expr):
    """
    Like eval() but limit to +-/*^** expressions for security
    """
    operators = {
        ast.Add: operator.add,
        ast.Sub: operator.sub,
        ast.Mult: operator.mul,
        ast.Div: operator.truediv,
        ast.Pow: operator.pow,
        ast.BitXor: operator.xor,
        ast.USub: operator.neg,
    }

    def _tracer(node):
        if isinstance(node, ast.Num):
            return node.n
        if isinstance(node, ast.BinOp):
            return operators[type(node.op)](_tracer(node.left), _tracer(node.right))
        if isinstance(node, ast.UnaryOp):
            return operators[type(node.op)](_tracer(node.operand))
        raise TypeError("Invalid characters on math expression")

    return _tracer(ast.parse(expr, mode="eval").body)


def _normalize_sample(sample):
    if isinstance(sample, str):
        if sample == "all":
            return None
        return int(_math_eval(sample))
    return sample


def _normalize_config(config: _RunConfig) -> _RunConfig:
    sample = _normalize_sample(config.sample)
    scaling_matrix = np.array(config.scaling_matrix)
    write_ewald = config.write_ewald or config.max_ewald is not None
    return _RunConfig(
        **{
            **{k: getattr(config, k) for k in config.__dataclass_fields__},
            "sample": sample,
            "scaling_matrix": scaling_matrix,
            "write_ewald": write_ewald,
        }
    )


def _validate_config(config: _RunConfig) -> None:
    if len(config.from_species) != len(config.to_species):
        raise RuntimeError("from_species and to_species must have the same length.")


def _prepare_structures(config: _RunConfig):
    structure = LabeledStructure.from_file(
        config.structure_file,
        symmetrize=config.symmetrize,
    )

    sym_mode = getattr(structure, "properties", {}).get("_shry_symmetry_mode", "sg")
    if sym_mode == "msg":
        print("[shry] Using MSG (magnetic symmetry)")
    else:
        print("[shry] Using SG (space-group symmetry)")

    modified_structure = structure.copy()
    # Note: since we don't limit the allowable scaling_matrix,
    # changing the order of enlargement vs. replace can have different meaning.
    modified_structure.replace_species(dict(zip(config.from_species, config.to_species)))
    modified_structure *= config.scaling_matrix
    return structure, modified_structure


def _build_substitutor(config: _RunConfig, modified_structure):
    cache = False if config.no_cache else None
    return Substitutor(
        modified_structure,
        symprec=config.symprec,
        atol=config.atol,
        angle_tolerance=config.angle_tolerance,
        shuffle=config.sample is not None,
        seed=config.seed,
        no_dmat=config.no_dmat,
        t_kind=config.t_kind,
        cache=cache,
    )


class _ScriptHelper:
    """
    Combine configurations into typical workflow in single methods
    """

    def __init__(
        self,
        structure_file,
        from_species=const.DEFAULT_FROM_SPECIES,
        to_species=const.DEFAULT_TO_SPECIES,
        scaling_matrix=const.DEFAULT_SCALING_MATRIX,
        symmetrize=const.DEFAULT_SYMMETRIZE,
        sample=const.DEFAULT_SAMPLE,
        seed=const.DEFAULT_SEED,
        symprec=const.DEFAULT_SYMPREC,
        atol=const.DEFAULT_ATOL,
        angle_tolerance=const.DEFAULT_ANGLE_TOLERANCE,
        dir_size=const.DEFAULT_DIR_SIZE,
        write_symm=const.DEFAULT_WRITE_SYMM,
        write_ewald=const.DEFAULT_WRITE_EWALD,
        max_ewald=const.DEFAULT_MAX_EWALD,
        no_write=const.DEFAULT_NO_WRITE,
        no_dmat=const.DEFAULT_NO_DMAT,
        no_cache=False,
        t_kind=const.DEFAULT_T_KIND,
    ):
        now = datetime.datetime.now()
        self._timestamp = now.timestamp()
        # Human-readable time tag for output directory names.
        # Format: HHMMSS-DMY, where month is an English 3-letter abbreviation.
        # Example: 145912-12Feb26
        mon = calendar.month_abbr[now.month]
        self._time_tag = f"{now:%H%M%S}-{now.day:02d}{mon}{now:%y}"
        config = _RunConfig(
            structure_file=structure_file,
            from_species=from_species,
            to_species=to_species,
            scaling_matrix=np.array(scaling_matrix),
            symmetrize=symmetrize,
            sample=sample,
            seed=seed,
            symprec=symprec,
            atol=atol,
            angle_tolerance=angle_tolerance,
            dir_size=dir_size,
            write_symm=write_symm,
            write_ewald=write_ewald,
            max_ewald=max_ewald,
            no_write=no_write,
            no_dmat=no_dmat,
            t_kind=t_kind,
            no_cache=no_cache,
        )
        config = _normalize_config(config)
        _validate_config(config)
        self.config = config
        self._seed = config.seed

        for key in config.__dataclass_fields__:
            setattr(self, key, getattr(config, key))

        logging.info("\nRun configurations:")
        logging.info(const.HLINE)
        logging.info(self)
        logging.info(const.HLINE)

        self.structure, self.modified_structure = _prepare_structures(config)
        self.substitutor = _build_substitutor(config, self.modified_structure)

        # Freeze output directory name for this run (avoid changing if accessed later).
        structure_file_basename = os.path.basename(self.structure_file).split(".")[0]
        base = f"shry-{structure_file_basename}-{self._time_tag}"
        if os.path.exists(base):
            raise FileExistsError(f"Output directory already exists: {base}")
        self._outdir_name = base

    def __str__(self):
        string = ""
        print_format = "  * {} = {}\n"
        string += print_format.format("SHRY version", shry_version)
        string += print_format.format("structure_file", self.structure_file)
        string += print_format.format("from_species", ", ".join(map(str, self.from_species)))
        string += print_format.format("to_species", ", ".join(map(str, self.to_species)))
        string += print_format.format("scaling_matrix", np.array(self.scaling_matrix).flatten())
        string += print_format.format("symmetrize", self.symmetrize)
        string += print_format.format("sample", self.sample)
        string += print_format.format("symprec", self.symprec)
        string += print_format.format("angle_tolerance", self.angle_tolerance)
        string += print_format.format("dir_size", self.dir_size)
        string += print_format.format("write_symm", self.write_symm)
        string += print_format.format("t-kind", self.t_kind)
        return string

    @property
    def _outdir(self):
        """
        Output directory name based on input structure filename and timestamp.

        Returns:
            str: Output directory name.
        """
        return self._outdir_name

    def save_modified_structure(self):
        """
        Save adjusted structure
        """
        if self.no_write:
            return
        os.makedirs(self._outdir, exist_ok=True)
        structure_file_basename = os.path.basename(self.structure_file).split(".")[0]
        scaling_matrix_string = "-".join(self.scaling_matrix.flatten().astype(str))
        sym_mode = (getattr(self.modified_structure, "properties", None) or {}).get("_shry_symmetry_mode", "sg")
        ext = "mcif" if sym_mode == "msg" else "cif"
        filename = os.path.join(
            self._outdir,
            structure_file_basename + "-" + scaling_matrix_string + f".{ext}",
        )
        # Avoid pymatgen's CIF writer here: it runs SpacegroupAnalyzer which may
        # raise `TypeError: unhashable type: 'numpy.ndarray'` in some environments.
        if sym_mode == "msg":
            cif_text = _structure_to_mcif_string(self.modified_structure)
        else:
            cif_text = _structure_to_cif_string(self.modified_structure)
        with open(filename, "w", encoding="utf-8") as handle:
            handle.write(cif_text)

    def write(self):
        """
        Save the irreducible structures
        """
        npatterns = self.substitutor.count()
        if not npatterns:
            logging.warning("No expected patterns.")
            return
        if self.sample is not None:
            npatterns = self.sample

        if self.no_write:
            # Parallel pattern generation benchmark (no I/O)
            n_cores = multiprocessing.cpu_count()
            chunk_size = max(100, npatterns // (n_cores * 20))
            logging.info(f"Parallel pattern generation (no-write): chunk_size={chunk_size}, n_cores={n_cores}")

            pbar = tqdm.tqdm(
                total=npatterns,
                desc=f"Generating {npatterns} order structures",
                **const.TQDM_CONF,
                disable=const.DISABLE_PROGRESSBAR,
            )

            with multiprocessing.Pool(processes=n_cores) as pool:
                pattern_batch = []
                batch_results = []

                for pattern in self.substitutor.make_patterns():
                    pattern_batch.append(pattern)

                    # Submit batch when it reaches chunk_size
                    if len(pattern_batch) >= chunk_size:
                        batch_results.append(pool.apply_async(_consume_patterns_chunk, (pattern_batch,)))
                        pattern_batch = []

                    # Update progress in real-time
                    pbar.update()

                # Submit remaining patterns
                if pattern_batch:
                    batch_results.append(pool.apply_async(_consume_patterns_chunk, (pattern_batch,)))

                # Wait for all processing to complete
                for result in batch_results:
                    result.get()

            pbar.close()
            return

        # Log file stream
        logio = io.StringIO()

        def dump_log():
            """
            Dump log
            """
            logfile = os.path.join(self._outdir, "sub.log")
            encoding = getattr(io, "LOCALE_ENCODING", "utf8")
            with open(logfile, "w", encoding=encoding, errors="surrogateescape") as f:
                logio.seek(0)
                shutil.copyfileobj(logio, f)

        def signal_handler(sig, frame):
            """
            Write log when interrupted
            """
            del sig, frame
            dump_log()
            now = datetime.datetime.now()
            tz = now.astimezone().tzname()
            time_string = now.strftime("%c ") + tz
            logging.info(const.HLINE)
            logging.info("Aborted %s", time_string)
            logging.info(const.HLINE)
            sys.exit(0)

        signal.signal(signal.SIGINT, signal_handler)

        # Filenames
        # Formatting.
        def outdir(i):
            return os.path.join(self._outdir, f"slice{i // self.dir_size}")

        formula = "".join(self.modified_structure.formula.split())
        ndigits = int(math.log10(npatterns)) + 1
        index_f = "_{:0" + str(ndigits) + "d}"
        sym_mode = (getattr(self.modified_structure, "properties", None) or {}).get("_shry_symmetry_mode", "sg")
        out_ext = "mcif" if sym_mode == "msg" else "cif"
        filenames = [os.path.join(outdir(i), formula + index_f.format(i)) for i in range(npatterns)]

        # Make directories
        ndirs = npatterns // self.dir_size + 1
        for i in range(ndirs):
            os.makedirs(os.path.join(self._outdir, f"slice{i}"), exist_ok=True)

        # Save the structures
        pbar = tqdm.tqdm(
            total=npatterns,
            desc=f"Writing {npatterns} order structures",
            **const.TQDM_CONF,
            disable=const.DISABLE_PROGRESSBAR,
        )
        os.makedirs(os.path.join(self._outdir, "slice0"), exist_ok=True)

        # TODO: Refactor
        header = "N Weight Configuration"
        if self.write_ewald:
            header = header + " EwaldEnergy"
            quantities = ("cifwriter", "weight", "letter", "ewald")
        else:
            quantities = ("cifwriter", "weight", "letter")

        if self.write_symm:
            header = header + " GroupName"
            symprec = self.symprec
        else:
            symprec = None

        if self.sample is not None:
            header += f" (seed={self._seed})"

        quantities_generator = enumerate(self.substitutor.quantities(quantities, symprec))
        # Limit amount if sampled
        if self.sample is not None:
            quantities_generator = itertools.takewhile(
                lambda x: x[0] < npatterns,
                quantities_generator,
            )

        print(header, file=logio)

        # Parallel CIF writing
        n_cores = multiprocessing.cpu_count()
        chunk_size = max(50, npatterns // (n_cores * 20))  # Smaller chunks for better load balancing
        logging.info(f"Parallel CIF writing: chunk_size={chunk_size}, n_cores={n_cores}")

        with multiprocessing.Pool(processes=n_cores) as pool:
            file_batch = []
            batch_results = []

            for i, packet in quantities_generator:
                cifwriter = packet["cifwriter"]
                ewald = packet["ewald"]
                weight = packet["weight"]
                letter = packet["letter"]

                if self.max_ewald is not None and ewald is not None and ewald > self.max_ewald:
                    continue

                # Build log line
                line = f"{i} {weight} {letter}"
                if ewald is not None:
                    line = line + f" {ewald}"
                if self.write_symm:
                    if hasattr(cifwriter, "cif_file"):
                        space_group = list(cifwriter.cif_file.data.values())[0]["_symmetry_space_group_name_H-M"]
                    else:
                        space_group = "MSG" if out_ext == "mcif" else "SG"
                    line = line + f" {space_group}"
                print(line, file=logio)

                # Convert CifWriter to string and prepare for parallel writing
                cif_string = str(cifwriter)
                filename = filenames[i] + f"_{weight}.{out_ext}"
                file_batch.append((filename, cif_string))

                # Submit batch to worker pool when it reaches chunk_size
                if len(file_batch) >= chunk_size:
                    batch_results.append(pool.apply_async(_write_cif_files, (file_batch,)))
                    file_batch = []

                # Update progress as we generate (not waiting for writes)
                pbar.update()

                if i >= npatterns - 1:
                    break

            # Submit remaining files
            if file_batch:
                batch_results.append(pool.apply_async(_write_cif_files, (file_batch,)))

            # Wait for all writes to complete
            for result in batch_results:
                result.get()

        # i should be set from the loop above
        pbar.close()
        dump_log()

        # Warn if too little structures are generated
        if i < npatterns - 1:
            logging.warning(f"Too little structures generated ({i + 1}/{npatterns}). Check `atol` value.")

    def count(self) -> None:
        """
        Count the number of unique substituted structures
        """
        count = self.substitutor.count()
        total_count = self.substitutor.total_count()  # pylint: disable=assignment-from-no-return
        logging.info(const.HLINE)
        logging.info(f"Total number of combinations is {total_count}")
        logging.info(f"Expected unique patterns is {count}")
        logging.info(const.HLINE)
        return count


class LabeledStructure(Structure):
    """Structure subclass that preserves CIF label and magnetic-moment metadata.

    In addition to standard ``pymatgen.Structure`` behavior, this class stores
    per-site CIF-derived properties such as ``_atom_site_label`` and magnetic
    moment mappings. These properties are kept consistent across common
    structure operations (copy, supercell expansion, species replacement), so
    downstream SHRY substitution logic can track site identity robustly.
    """

    def __str__(self):
        """Return a human-readable representation with subclass header."""
        return "LabeledStructure\n" + super().__str__()

    def __mul__(self, scaling_matrix):
        """Return a supercell while preserving the ``LabeledStructure`` type.

        This overrides ``Structure.__mul__`` to keep the subclass and its
        structure-level properties (including SG/MSG mode and magnetic metadata)
        after supercell creation.

        Args:
            scaling_matrix (array-like): Supercell scaling in pymatgen-compatible
                shape ``(1,)``, ``(3,)``, or ``(3, 3)``.

        Returns:
            LabeledStructure: Expanded structure with copied site and structure
                properties.

        Raises:
            ValueError: If ``scaling_matrix`` has unsupported shape.
        """

        scale_matrix = np.array(scaling_matrix, np.int16)

        # check the shape of the scaling matrix
        if scale_matrix.shape not in {(1,), (3,), (3, 3)}:
            logging.warning("The scale_matrix.shape should be (1,), (3,), or (3, 3)")
            raise ValueError
        else:
            if scale_matrix.shape != (3, 3):
                scale_matrix = np.array(scale_matrix * np.eye(3), np.int16)
            else:
                pass

        new_lattice = Lattice(np.dot(scale_matrix, self._lattice.matrix))

        f_lat = lattice_points_in_supercell(scale_matrix)
        c_lat = new_lattice.get_cartesian_coords(f_lat)

        new_sites = []
        for site in self:
            for v in c_lat:
                s = PeriodicSite(
                    site.species,
                    site.coords + v,
                    new_lattice,
                    properties=site.properties,
                    coords_are_cartesian=True,
                    to_unit_cell=False,
                    skip_checks=True,
                )
                new_sites.append(s)

        new_charge = self._charge * np.linalg.det(scale_matrix) if self._charge else None
        # This line.
        new_struct = self.__class__.from_sites(new_sites, charge=new_charge)

        # Preserve Structure-level properties (pymatgen doesn't always keep these
        # through from_sites / __mul__). This is required for MSG/SG mode tracking
        # and for global magCIF label->moment maps.
        old_props = getattr(self, "properties", None) or {}
        if getattr(new_struct, "properties", None) is None:
            new_struct.properties = {}
        new_struct.properties.update(dict(old_props))
        return new_struct

    def copy(self, *args, **kwargs):  # pylint: disable=signature-differs
        """Return a deep copy preserving the ``LabeledStructure`` type.

        pymatgen's Structure.copy() returns a Structure instance; SHRY relies on
        LabeledStructure methods (e.g. replace_species) to keep CIF labels and
        magnetic moments consistent after modifications.

        Returns:
            LabeledStructure: Copied instance with duplicated site and structure
                properties.
        """
        base = super().copy(*args, **kwargs)
        # base is typically a plain Structure; convert back to our subclass.
        new_struct = self.__class__.from_sites(base.sites, charge=getattr(base, "_charge", None))
        new_struct._lattice = base.lattice
        if getattr(base, "properties", None) is None:
            new_struct.properties = {}
        else:
            new_struct.properties = copy.deepcopy(dict(base.properties))
        return new_struct

    @classmethod
    def from_file(  # pylint: disable=arguments-differ
        cls,
        filename,
        primitive=False,
        sort=False,
        merge_tol=0.0,
        symmetrize=False,
    ):
        """Create a ``LabeledStructure`` from CIF/magCIF and attach labels.

        Args:
            filename (str): Path to ``.cif`` or ``.mcif`` file.
            primitive (bool): Whether to return a primitive structure.
            sort (bool): Whether to sort sites during parsing.
            merge_tol (float): Site merge tolerance used by parser.
            symmetrize (bool): If True, match labels using spglib symmetry
                operations instead of CIF symmetry operators.

        Returns:
            LabeledStructure: Parsed structure populated with
                ``_atom_site_label`` and related properties.

        Raises:
            ValueError: If file extension is not CIF/magCIF.
        """
        fname = os.path.basename(filename)
        if not fnmatch(fname.lower(), "*.cif*") and not fnmatch(fname.lower(), "*.mcif*"):
            raise ValueError("LabeledStructure only accepts CIFs.")

        sym_mode = "msg" if fnmatch(fname.lower(), "*.mcif*") else "sg"

        # Use pycifrw instead of pymatgen.Structure.from_file
        structure = parse_cif_to_structure(filename, primitive=primitive, sort=sort, merge_tol=merge_tol)

        # Convert to LabeledStructure
        instance = cls.from_sites(structure.sites)
        instance._lattice = structure.lattice

        # Record symmetry mode decision on Structure.
        # This will be used by spglib wrappers to switch SG/MSG behavior.
        if getattr(instance, "properties", None) is None:
            instance.properties = {}
        instance.properties["_shry_symmetry_mode"] = sym_mode

        instance.read_label(filename, symmetrize=symmetrize)
        return instance

    def replace_species(self, species_mapping):
        """Replace species/labels while preserving CIF and magnetic metadata.

        Args:
            species_mapping (dict[str, str]): Mapping from source label/element to
                target species expression. Keys are first matched against
                ``_atom_site_label`` entries; if not found, element-symbol matches
                are attempted.

        Raises:
            RuntimeError: If no sites match at least one source key in
                specification.
        """

        def _derive_label(from_label: str, to_species: str) -> str:
            """Derive a reasonable target CIF label.

            If the user specifies an element symbol (e.g. "Cu") but the input CIF
            uses numbered labels (e.g. "Co1"), keep the suffix: "Cu1".
            """
            to_species_s = str(to_species).strip()
            from_label_s = str(from_label).strip()
            if re.fullmatch(r"[A-Za-z]+", to_species_s):
                m = re.match(r"^[A-Za-z]+(.*)$", from_label_s)
                suffix = m.group(1) if m else ""
                if suffix:
                    return f"{to_species_s}{suffix}"
            return to_species_s

        def _lookup_moment(label: str, element_symbol: str, moments_by_label: dict):
            if label in moments_by_label:
                return moments_by_label[label]
            if element_symbol in moments_by_label:
                return moments_by_label[element_symbol]
            # Fallback: first matching label like "Cu1", "Cu2", ...
            for k in sorted(moments_by_label.keys()):
                if isinstance(k, str) and k.startswith(element_symbol):
                    return moments_by_label[k]
            return (0.0, 0.0, 0.0)

        def _transform_base_moment_for_site(site, base_vec: tuple[float, float, float]) -> tuple[float, float, float]:
            """Transform a label's *base* moment vector into this site's orientation.

            For magCIF-derived structures we store the symmetry op used to generate each
            site (rotation + translation) and time reversal in site.properties.
            """
            op_aff = site.properties.get("_shry_magn_op_affine")
            tr = site.properties.get("_shry_magn_time_reversal")
            if op_aff is None or tr is None:
                return base_vec
            try:
                r = np.array(op_aff, dtype=float)[:3, :3]
                det = np.linalg.det(r)
                det_sign = 1.0 if det >= 0 else -1.0
                v = np.array(base_vec, dtype=float)
                vv = det_sign * (r @ v)
                vv = float(int(tr)) * vv
                return (float(vv[0]), float(vv[1]), float(vv[2]))
            except Exception:
                return base_vec

        for map_num, mapping in enumerate(species_mapping.items()):
            from_species, to_species = mapping
            replace = False
            to_composition = Composition.from_string(to_species)
            for site in self:
                # Find matching site label
                if from_species in site.properties["_atom_site_label"]:
                    replace = True

                    to_label = _derive_label(from_species, to_species)

                    # If match by label, also replace the label-associated magnetic moment.
                    old_labels = tuple(site.properties.get("_atom_site_label", ()))
                    string_labels = [to_label if x == from_species else x for x in old_labels if isinstance(x, str)]
                    string_labels = sorted(set(string_labels))
                    other_labels = [x for x in old_labels if not isinstance(x, str)]
                    site.properties["_atom_site_label"] = tuple(string_labels + other_labels)

                    old_mom_map = dict(site.properties.get("_atom_site_moment_by_label", {}))
                    old_occ_map = dict(site.properties.get("_atom_site_occupancy_by_label", {}))

                    # Merge occupancies if from->to collides with an existing label.
                    new_occ_map: dict[str, float] = {}
                    for lab, occ in old_occ_map.items():
                        new_lab = to_label if lab == from_species else lab
                        new_occ_map[new_lab] = float(new_occ_map.get(new_lab, 0.0)) + float(occ)

                    moments_by_label = getattr(self, "properties", {}).get("_atom_site_moments_by_label", {})

                    # Determine the transformed moment for to_label.
                    if to_label in old_mom_map:
                        to_vec = old_mom_map[to_label]
                    else:
                        base = _lookup_moment(to_label, str(to_species).strip(), moments_by_label)
                        to_vec = _transform_base_moment_for_site(site, base)

                    # Build new mom_map keeping existing labels (except removed one).
                    new_mom_map: dict[str, tuple[float, float, float]] = {}
                    for lab, vec in old_mom_map.items():
                        if lab == from_species:
                            continue
                        # If some other label already equals to_label, keep it (it is already transformed).
                        if lab == to_label:
                            new_mom_map[lab] = vec
                        else:
                            new_mom_map[lab] = vec
                    new_mom_map[to_label] = to_vec

                    site.properties["_atom_site_moment_by_label"] = new_mom_map
                    site.properties["_atom_site_occupancy_by_label"] = new_occ_map
                    site.properties["magmom"] = self._representative_magmom(new_mom_map, new_occ_map)

                    try:
                        site.species = to_composition
                    except ValueError:
                        site.species = to_composition.fractional_composition
                # If failed, try to find matching Element
                elif any(e.symbol == from_species for e in site.species.elements):
                    replace = True
                    # Since all sites are replaced, merge them under one label.
                    # But give a distinct ID in case two or more sites are
                    # replaced into the same species.
                    new_label = tuple(sorted({to_species}) + [map_num])
                    site.properties["_atom_site_label"] = new_label

                    moments_by_label = getattr(self, "properties", {}).get("_atom_site_moments_by_label", {})
                    mom = _lookup_moment(str(to_species).strip(), str(to_species).strip(), moments_by_label)
                    site.properties["_atom_site_moment_by_label"] = {to_species: mom}
                    site.properties["_atom_site_occupancy_by_label"] = {to_species: 1.0}
                    site.properties["magmom"] = self._representative_magmom(
                        site.properties["_atom_site_moment_by_label"], site.properties["_atom_site_occupancy_by_label"]
                    )

                    try:
                        site.species = to_composition
                    except ValueError:
                        site.species = to_composition.fractional_composition

            if not replace:
                raise RuntimeError(f"Can't find the specified site ({from_species}).")

    def read_label(
        self,
        cif_filename,
        symprec=const.DEFAULT_SYMPREC,
        angle_tolerance=const.DEFAULT_ANGLE_TOLERANCE,
        symmetrize=False,
    ):
        """Read CIF/magCIF labels and map them onto current structure sites.

        Adds ``_atom_site_label`` and associated magnetic/occupancy properties to
        each site. For magCIF input, magnetic moments are transformed according to
        magnetic symmetry operations when required.

        Args:
            cif_filename (str): Source CIF file.
            symprec (float): Symmetry precision used when ``symmetrize=True``.
            angle_tolerance (float): Angle tolerance for symmetry operations.
            symmetrize (bool): If True, use spglib-derived symmetry operations for
                label matching.

        Raises:
            RuntimeError: If any structure site cannot be matched to CIF-defined
                positions/labels.
        """
        logging.info(f"Reading _atom_site_label from {cif_filename}")

        def ufloat(string):
            """Remove uncertainties notion from floats (if any)."""
            try:
                return float(string)
            except ValueError:
                return float(string.split("(")[0])

        # Use pycifrw instead of pymatgen CifParser
        cif_dict_all = _get_cif_dict(cif_filename)

        # Get first block (same behavior as original)
        cif_dict = list(cif_dict_all.values())[0]

        def _as_list(v):
            if isinstance(v, (list, tuple)):
                return list(v)
            return [v]

        x_raw = _as_list(cif_dict["_atom_site_fract_x"])
        y_raw = _as_list(cif_dict["_atom_site_fract_y"])
        z_raw = _as_list(cif_dict["_atom_site_fract_z"])

        x_list = map(ufloat, x_raw)
        y_list = map(ufloat, y_raw)
        z_list = map(ufloat, z_raw)
        coords = [(x, y, z) for x, y, z in zip(x_list, y_list, z_list)]

        labels = _as_list(cif_dict["_atom_site_label"])

        # Optional: per-label occupancy (needed for mixed sites).
        occupancies = _as_list(cif_dict.get("_atom_site_occupancy", []))
        if not occupancies:
            occupancies = [1.0] * len(labels)

        # Parse magCIF moments if present.
        moments_by_label = {}
        if "_atom_site_moment.label" in cif_dict:
            m_labels = _as_list(cif_dict.get("_atom_site_moment.label"))
            mxs = _as_list(cif_dict.get("_atom_site_moment.crystalaxis_x", []))
            mys = _as_list(cif_dict.get("_atom_site_moment.crystalaxis_y", []))
            mzs = _as_list(cif_dict.get("_atom_site_moment.crystalaxis_z", []))
            for lab, mx, my, mz in zip(m_labels, mxs, mys, mzs, strict=False):
                try:
                    moments_by_label[str(lab).strip()] = (_str2float(mx), _str2float(my), _str2float(mz))
                except Exception:
                    continue

        if getattr(self, "properties", None) is None:
            self.properties = {}
        if moments_by_label:
            self.properties["_atom_site_moments_by_label"] = moments_by_label

        def _transform_moment(vec, op, time_reversal: int) -> tuple[float, float, float]:
            v = np.array(vec, dtype=float)
            if v.shape == ():
                v = np.array([float(v), 0.0, 0.0], dtype=float)
            if v.shape != (3,):
                return (0.0, 0.0, 0.0)
            r = np.array(op.rotation_matrix, dtype=float)
            det = np.linalg.det(r)
            det_sign = 1.0 if det >= 0 else -1.0
            vv = det_sign * (r @ v)
            vv = float(time_reversal) * vv
            return (float(vv[0]), float(vv[1]), float(vv[2]))

        # Merge labels to allow multiple references.
        # Keep per-label occupancy so we can compute representative magmom as
        # occupancy-weighted average for mixed/disordered sites.
        merged_rows = list(zip(coords, labels, occupancies, strict=False))
        merged_rows.sort(key=lambda t: t[0])
        cif_sites = []
        for coord, zipgroup in itertools.groupby(merged_rows, key=lambda x: x[0]):
            items = list(zipgroup)
            merged_labels = tuple(sorted({x[1] for x in items}))
            mom_map = {lab: moments_by_label.get(lab, (0.0, 0.0, 0.0)) for lab in merged_labels}
            occ_map = {}
            for _c, lab, occ in items:
                try:
                    occ_map[str(lab)] = occ_map.get(str(lab), 0.0) + _str2float(occ)
                except Exception:
                    occ_map[str(lab)] = occ_map.get(str(lab), 0.0)
            cif_sites.append(
                PeriodicSite(
                    "X",
                    coord,
                    self.lattice,
                    properties={
                        "_atom_site_label": merged_labels,
                        "_atom_site_moment_by_label": mom_map,
                        "_atom_site_occupancy_by_label": occ_map,
                    },
                )
            )

        # Find equivalent sites.
        # For magCIF: we must transform magnetic moments under symmetry operations.
        sym_mode = (getattr(self, "properties", None) or {}).get("_shry_symmetry_mode", "sg")
        use_mag_ops = (sym_mode == "msg") and (
            "_space_group_symop_magn_operation.xyz" in cif_dict or "_space_group_symop_magn_operation_xyz" in cif_dict
        )

        if symmetrize or not use_mag_ops:
            # Keep the previous (faster) logic for non-magnetic or spglib-symmetrized matching.
            if symmetrize:
                symm_ops = _get_symmetry_operations_from_spglib(self, symprec, angle_tolerance)
            else:
                symops_list = _get_symops_from_cif(cif_dict)
                symm_ops = SpacegroupOperations(0, 0, symops_list)

            cif_coords = [x.frac_coords for x in cif_sites]
            o_cif_coords = [symmop.operate_multi(cif_coords) for symmop in symm_ops]
            o_cif_coords = np.swapaxes(np.stack(o_cif_coords), 0, 1)
            o_cif_coords = [np.unique(np.mod(x, 1), axis=0) for x in o_cif_coords]

            for site in tqdm.tqdm(
                self.sites,
                desc="Matching CIF labels",
                **const.TQDM_CONF,
                disable=const.DISABLE_PROGRESSBAR,
            ):
                equivalent = [in_coord_list_pbc(o, site.frac_coords, atol=const.DEFAULT_ATOL) for o in o_cif_coords]
                try:
                    equivalent_site = cif_sites[equivalent.index(True)]
                    site.properties["_atom_site_label"] = equivalent_site.properties["_atom_site_label"]
                    site.properties["_atom_site_moment_by_label"] = equivalent_site.properties.get(
                        "_atom_site_moment_by_label", {}
                    )
                    site.properties["_atom_site_occupancy_by_label"] = equivalent_site.properties.get(
                        "_atom_site_occupancy_by_label", {}
                    )
                    mom_map = site.properties.get("_atom_site_moment_by_label", {})
                    occ_map = site.properties.get("_atom_site_occupancy_by_label", {})
                    site.properties["magmom"] = self._representative_magmom(mom_map, occ_map)
                except ValueError as exc:
                    raise RuntimeError("CIF-Structure mismatch.") from exc
            return

        # Magnetic matching with moment transformation.
        mag_ops = _get_magnetic_symops_from_cif(cif_dict)
        base_sites = []
        for cs in cif_sites:
            base_sites.append(
                (
                    np.array(cs.frac_coords, dtype=float),
                    tuple(cs.properties.get("_atom_site_label", ())),
                    dict(cs.properties.get("_atom_site_occupancy_by_label", {})),
                )
            )

        for site in tqdm.tqdm(
            self.sites,
            desc="Matching CIF labels (magCIF)",
            **const.TQDM_CONF,
            disable=const.DISABLE_PROGRESSBAR,
        ):
            found = False
            for base_coord, base_labels, occ_map in base_sites:
                for op, tr in mag_ops:
                    new_coord = op.operate(base_coord) % 1.0
                    if in_coord_list_pbc([new_coord], site.frac_coords, atol=const.DEFAULT_ATOL):
                        site.properties["_atom_site_label"] = base_labels
                        site.properties["_atom_site_occupancy_by_label"] = occ_map

                        # Record which magnetic symmetry operation generated this site,
                        # so we can transform moments for substituted labels later.
                        site.properties["_shry_magn_op_affine"] = np.array(op.affine_matrix, dtype=float)
                        site.properties["_shry_magn_time_reversal"] = int(tr)

                        mom_map = {}
                        for lab in base_labels:
                            vec = moments_by_label.get(str(lab).strip(), (0.0, 0.0, 0.0))
                            mom_map[str(lab)] = _transform_moment(vec, op, tr)
                        site.properties["_atom_site_moment_by_label"] = mom_map
                        site.properties["magmom"] = self._representative_magmom(mom_map, occ_map)
                        found = True
                        break
                if found:
                    break
            if not found:
                raise RuntimeError("CIF-Structure mismatch (magCIF).")

    @staticmethod
    def _representative_magmom(mom_map: dict, occ_map: dict | None = None) -> np.ndarray:
        """Compute representative magmom for spglib magnetic symmetry.

        For mixed sites (multiple labels with moments), use an occupancy-weighted
        average when per-label occupancies are available.
        """
        if not mom_map:
            return np.zeros(3, dtype=float)

        if occ_map:
            total = 0.0
            vec = np.zeros(3, dtype=float)
            for lab, m in mom_map.items():
                w = float(occ_map.get(lab, 0.0))
                if w == 0.0:
                    continue
                vec += w * np.array(m, dtype=float)
                total += w
            if total > 0.0:
                return vec / total

        for m in mom_map.values():
            v = np.array(m, dtype=float)
            if not np.allclose(v, 0.0):
                return v
        return np.zeros(3, dtype=float)
