# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
Main task abstraction
"""

# python modules
import ast
import datetime
import io
import itertools
import logging
import math
import multiprocessing
import operator
import os
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
from .cif_io import get_cif_block, get_cif_dict, get_symops_from_cif, parse_cif_to_structure
from .core import Substitutor, get_symmetry_operations_from_spglib
from .patches import apply_pymatgen_patches

# shry version control
try:
    from ._version import version as shry_version
except (ModuleNotFoundError, ImportError):
    shry_version = "unknown (not installed)"

apply_pymatgen_patches()


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
class RunConfig:
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


def _normalize_config(config: RunConfig) -> RunConfig:
    sample = _normalize_sample(config.sample)
    scaling_matrix = np.array(config.scaling_matrix)
    write_ewald = config.write_ewald or config.max_ewald is not None
    return RunConfig(
        **{
            **{k: getattr(config, k) for k in config.__dataclass_fields__},
            "sample": sample,
            "scaling_matrix": scaling_matrix,
            "write_ewald": write_ewald,
        }
    )


def _validate_config(config: RunConfig) -> None:
    if len(config.from_species) != len(config.to_species):
        raise RuntimeError("from_species and to_species must have the same length.")


def prepare_structures(config: RunConfig):
    structure = LabeledStructure.from_file(
        config.structure_file,
        symmetrize=config.symmetrize,
    )
    modified_structure = structure.copy()
    # Note: since we don't limit the allowable scaling_matrix,
    # changing the order of enlargement vs. replace can have different meaning.
    modified_structure.replace_species(dict(zip(config.from_species, config.to_species)))
    modified_structure *= config.scaling_matrix
    return structure, modified_structure


def build_substitutor(config: RunConfig, modified_structure):
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


class ScriptHelper:
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
        self._timestamp = datetime.datetime.now().timestamp()
        config = RunConfig(
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

        self.structure, self.modified_structure = prepare_structures(config)
        self.substitutor = build_substitutor(config, self.modified_structure)

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
        structure_file_basename = os.path.basename(self.structure_file).split(".")[0]
        return f"shry-{structure_file_basename}-{self._timestamp}"

    def save_modified_structure(self):
        """
        Save adjusted structure
        """
        if self.no_write:
            return
        os.makedirs(self._outdir, exist_ok=True)
        structure_file_basename = os.path.basename(self.structure_file).split(".")[0]
        scaling_matrix_string = "-".join(self.scaling_matrix.flatten().astype(str))
        filename = os.path.join(
            self._outdir,
            structure_file_basename + "-" + scaling_matrix_string + ".cif",
        )
        self.modified_structure.to(filename=filename, symprec=self.symprec)

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
                    space_group = list(cifwriter.cif_file.data.values())[0]["_symmetry_space_group_name_H-M"]
                    line = line + f" {space_group}"
                print(line, file=logio)

                # Convert CifWriter to string and prepare for parallel writing
                cif_string = str(cifwriter)
                filename = filenames[i] + f"_{weight}.cif"
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
    """
    Structure + CIF's _atom_site_label
    """

    def __str__(self):
        return "LabeledStructure\n" + super().__str__()

    def __mul__(self, scaling_matrix):
        """
        The parent method returns Structure instance!
        Overwrite the offending line.
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
        return self.__class__.from_sites(new_sites, charge=new_charge)

    @classmethod
    def from_file(  # pylint: disable=arguments-differ
        cls,
        filename,
        primitive=False,
        sort=False,
        merge_tol=0.0,
        symmetrize=False,
    ):
        fname = os.path.basename(filename)
        if not fnmatch(fname.lower(), "*.cif*") and not fnmatch(fname.lower(), "*.mcif*"):
            raise ValueError("LabeledStructure only accepts CIFs.")

        # Use pycifrw instead of pymatgen.Structure.from_file
        structure = parse_cif_to_structure(filename, primitive=primitive, sort=sort, merge_tol=merge_tol)

        # Convert to LabeledStructure
        instance = cls.from_sites(structure.sites)
        instance._lattice = structure.lattice

        instance.read_label(filename, symmetrize=symmetrize)
        return instance

    def replace_species(self, species_mapping):
        """
        Replace sites from species to another species,
        or labels from _atom_site_label.

        Args:
            species_mapping (dict): from-to map of the species to be replaced.
                (string) to (string).

        Raises:
            RuntimeError: If no sites matches at least one of the from_species
                specification.
        """
        for map_num, mapping in enumerate(species_mapping.items()):
            from_species, to_species = mapping
            replace = False
            to_composition = Composition.from_string(to_species)
            for site in self:
                # Find matching site label
                if from_species in site.properties["_atom_site_label"]:
                    replace = True
                    # If match by label, don't change the label!
                    # site.properties["_atom_site_label"] = tuple(sorted({to_species}))
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
        """
        Add _atom_site_label as site_properties.
        This is useful for enforcing a certain concentration over
        a group of sites, which may not necessarily consists
        of a single orbit after enlargement to a supercell.

        Args:
            cif_filename (str): Source CIF file.
            symprec (float): Precision for the symmetry operations.

        Raises:
            RuntimeError: If any sites couldn't be matched
                to one any sites defined within the CIF.
        """
        logging.info(f"Reading _atom_site_label from {cif_filename}")

        def ufloat(string):
            """Remove uncertainties notion from floats (if any)."""
            try:
                return float(string)
            except ValueError:
                return float(string.split("(")[0])

        # Use pycifrw instead of pymatgen CifParser
        cif_dict_all = get_cif_dict(cif_filename)

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

        # Merge labels to allow multiple references.
        cif_sites = []
        for coord, zipgroup in itertools.groupby(zip(coords, labels), key=lambda x: x[0]):
            labels = tuple(sorted({x[1] for x in zipgroup}))
            cif_sites.append(
                PeriodicSite(
                    "X",
                    coord,
                    self.lattice,
                    properties={"_atom_site_label": labels},
                )
            )

        # Find equivalent sites.
        if symmetrize:
            symm_ops = get_symmetry_operations_from_spglib(self, symprec, angle_tolerance)
        else:
            # Get symmetry operations from CIF
            symops_list = get_symops_from_cif(cif_dict)
            # Spacegroup symbol and number are not important here.
            symm_ops = SpacegroupOperations(0, 0, symops_list)

        coords = [x.frac_coords for x in self.sites]
        cif_coords = [x.frac_coords for x in cif_sites]

        # List of coordinates that are equivalent to this site
        o_cif_coords = [symmop.operate_multi(cif_coords) for symmop in symm_ops]
        o_cif_coords = np.swapaxes(np.stack(o_cif_coords), 0, 1)
        o_cif_coords = [np.unique(np.mod(x, 1), axis=0) for x in o_cif_coords]

        for site in tqdm.tqdm(
            self.sites,
            desc="Matching CIF labels",
            **const.TQDM_CONF,
            disable=const.DISABLE_PROGRESSBAR,
        ):
            # Use slightly larger tolerance than pymatgen default (1e-8) to handle
            # coordinate rounding by CifParser (typically rounds to ~1e-4 precision)
            equivalent = [in_coord_list_pbc(o, site.frac_coords, atol=const.DEFAULT_ATOL) for o in o_cif_coords]

            try:
                equivalent_site = cif_sites[equivalent.index(True)]
                site.properties["_atom_site_label"] = equivalent_site.properties["_atom_site_label"]
            except ValueError as exc:
                raise RuntimeError("CIF-Structure mismatch.") from exc
