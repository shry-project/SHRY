# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""
Command line interface.
"""

# python modules
import argparse
import datetime
import fnmatch
import logging

# python modules
import tqdm
import numpy as np

# shry modules
from . import const

# shry version control
try:
    from ._version import version as shry_version
except (ModuleNotFoundError, ImportError):
    shry_version = "unknown (not installed)"


def _setup_tqdm_logging(level=logging.INFO, fmt="%(message)s"):
    """
    Configure logging to work with tqdm without overwriting the progress bar.
    Returns a tuple needed to restore previous logging handlers.
    """

    class _TqdmLoggingHandler(logging.Handler):
        def emit(self, record):
            try:
                msg = self.format(record)
                tqdm.tqdm.write(msg)
                self.flush()
            except (KeyboardInterrupt, SystemExit):
                raise
            except Exception:  # pylint: disable=bare-except
                self.handleError(record)

    handler = _TqdmLoggingHandler()
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter(fmt))

    root = logging.getLogger()
    prev_handlers = list(root.handlers)
    prev_level = root.level
    root.handlers = []
    root.setLevel(level)
    root.addHandler(handler)
    return root, prev_handlers, prev_level


def _restore_logging(state):
    """
    Restore logging handlers and level from setup_tqdm_logging.
    """
    root, prev_handlers, prev_level = state
    root.handlers = prev_handlers
    root.setLevel(prev_level)


def _print_header():
    """
    Print header text
    """
    tz = const.NOW.astimezone().tzname()
    time_string = const.NOW.strftime("%c ") + tz

    logging.info("********************************\n")
    logging.info("SHRY: Suite for High-throughput generation of modelswith atomic substitutions implemented by Python")
    logging.info("")
    logging.info("Please cite the following paper::\n")
    logging.info("G.I. Prayogo, et al., J. Chem. Inf. Model. 62, 2909-2915 (2022)")
    logging.info("https://doi.org/10.1021/acs.jcim.2c00389")
    logging.info("\n********************************")
    logging.info(f"Begin {time_string} (unixtime: {const.DEFAULT_SEED})")


def _print_footer():
    """
    Print footer text
    """
    now = datetime.datetime.now()
    tz = now.astimezone().tzname()
    time_string = now.strftime("%c ") + tz
    logging.info(const.HLINE)
    logging.info("When you publish a paper, please cite the following paper!!")
    logging.info(const.HLINE)
    logging.info("SHRY: Application of Canonical Augmentation to the Atomic Substitution Problem")
    logging.info("G.I. Prayogo*, A. Tirelli, K. Utimula, K. Hongo, R. Maezono, and K. Nakano*,")
    logging.info("J. Chem. Inf. Model. 62, 2909-2915 (2022)")
    logging.info("")
    logging.info("@article{doi:10.1021/acs.jcim.2c00389,")
    logging.info(
        "  author = {Prayogo, Genki and Tirelli, Andrea and Utimula, Keishu and Hongo, Kenta and Maezono, Ryo and Nakano, Kousuke},"
    )
    logging.info("  title = {SHRY: Application of Canonical Augmentation to the Atomic Substitution Problem},")
    logging.info("  journal = {J. Chem. Inf. Model.},")
    logging.info("  volume = {62},")
    logging.info("  number = {12},")
    logging.info("  pages = {2909-2915},")
    logging.info("  year = {2022},")
    logging.info("  doi = {10.1021/acs.jcim.2c00389},")
    logging.info("}")
    logging.info("")
    logging.info("If SHRY helps your work, please consider giving us a star on GitHub:")
    logging.info("https://github.com/shry-project/SHRY")
    logging.info(const.HLINE)
    logging.info("Ends " + time_string)
    logging.info(const.HLINE)


def _get_parser():
    parser = argparse.ArgumentParser(
        epilog=f"SHRY version {shry_version}",
        description="Quick use: `shry STRUCTURE_CIF`. See `shry -h` for more options.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # group = parser.add_argument_group("Input")
    parser.add_argument(
        "input",
        type=str,
        help=("A CIF containing command line options as keys.\n(See `$SHRY_INSTALLDIR/examples`)"),
    )

    group = parser.add_argument_group("structure modification")
    group.add_argument(
        "--from-species",
        "-f",
        nargs="*",
        type=str,
        help=(
            "Replace FROM_SPECIES from the CIF file into TO_SPECIES. "
            "Matches either `_atom_site_label` or `atom_site_type_symbol`. "
            "Use comma or space (preferred) for multiple substitutions."
        ),
        default=const.DEFAULT_FROM_SPECIES,
    )
    group.add_argument(
        "--to-species",
        "-t",
        nargs="*",
        type=str,
        help=(
            "Final chemical formula of the replaced FROM_SPECIES. "
            "Accepts various formats for the concentration and oxidation states "
            "such as SmFe12, Sm0.5Fe0.5, Sm3+Sm2+2Fe3, etc. "
            "The number of entry must be the same as FROM_SPECIES."
        ),
        default=const.DEFAULT_TO_SPECIES,
    )
    group.add_argument(
        "--scaling-matrix",
        "-s",
        nargs="*",
        type=str,  # To allow flexible separator
        help=(
            "Three (for diagonal supercells) or nine (for non-diagonal supercells) integers specifying "
            "the scaling matrix for constructing a supercell. One scalar value is also accepted, "
            "which is converted to three scalar values (e.g., 1 -> 1 1 1)"
        ),
        default=const.DEFAULT_SCALING_MATRIX_STR,
    )

    group = parser.add_argument_group("input and output")
    group.add_argument(
        "--count-only",
        action="store_true",
        help=("Enumerate the total number of ordered configurations from the given CIF and quit."),
    )
    group.add_argument(
        "--mod-only",
        action="store_true",
        help="Write a modified CIF without generating the ordered structures.",
    )
    group.add_argument(
        "--no-write",
        action="store_true",
        help="Generate ordered structures, but do not store them into disk.",
    )
    group.add_argument(
        "--dir-size",
        type=int,
        default=const.DEFAULT_DIR_SIZE,
        help="Number of output CIFs written to each output directories.",
    )
    group.add_argument(
        "--write-symm",
        action="store_true",
        help="Write symmetries for all output CIFs (slower).",
    )
    group.add_argument(
        "--write-ewald",
        action="store_true",
        help="Write Ewald energy for all output CIFs.",
    )
    group.add_argument(
        "--symmetrize",
        action="store_true",
        help=(
            "Use Wyckoff labels from a symmetry search to label the input CIF's sites. "
            "By default, the label given within the CIF is used."
        ),
    )
    group.add_argument(
        "--max-ewald",
        type=float,
        help="Set maximum Ewald energy for the output CIFs",
    )

    group = parser.add_argument_group("run configuration")
    group.add_argument(
        "--symprec",
        type=float,
        default=const.DEFAULT_SYMPREC,
        help="Symmetry search precision (simulation cell fraction).",
    )
    group.add_argument(
        "--atol",
        type=float,
        default=const.DEFAULT_ATOL,
        help="Discretization absolute tolerance (angstrom).",
    )
    group.add_argument(
        "--angle-tolerance",
        type=float,
        default=const.DEFAULT_ANGLE_TOLERANCE,
        help="Symmetry search angle tolerance (degrees).",
    )
    group.add_argument(
        "--sample",
        type=int,
        default=const.DEFAULT_SAMPLE,
        help="Sample N number of ordered structures.",
    )
    group.add_argument(
        "--seed",
        type=int,
        default=const.DEFAULT_SEED,
        help="Random seed for the random number generator."
        " Used only when --sample is set."
        " Defaults to the (integer rounded) unixtime.",
    )
    group.add_argument(
        "--disable-progressbar",
        action="store_true",
        help="Disable progress bar. May stabilize longer runs.",
    )
    group.add_argument(
        "--no-dmat",
        action="store_true",
        # help="(devel/algo) Alternative algorithm without distance matrix (slower).",
        help=argparse.SUPPRESS,
    )
    group.add_argument(
        "--t-kind",
        default="sum",
        choices=("sum", "plsum", "det"),
        # help=(
        #     "(devel/algo) Type of T function applied to "
        #     "distance matrix (sum, plsum, det)."
        # ),
        help=argparse.SUPPRESS,
    )
    return parser


def _main():  # pylint: disable=missing-function-docstring
    parser = _get_parser()
    args = parser.parse_args()
    const.DISABLE_PROGRESSBAR = args.disable_progressbar
    logging_state = _setup_tqdm_logging()
    try:
        # Print header first for faster perceived response
        _print_header()

        # Late patch for count/mod/nowrite
        if args.count_only:
            args.no_write = True

        from .main import _ScriptHelper  # pylint:disable=import-outside-toplevel

        # Trick to allow ",", ";", and whiteline as separator
        from_species = const.FLEXIBLE_SEPARATOR.split(",".join(args.from_species))
        to_species = const.FLEXIBLE_SEPARATOR.split(",".join(args.to_species))
        scaling_matrix = [int(x) for x in const.FLEXIBLE_SEPARATOR.split(",".join(args.scaling_matrix))]

        # check the dimension of the scaling matrix
        if not len(scaling_matrix) in {1, 3, 9}:
            logging.warning("The scaling_matrix should be 1, 3, or 9 scalar values.")
            raise ValueError
        else:
            if len(scaling_matrix) == 9:
                scaling_matrix = np.array(scaling_matrix).reshape(3, 3)
            else:
                scaling_matrix = np.array(scaling_matrix)

        from_species = list(filter(None, from_species))
        to_species = list(filter(None, to_species))
        # scaling_matrix = list(filter(None, scaling_matrix)) #here! the problem is that filter removes 0!

        helper = _ScriptHelper(
            structure_file=args.input,
            from_species=from_species,
            to_species=to_species,
            scaling_matrix=scaling_matrix,
            symmetrize=args.symmetrize,
            sample=args.sample,
            seed=args.seed,
            symprec=args.symprec,
            atol=args.atol,
            angle_tolerance=args.angle_tolerance,
            dir_size=args.dir_size,
            write_symm=args.write_symm,
            write_ewald=args.write_ewald,
            max_ewald=args.max_ewald,
            no_write=args.no_write,
            n_jobs=-1,
            no_dmat=args.no_dmat,
            t_kind=args.t_kind,
        )
        helper.count()
        helper.save_modified_structure()
        if not args.count_only and not args.mod_only:
            helper.write()
        _print_footer()
    finally:
        _restore_logging(logging_state)


if __name__ == "__main__":
    _main()
