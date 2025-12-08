import sys
from argparse import Namespace
from .utils import init_logging_from_cut, timer
from .ModalAnalysis import ModalAnalysis, parse_arguments

def main(argobj=None):
    """
    Entry point for source_localization_core.

    - If argobj is a Namespace -> use it directly (programmatic use).
    - If argobj is a list/tuple of strings -> parse_arguments(argobj).
    - If argobj is None -> parse_arguments() from sys.argv[1:].
    """
    # Decide how to get `args`
    if isinstance(argobj, Namespace):
        args = argobj
    elif isinstance(argobj, (list, tuple)):
        args = parse_arguments(argobj)
    elif argobj is None:
        args = parse_arguments()
    else:
        raise TypeError(
            f"main() expected Namespace, list/tuple, or None, got {type(argobj)}"
        )

    if args is None:
        raise RuntimeError("parse_arguments() returned None – please fix its implementation.")

    # Initialize logging once, in serial
    init_logging_from_cut(args.compute_mode, args.var, args.dmd_output_dir)
    # Create and run ModalAnalysis
    modal_analysis = ModalAnalysis(args)
    # Pre-process the data
    modal_analysis.pre_process()
    if not args.extract_only:
        # Perform modal analysis
        modal_analysis.run_modal_analysis()
    

if __name__ == "__main__":
    main()