import argparse, os
from typing import List, Optional
from .core import run_MixMHCpred

def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse and validate command-line arguments.

    Validates required arguments and file paths. Raises SystemExit if arguments are invalid.
    """
    parser = argparse.ArgumentParser(
        description='MixMHCpred: Predict MHC-I ligands across alleles and species.',
        prog='MixMHCpred',
        add_help=True
    )
    parser.add_argument(
        "-i", "--InputFilePath",
        required=True,
        help="file path to peptide input file (required)"
    )
    parser.add_argument(
        "-o", "--OutputFilePath",
        required=True,
        help="file path where to save results (required)"
    )
    parser.add_argument(
        "-a", "--Alleles",
        required=True,
        help="Alleles to compute (comma-separated), e.g. A0101,B0702,C0107 (required)"
    )
    parser.add_argument(
        "-l", "--lib_path",
        required=True,
        help="Library directory path (required)"
    )
    parser.add_argument(
        "-m", "--output_motifs",
        default='0',
        choices=['0', '1'],
        help="Create binding motifs (default: 0)"
    )
    parser.add_argument(
        "-p", "--score_peptides",
        default='1',
        choices=['0', '1'],
        help="Compute binding scores (default: 1)"
    )

    args = parser.parse_args(argv)

    # Additional validation
    if not os.path.isfile(args.InputFilePath):
        parser.error(f"Input file does not exist: {args.InputFilePath}")

    if not os.path.isdir(args.lib_path):
        parser.error(f"Library directory does not exist: {args.lib_path}")

    if os.path.isfile(args.OutputFilePath):
        parser.error(f"Output file '{args.OutputFilePath}' already exists. Please choose a different path.")

    if not args.Alleles or not args.Alleles.strip():
        parser.error("Alleles argument cannot be empty.")

    try:
        _ = int(args.output_motifs)
        if int(args.output_motifs) not in (0, 1):
            raise ValueError
    except ValueError:
        parser.error("output_motifs (-m) must be '0' or '1'")

    try:
        _ = int(args.score_peptides)
        if int(args.score_peptides) not in (0, 1):
            raise ValueError
    except ValueError:
        parser.error("score_peptides (-p) must be '0' or '1'")

    return args


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)
    
    file_input = args.InputFilePath
    file_output = args.OutputFilePath
    lib_path = args.lib_path
    alleles_raw = args.Alleles or ''

    header_comments, ligands = run_MixMHCpred(file_input, lib_path, alleles_raw)
    
    with open(file_output, 'w') as f:
        f.write('\n'.join(header_comments) + '\n')
        ligands.to_csv(f, sep='\t', index=False, lineterminator='\n')

    print(f"\n{file_output} is created\n")
    print('\nDONE\n')

if __name__ == "__main__":
    main()