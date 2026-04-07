#BLAST101 CLI
# Code adapted from lecture material and class practical unless otherwise stated in comments

import os
import argparse
import sys
import unittest

import programme_settings
import blast_101_search
import smith_waterman_search
import test_blast101

#Validate query sequence
def validate_query_sequence(seq: str, valid_residues: str) -> str:
    if seq is None:
        raise ValueError("No query sequence provided")

    seq = seq.strip().upper()

    if len(seq) == 0:
        raise ValueError("Sequence cannot be empty")

    valid_set = set(valid_residues.upper())
    invalid_chars = sorted(set(seq) - valid_set)

    if invalid_chars:
        raise ValueError(
            f"Sequence contains invalid character(s): {', '.join(invalid_chars)}"
        )

    # reject likely DNA input
    dna_set = set("ACGT")
    if set(seq).issubset(dna_set):
        raise ValueError("Sequence appears to be DNA, not protein")

    return seq
#Validate database file
def validate_database_file(path: str) -> str:
    if path is None:
        raise ValueError("No database path provided")
    if len(path) == 0:
        raise ValueError("Database file cannot be empty")
    if not os.path.isfile(path):
        raise ValueError(f"Database file cannot be found: {path}")
    if os.path.getsize(path) == 0:
        raise ValueError(f"Database file cannot be empty: {path}")

    nonempty = None
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                nonempty = line
                break
    if nonempty is None:
        raise ValueError("Database file contains no readable content")
    if not nonempty.startswith(">"):
        raise ValueError("Database file does not appear to be in FASTA format")
    return path
#Print state of inputs before running
def print_run_summary(mode: str, query_seq: str | None, database: str | None) -> None:
    print(f"=========================================================")
    print("BLAST101 Command Line Interface")
    print("==========================================================")
    print(f"Mode: {mode}")
    if query_seq is not None:
        print(f"Query sequence: {query_seq}")
    if database is not None:
        print(f"Database path: {database}")
    print("Using remaining parameters from settings.ini")
    print("=========================================================")
    if mode == "blast":
        print("Running BLAST-like local search...")
    elif mode == "sw":
        print("Running Smith-Waterman local alignment search...")
    elif mode == "test":
        print("Running automated unit tests...")
    print("=========================================================")
#Run tests
def run_tests() -> None:
    import test_blast101
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromModule(test_blast101)
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
    if not result.wasSuccessful():
        sys.exit(1)

# argparse module for command line - https://docs.python.org/3/library/argparse.html
def main():
    parser = argparse.ArgumentParser(description="BLAST101 Command Line Interface")
    parser.add_argument("--mode", required=True, choices=["blast", "sw", "test"],
                        help="Choose between blast, sw or test")
    parser.add_argument("--database", help="FASTA database file to search")
    parser.add_argument("--query-seq", help="Protein query sequence as a raw amino acid string")

    args = parser.parse_args()

    #loads settings first to read valid residues and defaults
    programme_settings.read()

    if args.mode == "test":
        print_run_summary("test", query_seq=None, database=None)
        run_tests()
        return

    if args.query_seq is None:
        parser.error("--query-seq is required for blast and sw modes")

    if args.database is None:
        parser.error("--database is required for blast and sw modes")

    valid_residues = programme_settings.settings["BLAST"]["valid_residues"]

    try:
        query_seq = validate_query_sequence(args.query_seq, valid_residues)
        database = validate_database_file(args.database)
    except ValueError as exc:
        print(f'Input error: {exc}')
        sys.exit(1)

    #overriding selected settings in memory only
    programme_settings.settings["DEFAULT"]["query_sequence"] = query_seq
    programme_settings.settings["DEFAULT"]["database"] = database

    print_run_summary(mode=args.mode, query_seq=query_seq, database=database)

    if args.mode == "blast":
        blast_101_search.blast101_run()
    elif args.mode == "sw":
        smith_waterman_search.smith_waterman_run(read_settings=False)

if __name__ == "__main__":
    main()

