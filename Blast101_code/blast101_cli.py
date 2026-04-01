
# 1. run BLAST-like local search using existing blast_101_search.py
# 2. run smith-waterman search using existing smith_waterman_search.py
# 3. run tests using a test file you create
# 4. print settings/ execution summary before running

# this is the controller for existing code

#PARSE ARGUMENTS
#VALIDATE INPUTS
#LOAD SETTINGS
#OVERRIDE SELECTED SETTINGS IN MEMORY
#CALL EITHER BLAST_101_SEARCH.BLAST101_RUN() OR SMITH_WATERMAN_SEARCH.. OR TESTS

#python blast101_cli.py --mode blast --query-seq PWNAAPLHNFGEDFLQPYVQLQQNFSASDLEVNLEATRESHAHFSTPQALELFLNYSVTP --database uniprot_bit2.fasta

import os
import argparse
import sys
import unittest

import programme_settings
import blast_101_search
import smith_waterman_search

def validate_query_sequence(seq: str, valid_residues: str) -> str:
    if len(seq) == 0:
        raise ValueError("Sequence cannot be empty")

    valid_set = set(valid_residues.upper())
    invalid_chars = sorted(set(seq)- valid_set)
    if invalid_chars:
        raise ValueError("Sequence contains invalid characters")

    if seq is None:
        raise ValueError("Sequence cannot be empty")
    seq = seq.strip().upper()

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

def print_run_summary(mode: str, query_seq: str | None, database_path: str) -> None:
    print(f"=========================================================")
    print("BLAST101 Command Line Interface")
    print("==========================================================")
    print(f"Mode: {mode}")
    if query_seq is not None:
        print(f"Query: {query_seq}")
    if database_path is not None:
        print(f"Database path: {database_path}")
    print("Using remaining parameters from settings.ini")
    print("=========================================================")

def run_tests() -> None:
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
    if not(result.wasSuccessful()):
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="BLAST101 Command Line Interface")
    parser.add_argument("--mode", required=True, choices="blast", "sw", "test",
                        help="BLAST101 Command Line Interface")
    parser.add_argument("--database", help="FASTA database file to search")
    parser.add_argument("--query", help="FASTA query file to search")

    args = parser.parse_args()
    programme_settings.read(args.database)
    programme_settings.read(args.database)
