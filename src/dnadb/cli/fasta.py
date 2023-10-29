import argparse
from pathlib import Path
from tqdm import tqdm
from typing import Dict

def define_arguments(parser: argparse.ArgumentParser):
    subparsers = parser.add_subparsers(required=True, dest="command")

    import_args = subparsers.add_parser(
        "import",
        help="Convert a FASTA file to a FASTA DB.")
    import_args.add_argument(
        "input_path",
        type=Path,
        help="The FASTA file to convert.")
    import_args.add_argument(
        "output_path",
        type=Path,
        nargs='?',
        default=None,
        help="The destination FASTA DB file.")

    export_args = subparsers.add_parser(
        "export",
        help="Convert a FASTA DB to a FASTA file.")
    export_args.add_argument(
        "input_path",
        type=Path,
        help="The FASTA DB file to convert.")
    export_args.add_argument(
        "output_path",
        type=Path,
        nargs='?',
        default=None,
        help="The destination FASTA file.")

    info_args = subparsers.add_parser(
        "info",
        help="Display information about a FASTA or FASTA DB file.")
    info_args.add_argument(
        "input_path",
        type=Path,
        help="The FASTA or FASTA DB file to display information about.")

    lookup_args = subparsers.add_parser(
        "lookup",
        help="Lookup A FASTA ID in the given FASTA or FASTA DB.")
    lookup_args.add_argument(
        "input_path",
        type=Path,
        help="The FASTA or FASTA DB file to lookup the FASTA ID in.")
    lookup_args.add_argument(
        "ids",
        nargs='+',
        type=str,
        help="The FASTA ID(s) to lookup.")


def command_import(config: argparse.Namespace):
    print("Importing FASTA...")
    from dnadb import fasta
    output_path = config.output_path
    if output_path is None:
        output_path = config.input_path.with_suffix(".fasta.db")
    with fasta.FastaDbFactory(output_path) as factory:
        factory.write_entries(tqdm(fasta.entries(config.input_path)))


def command_export(config: argparse.Namespace):
    print("Exporting FASTA DB...")
    from dnadb import fasta
    output_path = config.output_path
    if output_path is None:
        output_path = config.input_path.with_suffix("") # Remove .db suffix
    with open(output_path, "w") as output, fasta.FastaDb(config.input_path) as db:
        fasta.write(output, tqdm(db))


def command_info(config: argparse.Namespace):
    from dnadb import fasta
    if config.input_path.suffix == ".db":
        entries = iter(fasta.FastaDb(config.input_path))
    else:
        entries = fasta.entries(config.input_path)
    min_length = float("inf")
    max_length = float("-inf")
    count = 0
    for entry in entries:
        count += 1
        length = len(entry.sequence)
        min_length = min(min_length, length)
        max_length = max(max_length, length)
    print(f"Info for: {config.input_path}")
    print(f"             Length: {count:,}")
    print(f"  Shortest Sequence: {min_length:,}")
    print(f"   Longest Sequence: {max_length:,}")


def command_lookup(config: argparse.Namespace):
    from dnadb import fasta
    if config.input_path.suffix == ".db":
        db = fasta.FastaDb(config.input_path)
        for id in config.ids:
            if id not in db:
                print(f"'>{id}' not found.")
            else:
                print(db[id])
    else:
        entries = fasta.entries(config.input_path)
        entries_to_print: Dict[str, fasta.FastaEntry|None] = {id: None for id in config.ids}
        found = 0
        for entry in entries:
            if entry.identifier in entries_to_print:
                found += 1
                entries_to_print[entry.identifier] = entry
                if found == len(entries_to_print):
                    break
        for id, entry in entries_to_print.items():
            if entry is None:
                print(f"'{id}' not found.")
            else:
                print(entry)
