import argparse
from dataclasses import replace
from pathlib import Path
from tqdm import tqdm
from typing import Callable, Dict

def define_arguments(parser: argparse.ArgumentParser):
    subparsers = parser.add_subparsers(required=True, dest="command")

    import_args = subparsers.add_parser(
        "import",
        help="Convert a taxonomy TSV file to a taxonomy TSV DB.")
    import_args.add_argument(
        "--depth",
        type=int,
        default=None)
    import_args.add_argument(
        "--fasta-db",
        type=Path,
        required=True)
    import_args.add_argument(
        "input_path",
        type=Path,
        help="The taxonomy TSV file to convert.")
    import_args.add_argument(
        "output_path",
        type=Path,
        nargs='?',
        default=None,
        help="The destination taxonomy TSV DB file.")

    export_args = subparsers.add_parser(
        "export",
        help="Convert a taxonomy TSV DB to a taxonomy TSV file.")
    export_args.add_argument(
        "input_path",
        type=Path,
        help="The taxonomy TSV DB file to convert.")
    export_args.add_argument(
        "output_path",
        type=Path,
        nargs='?',
        default=None,
        help="The destination taxonomy TSV file.")

    info_args = subparsers.add_parser(
        "info",
        help="Display information about a taxonomy TSV or TSV DB file.")
    info_args.add_argument(
        "input_path",
        type=Path,
        help="The taxonomy TSV or TSV DB file to display information about.")

    lookup_args = subparsers.add_parser(
        "lookup",
        help="Lookup A FASTA ID in the given taxonomy TSV or TSV DB.")
    lookup_args.add_argument(
        "input_path",
        type=Path,
        help="The taxonomy TSV or TSV DB file to lookup the FASTA ID in.")
    lookup_args.add_argument(
        "ids",
        nargs='+',
        type=str,
        help="The FASTA ID(s) to lookup.")


def command_import(config: argparse.Namespace):
    print("Importing taxonomy TSV...")
    from dnadb import fasta, taxonomy
    output_path = config.output_path
    if output_path is None:
        output_path = config.input_path.with_suffix(".tsv.db")
    fasta_db = fasta.FastaDb(config.fasta_db)
    num_skipped = 0
    num_processed = 0
    with taxonomy.TaxonomyDbFactory(output_path, fasta.FastaDb(config.fasta_db), depth=config.depth) as factory:
        for entry in tqdm(taxonomy.entries(config.input_path)):
            if entry.sequence_id not in fasta_db:
                num_skipped += 1
                continue
            if config.depth is not None:
                taxons = taxonomy.split_taxonomy(entry.label, keep_empty=True)[:config.depth]
                entry = replace(entry, label=taxonomy.join_taxonomy(taxons, depth=config.depth))
            factory.write_entry(entry)
            num_processed += 1
    print(f"Done. Imported {num_processed:,} sequences. Skipped {num_skipped:,} sequences.")


def command_export(config: argparse.Namespace):
    print("Exporting taxonomy TSV DB...")
    from dnadb import taxonomy
    output_path = config.output_path
    if output_path is None:
        output_path = config.input_path.with_suffix("") # Remove .db suffix
    with open(output_path, "w") as output, taxonomy.TaxonomyDb(config.input_path) as db:
        output.write("Feature ID\tTaxon\n")
        taxonomy.write(output, tqdm(db))


def command_info(config: argparse.Namespace):
    from dnadb import taxonomy
    uuid = None
    if config.input_path.suffix == ".db":
        db = taxonomy.TaxonomyDb(config.input_path)
        uuid = db.uuid
        count = sum(db.count(label_index) for label_index in range(db.num_labels))
        unique_labels = len(db)
    else:
        entries = taxonomy.entries(config.input_path)
        count = 0
        unique_labels = set()
        for entry in entries:
            count += 1
            unique_labels.add(entry.label)
        unique_labels = len(unique_labels)

    print(f"Info for: {config.input_path}")
    if uuid is not None:
        print(f"               UUID: {uuid}")
    print(f"             Length: {count:,}")
    print(f"  Unique Taxonomies: {unique_labels:,}")


def command_lookup(config: argparse.Namespace):
    from dnadb import taxonomy
    if config.input_path.suffix == ".db":
        db = taxonomy.TaxonomyDb(config.input_path)
        for id in config.ids:
            if id not in db:
                print(f"'>{id}' not found.")
            else:
                print(f"{id}\t{db[id]}")
    else:
        entries = taxonomy.entries(config.input_path)
        entries_to_print: Dict[str, str|None] = {id: None for id in config.ids}
        found = 0
        for entry in entries:
            if entry.sequence_id in entries_to_print:
                found += 1
                entries_to_print[entry.sequence_id] = entry.label
                if found == len(entries_to_print):
                    break
        for id, entry in entries_to_print.items():
            if entry is None:
                print(f"'{id}' not found.")
            else:
                print(entry)
