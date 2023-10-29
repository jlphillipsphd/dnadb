import argparse
from pathlib import Path
from tqdm import tqdm
from typing import Dict

def define_arguments(parser: argparse.ArgumentParser):
    subparsers = parser.add_subparsers(required=True, dest="command")

    import_args = subparsers.add_parser(
        "import",
        help="Convert a taxonomy TSV file to a taxonomy TSV DB.")
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
    from dnadb import taxonomy
    output_path = config.output_path
    if output_path is None:
        output_path = config.input_path.with_suffix(".tsv.db")
    with taxonomy.TaxonomyDbFactory(output_path) as factory:
        factory.write_entries(tqdm(taxonomy.entries(config.input_path)))


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
    if config.input_path.suffix == ".db":
        db = taxonomy.TaxonomyDb(config.input_path)
        count = sum(db.count(label_index) for label_index in range(len(db)))
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
    print(f"              Length: {count:,}")
    print(f"  Unique Taxonomies: {unique_labels:,}")


def command_lookup(config: argparse.Namespace):
    from dnadb import taxonomy
    if config.input_path.suffix == ".db":
        db = taxonomy.TaxonomyDb(config.input_path)
        for id in config.ids:
            if not db.contains_fasta_id(id):
                print(f"'>{id}' not found.")
            else:
                print(f"{id}\t{db.fasta_id_to_label(id)}")
    else:
        entries = taxonomy.entries(config.input_path)
        entries_to_print: Dict[str, str|None] = {id: None for id in config.ids}
        found = 0
        for entry in entries:
            if entry.identifier in entries_to_print:
                found += 1
                entries_to_print[entry.identifier] = entry.label
                if found == len(entries_to_print):
                    break
        for id, entry in entries_to_print.items():
            if entry is None:
                print(f"'{id}' not found.")
            else:
                print(entry)
