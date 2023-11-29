import argparse
from pathlib import Path
import tempfile
from tqdm import tqdm
from typing import Dict

def define_arguments(parser: argparse.ArgumentParser):
    subparsers = parser.add_subparsers(required=True, dest="command")

    import_args = subparsers.add_parser(
        "import",
        help="Convert a FASTA file to a FASTA DB.")
    import_args.add_argument(
        "--min-length",
        type=int,
        default=0,
        help="The minimum length of a sequence to import.")
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

    import_multiplexed_args = subparsers.add_parser(
        "import-multiplexed",
        help="Import sequences from FASTQ files into a single FASTA DB with a corresonding sample mapping.")
    import_multiplexed_args.add_argument(
        "--min-length",
        type=int,
        default=0,
        help="The minimum length of a sequence to import.")
    import_multiplexed_args.add_argument(
        "--output-sequences-path",
        type=Path,
        required=True,
        help="The destination FASTA DB sequences file.")
    import_multiplexed_args.add_argument(
        "--output-mapping-path",
        type=Path,
        required=True,
        help="The destination FASTA mapping file.")
    import_multiplexed_args.add_argument(
        "input_sequences",
        type=Path,
        nargs="+",
        help="Sample files in either FASTA or FASTQ format.")

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
    num_skipped = 0
    num_processed = 0
    with fasta.FastaDbFactory(output_path) as factory:
        for entry in tqdm(fasta.entries(config.input_path)):
            if len(entry.sequence) < config.min_length:
                num_skipped += 1
                continue
            factory.write_entry(entry)
            num_processed += 1
    print(f"Done. Imported {num_processed:,} sequences. Skipped {num_skipped:,} sequences.")


def command_import_multiplexed(config: argparse.Namespace):
    print("Importing FASTA...")
    import gzip
    import heapq
    from itertools import count, repeat
    import re
    from dnadb import fasta

    num_skipped = 0
    num_sequences = 0
    num_processed = 0
    with tempfile.TemporaryDirectory() as scratch_path:
        scratch_path = Path(scratch_path)
        sequence_files: list[Path] = config.input_sequences
        print(f"Gathering Sequences from {len(sequence_files)} file(s)...")
        sequence_file_ids: dict[Path, Path] = {}
        for i, sequence_file in tqdm(enumerate(sorted(sequence_files))):
            scratch_file = scratch_path / f"sequences_{i}"
            sequence_file_ids[scratch_file] = sequence_file
            if sequence_file.name.endswith(".gz"):
                file_handle = gzip.open(sequence_file, "rt")
            else:
                file_handle = open(sequence_file)
            if sequence_file.name.endswith(".fastq") or sequence_file.name.endswith(".fastq.gz"):
                sequences = file_handle.readlines()[1::4]
            else:
                sequences = file_handle.readlines()[1::2]
            file_handle.close()
            num_processed += len(sequences)
            sequences = sorted(filter(lambda s: len(s) >= config.min_length, sequences))
            sequences.sort()
            num_sequences += len(sequences)
            with open(scratch_file, "w") as f:
                for sequence in sequences:
                    if len(sequence) < config.min_length:
                        num_skipped += 1
                        continue
                    f.write(sequence.rstrip() + '\n')

        # Create the output FASTA DB files
        print("Creating FASTA DB...")
        fasta_db = fasta.FastaDbFactory(config.output_sequences_path)
        prev: str = ""
        fasta_id_generator = count()
        for sequence in tqdm(heapq.merge(*map(open, sequence_file_ids.keys())), total=num_sequences):
            sequence = sequence.rstrip()
            if sequence == prev:
                continue
            identifier = str(next(fasta_id_generator))
            fasta_db.write_entry(fasta.FastaEntry(identifier, sequence))
            prev = sequence
        fasta_db.close()

        print("Creating FASTA Mapping DB...")
        fasta_db = fasta.FastaDb(config.output_sequences_path)
        fasta_mapping = fasta.FastaMappingDbFactory(config.output_mapping_path, fasta_db)
        mappings: dict[str, fasta.FastaMappingEntryFactory] = {}
        for fastq_file in sequence_file_ids.values():
            name = re.sub(r"(?:\.fasta|\.fastq)(?:\.gz)?$", "", fastq_file.name)
            mappings[name] = fasta_mapping.create_entry(name)

        print("Writing Sample Mappings...")
        fasta_entries = iter(fasta_db)
        fasta_entry = next(fasta_entries)
        for sequence, mapping in tqdm(heapq.merge(*[zip(f, repeat(mapping)) for f, mapping in zip(map(open, sequence_file_ids.keys()), mappings)]), total=num_sequences):
            sequence = sequence.rstrip()
            if sequence != fasta_entry.sequence:
                fasta_entry = next(fasta_entries)
            assert sequence == fasta_entry.sequence
            mappings[mapping].write_entry(fasta_entry)
        for mapping in mappings.values():
            fasta_mapping.write_entry(mapping)
        fasta_mapping.close()
    print(f"Done. Imported {num_sequences:,} sequences. Skipped {num_skipped:,} sequences.")


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
    uuid = None
    if config.input_path.suffix == ".db":
        fasta_db = fasta.FastaDb(config.input_path)
        uuid = fasta_db.uuid
        entries = iter(fasta_db)
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
    if uuid is not None:
        print(f"               UUID: {uuid}")
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
