import argparse
import importlib

_module_names = [
    "fasta",
    "taxonomy"
]

def main():
    modules = {}

    parser = argparse.ArgumentParser("dnadb")
    subparsers = parser.add_subparsers(title="module", required=True, dest="module")
    for module_name in _module_names:
        modules[module_name] = importlib.import_module(f".{module_name}", "dnadb.cli")
        modules[module_name].define_arguments(subparsers.add_parser(module_name))
    config = parser.parse_args()
    getattr(modules[config.module], f"command_{config.command.replace('-', '_')}")(config)
