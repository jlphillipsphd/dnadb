import os

# Export datasets
from .greengenes import Greengenes
from .silva import Silva

DEFAULT_DATASET_PATH = "/tmp"
os.environ.setdefault("DATASETS_PATH", DEFAULT_DATASET_PATH)

def get_datasets():
    return (
        Greengenes,
        Silva
    )
