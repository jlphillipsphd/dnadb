import abc
from typing import Generator, Generic, TypeVar

from .dna import AbstractSequenceWrapper

T = TypeVar("T", bound=AbstractSequenceWrapper)
class ISample(abc.ABC, Generic[T]):

    name: str

    def __contains__(self, sequence_index: int) -> bool:
        """
        Check whether a sequence entry is in the sample.
        """
        raise NotImplementedError()

    def __getitem__(self, sequence_index: int) -> T:
        """
        Get a sequence entry by its index.
        """
        raise NotImplementedError()

    def __iter__(self) -> Generator[T, None, None]:
        """
        Iterate over all sequence entries in the sample.
        """
        raise NotImplementedError()

    def __len__(self) -> int:
        """
        Get the number of sequence entries in the sample.
        """
        raise NotImplementedError()
