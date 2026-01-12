from typing import Tuple

import pysam


class PysamTxEff:
    """
    Class used to query reference genome fasta file.    
    """

    def __init__(self, filename: str, size: int=2000) -> None:
        self.filename = filename
        self.my_fasta = pysam.FastaFile(self.filename)
        self.size = size

    def direct_query(self, chrom: str, start_pos: int, end_pos: int) -> str:
        """
        Provide reference sequence at given coordinates.
        :param chrom:
        :param start_pos:
        :param end_pos:
        :return: Sequence at these coordinates.
        """
        return self.my_fasta.fetch(reference=chrom, start=start_pos, end=end_pos)

    def faidx_query(self, chrom: str, pos: int) -> str:
        """
        Get the sequence corresponding to a given position and chromosome.  The coordinates in start and end need
        to be 0-based, as per the pysam docs.
        :param chrom:
        :param pos:
        :return:
        """
        start_pos, end_pos = self._get_endpts(chrom, pos)
        return self.my_fasta.fetch(reference=chrom, start=start_pos, end=end_pos)

    def _get_endpts(self, chrom: str, pos: int) -> Tuple[int, int]:
        """
        Find the search coordinates we will feed to pysam.FastaFile.fetch.
        :param chrom:
        :param pos:
        :return:
        """
        # Min and max coordinates for given reference chrom.
        max_pos: int = self.my_fasta.get_reference_length(chrom)
        min_pos: int = 1

        # Get the search range
        start_pos, end_pos = self._set_range(pos)
        start_pos = max(start_pos, min_pos) - 1
        end_pos = min(end_pos, max_pos) - 1

        return start_pos, end_pos

    def _set_range(self, pos: int) -> Tuple[int, int]:
        """
        Set the end point coordinates, based on a given size or range.
        :param pos:
        :return:
        """
        start = pos - self.size
        end = pos + self.size
        return start, end
