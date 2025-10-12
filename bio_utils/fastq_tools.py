# bio_utils/fastq_tools.py
from typing import Dict, Tuple, Union

FastqDict = Dict[str, Tuple[str, str]]


def is_passing_gc_filter(
    seq: str,
    gc_bounds: Union[int, float, Tuple[Union[int, float], Union[int, float]]],
) -> bool:
    """Проверяет, проходит ли последовательность по GC-составу."""
    if not seq:
        return False
    gc_count = seq.upper().count("G") + seq.upper().count("C")
    gc_content = (gc_count / len(seq)) * 100

    if isinstance(gc_bounds, (int, float)):
        return gc_content <= gc_bounds

    # Распаковка кортежа
    low_bound, upper_bound = gc_bounds
    return low_bound <= gc_content <= upper_bound


def is_passing_length_filter(
    seq: str, length_bounds: Union[int, Tuple[int, int]]
) -> bool:
    """Проверяет, проходит ли последовательность по длине."""
    seq_len = len(seq)
    if isinstance(length_bounds, int):
        return seq_len <= length_bounds
    low_bound, upper_bound = length_bounds
    return low_bound <= seq_len <= upper_bound


def is_passing_quality_filter(quality_str: str, quality_threshold: int) -> bool:
    """Проверяет, проходит ли последовательность по среднему качеству Phred33."""
    if not quality_str:
        return False
    quality_scores = [ord(char) - 33 for char in quality_str]
    average_quality = sum(quality_scores) / len(quality_str)
    return average_quality >= quality_threshold
