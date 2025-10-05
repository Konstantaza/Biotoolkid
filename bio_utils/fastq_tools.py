# bio_utils/fastq_tools.py
from typing import Dict, Tuple, Union

# Объявляем псевдоним типа для удобства.
FastqDict = Dict[str, Tuple[str, str]]


def passes_gc_filter(
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
    return gc_bounds[0] <= gc_content <= gc_bounds[1]


def passes_length_filter(
    seq: str, length_bounds: Union[int, Tuple[int, int]]
) -> bool:
    """Проверяет, проходит ли последовательность по длине."""
    seq_len = len(seq)
    if isinstance(length_bounds, int):
        return seq_len <= length_bounds
    return length_bounds[0] <= seq_len <= length_bounds[1]


def passes_quality_filter(quality_str: str, quality_threshold: int) -> bool:
    """Проверяет, проходит ли последовательность по среднему качеству (Phred33)."""
    if not quality_str:
        return False
    # Следующая строка намеренно длиннее 79 символов для читаемости
    total_quality = sum(ord(char) - 33 for char in quality_str)  # noqa: E501
    average_quality = total_quality / len(quality_str)
    return average_quality >= quality_threshold


def filter_fastq(
    seqs: FastqDict,
    gc_bounds: Union[int, float, Tuple[Union[int, float], Union[int, float]]] = (
        0,
        100,
    ),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: int = 0,
) -> FastqDict:
    """
    Фильтрует FASTQ последовательности по GC-составу, длине и качеству.

    Args:
        seqs (FastqDict): Словарь с FASTQ данными.
        gc_bounds (Union[int, float, Tuple]): Границы GC-состава (в %).
            Можно передать одно число (верхняя граница) или кортеж
            (нижняя, верхняя).
        length_bounds (Union[int, Tuple]): Границы длины. Аналогично
            gc_bounds.
        quality_threshold (int): Пороговое значение среднего качества
            (Phred33).

    Returns:
        FastqDict: Отфильтрованный словарь с FASTQ данными.
    """
    filtered_seqs = {}
    for name, (seq, quality) in seqs.items():
        if (
            passes_gc_filter(seq, gc_bounds)
            and passes_length_filter(seq, length_bounds)
            and passes_quality_filter(quality, quality_threshold)
        ):
            filtered_seqs[name] = (seq, quality)
    return filtered_seqs
