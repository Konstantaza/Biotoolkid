# bio_utils/fastq_tools.py

from typing inport Dict, Tuple, Union

# Объявляем псевдоним типа для удобства.
# Теперь вместо конструкции Dict[str, Tuple[str, str]]
# можно писать FastqDict.

FastqDict = Dict[str, Tuple[str, str]]

def gc_filter(
    seq: str,
    gc_bounds: Union[int, float, Tuple[Union[int, float], Union[int, float]]]
) -> bool:
    """Функция проверяет, проходит ли последовательность по GC-составу."""
    # Избегаем деления на ноль для пустых последовательностей
    if not seq:
        return False
    gc_content = (seq.upper().count('G') + seq.upper().count('C')) / len(seq) * 100

    # Если передано одно число (верхняя граница)
    if isinstance(gc_bounds, (int, float)):
        return gc_content <= gc_bounds
    # Если передан кортеж (нижняя, верхняя граница)
    return gc_bounds[0] <= gc_content <= gc_bounds[1]


def length_filter(
    seq: str,
    length_bounds: Union[int, Tuple[int, int]]
) -> bool:
    """Функция проверяет, проходит ли последовательность по длине."""
    seq_len = len(seq)

    if isinstance(length_bounds, int):
        return seq_len <= length_bounds
    return length_bounds[0] <= seq_len <= length_bounds[1]


def quality_filter(
    quality_str: str,
    quality_threshold: int
) -> bool:
    """Функция проверяет, проходит ли последовательность по среднему качеству (Phred33)."""
    if not quality_str:
        return False
    # Суммируем ASCII-код каждого символа минус 33
    total_quality = sum(ord(char) - 33 for char in quality_str)
    average_quality = total_quality / len(quality_str)

    return average_quality >= quality_threshold


def filter_fastq(
    seqs: FastqDict,
    gc_bounds: Union[int, float, Tuple[Union[int, float], Union[int, float]]] = (0, 100),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: int = 0
) -> FastqDict:
    """
    Фильтрует FASTQ последовательности по GC-составу, длине и качеству.

    Args:
        seqs (FastqDict): Словарь с FASTQ данными.
        gc_bounds (Union[int, float, Tuple]): Границы GC-состава (в %).
                                              Можно передать одно число (верхняя граница)
                                              или кортеж (нижняя, верхняя).
        length_bounds (Union[int, Tuple]): Границы длины. Аналогично gc_bounds.
        quality_threshold (int): Пороговое значение среднего качества (Phred33).

    Returns:
        FastqDict: Отфильтрованный словарь с FASTQ данными.
    """
    filtered_seqs = {}
    for name, (seq, quality) in seqs.items():
        # Проверяем последовательность по всем трем условиям
        if (
            _passes_gc_filter(seq, gc_bounds) and
            _passes_length_filter(seq, length_bounds) and
            _passes_quality_filter(quality, quality_threshold)
        ):
            filtered_seqs[name] = (seq, quality)

    return filtered_seqs
