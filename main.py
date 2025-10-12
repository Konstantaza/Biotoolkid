# main.py
from typing import Dict, List, Tuple, Union

# Импортируем вспомогательные функции из модулей
from bio_utils.fastq_tools import (
    is_passing_gc_filter,
    is_passing_length_filter,
    is_passing_quality_filter,
)
from bio_utils.dna_rna_tools import (
    apply_actions,
    get_command_handler
)

# Определяем главную функцию
def filter_fastq(
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Union[int, float, Tuple[Union[int, float], Union[int, float]]] = (
        0, 100
    ),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: int = 0,
) -> Dict[str, Tuple[str, str]]:
    """
    Фильтрует FASTQ последовательности по GC-составу, длине и качеству.
    Эта функция является точкой входа и делегирует задачи
    вспомогательным функциям.
    """
    filtered_seqs = {}
    for name, (seq, quality) in seqs.items():
        if (
            is_passing_gc_filter(seq, gc_bounds)
            and is_passing_length_filter(seq, length_bounds)
            and is_passing_quality_filter(quality, quality_threshold)
        ):
            filtered_seqs[name] = (seq, quality)
    return filtered_seqs

# Определяем вторую главную функцию
def run_dna_rna_tools(
    *sequences: str, command: str
) -> Union[List[Union[str, bool, None]], Union[str, bool, None]]:
    """
    Главная функция для работы с ДНК/РНК последовательностями.
    Эта функция является точкой входа и делегирует задачи
    вспомогательным функциям.
    """
    if not sequences:
        raise ValueError("Требуется как минимум одна последовательность")

    mode, handler = get_command_handler(command)
    return apply_actions(sequences, mode, handler)


# Блок для демонстрации работы.
# В реальном проекте здесь может быть парсинг аргументов командной строки.
if __name__ == "__main__":
    from bio_utils.example_data import FASTQ_SEQS

    # Пример вызова filter_fastq
    filtered_results = filter_fastq(FASTQ_SEQS, gc_bounds=(40, 60))
    print("Отфильтрованные FASTQ:", filtered_results)

    # Пример вызова run_dna_rna_tools
    transcribed_results = run_dna_rna_tools(
        "GATTACA", "AGCTTGAAACGT", command="transcribe"
    )
    print("Транскрибированные последовательности:", transcribed_results)
