# main.py
from pathlib import Path
from typing import Dict, List, Tuple, Union

# Импортируем вспомогательные функции, включая новые
from bio_utils.fastq_tools import (
    is_passing_gc_filter,
    is_passing_length_filter,
    is_passing_quality_filter,
    read_fastq,  # Функция читатель
    write_fastq, # Функция писатель
)
from bio_utils.dna_rna_tools import (
    apply_actions,
    get_command_handler
)

# Главные функции

def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[int, float, Tuple[Union[int, float], Union[int, float]]] = (0, 100),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """
    Фильтрует FASTQ файл по GC-составу, длине и качеству и сохраняет результат.
    """
    # Читаем исходный файл в словарь, используя функцию из модуля
    sequences = read_fastq(Path(input_fastq))
    filtered_seqs = {}
    for name, (seq, quality) in sequences.items():
        # Фильтруем, используя вспомогательные функции
        if (
            is_passing_gc_filter(seq, gc_bounds)
            and is_passing_length_filter(seq, length_bounds)
            and is_passing_quality_filter(quality, quality_threshold)
        ):
            filtered_seqs[name] = (seq, quality)

    # Записываем отфильтрованный словарь в новый файл
    write_fastq(filtered_seqs, output_fastq)


def run_dna_rna_tools(
    *sequences: str, command: str
) -> Union[List[Union[str, bool, None]], Union[str, bool, None]]:
    """Главная функция для работы с ДНК/РНК последовательностями."""
    if not sequences:
        raise ValueError("Требуется как минимум одна последовательность")
    mode, handler = get_command_handler(command)
    return apply_actions(sequences, mode, handler)


# Блок для демонстрации работы
if __name__ == "__main__":
    # Указываем путь к нашему тестовому файлу
    input_file = "example_data/example_fastq.fastq"
    # Указываем имя для выходного файла
    output_file = "filtered_reads.fastq"

    print(f"Запускаем фильтрацию файла: {input_file}")

    # Вызываем главную функцию с путями к файлам и параметрами фильтрации
    filter_fastq(
        input_fastq=input_file,
        output_fastq=output_file,
        gc_bounds=(30, 70),  # Например, ищем риды с GC от 30% до 70%
        length_bounds=(50, 150), # и длиной от 50 до 150
        quality_threshold=20 # и средним качеством > 20
    )

    print("Фильтрация завершена. Результат сохранен в папке 'filtered'")
