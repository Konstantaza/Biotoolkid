# bio_utils/fastq_tools.py
from typing import Dict, Tuple, Union
from pathlib import Path
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


# Новая функция для чтения
def read_fastq(input_path: Path) -> FastqDict:
    """
    Читает FASTQ файл и преобразует его в словарь.
    """
    sequences = {}
    with open(input_path, 'r') as f:
        while True:
            name = f.readline().strip()
            if not name:
                break
            seq = f.readline().strip()
            f.readline()  # Пропускаем строку с '+'
            quality = f.readline().strip()
            sequences[name] = (seq, quality)
    return sequences

# Новая функция для записи
def write_fastq(sequences: FastqDict, output_path: str) -> None:
    """
    Записывает отфильтрованные последовательности в новый FASTQ файл.
    Создает папку 'filtered', если ее нет, и избегает перезаписи.
    """
    # Создаем путь к папке 'filtered' в текущей директории
    output_dir = Path.cwd() / "filtered"
    output_dir.mkdir(exist_ok=True)

    final_path = output_dir / output_path

    # Защита от случайной перезаписи
    if final_path.exists():
        print(f"Предупреждение: Файл {final_path} уже существует. Запись отменена.")
        return  # Прекращаем выполнение, чтобы не стереть данные

    with open(final_path, "w") as f:
        for name, (seq, quality) in sequences.items():
            f.write(f"{name}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{quality}\n")
