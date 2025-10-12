# bio_files_processor.py
import re
from pathlib import Path
from typing import Optional, Generator, Tuple

# Функция 1: Работа с FASTA


def _fasta_parser(input_fasta: str) -> Generator[Tuple[str, str], None, None]:
    """
    Вспомогательный генератор для чтения FASTA файлов.
    Выдает по одной записи (заголовок, последовательность) за раз.
    """
    header = None
    sequence = []
    with open(input_fasta, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    yield header, ''.join(sequence)
                header = line
                sequence = []
            else:
                sequence.append(line)
        if header:
            yield header, ''.join(sequence)


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: Optional[str] = None
) -> None:
    """
    Преобразует многострочный FASTA файл в однострочный.
    """
    output_dir = Path("processing_results")
    output_dir.mkdir(exist_ok=True)

    if output_fasta is None:
        base_name = Path(input_fasta).stem
        output_fasta = f"{base_name}_oneline.fasta"

    final_output_path = output_dir / output_fasta

    with open(final_output_path, 'w') as out_f:
        for header, sequence in _fasta_parser(input_fasta):
            out_f.write(f"{header}\n")
            out_f.write(f"{sequence}\n")
    print(f"Однострочный FASTA файл сохранен как: {final_output_path}")


# Функция 2: Работа с BLAST


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Извлекает описания лучших совпадений из текстового вывода BLAST.
    """
    output_dir = Path("processing_results")
    output_dir.mkdir(exist_ok=True)

    top_hits = []
    with open(input_file, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith('Sequences producing significant alignments:'):
            target_line_index = i + 3
            if target_line_index < len(lines):
                top_hit_line = lines[target_line_index].strip()
                parts = re.split(r'\s{2,}', top_hit_line)
                if parts:
                    top_hits.append(parts[0])

    top_hits.sort()
    final_output_path = output_dir / output_file

    with open(final_output_path, 'w') as out_f:
        for hit in top_hits:
            out_f.write(f"{hit}\n")
    print(f"Результаты парсинга BLAST сохранены в: {final_output_path}")


# Блок для демонстрации работы

if __name__ == "__main__":
    print("--- Тестирование конвертера FASTA ---")
    convert_multiline_fasta_to_oneline(
        "example_data/example_multiline_fasta.fasta"
    )

    print("\n--- Тестирование парсера BLAST ---")
    parse_blast_output(
        "example_data/example_blast_results.txt",
        output_file="blast_top_hits_sorted.txt"
    )
