# bio_utils/dna_rna_tools.py

from typing import Callable, List, Tuple, Union, Sequence


def is_nucleic_acid(seq: str) -> bool:
    """
    Проверяет, является ли строка
    нуклеиновой кислотой (ДНК или РНК),
    используя итеративный подход.

    Args:
        seq (str): Входная последовательность.

    Returns:
        bool: True, если последовательность валидна,
        иначе False.
    """
    if not isinstance(seq, str) or not seq:
        return False
    seq_upper = seq.upper()
    valid_nucleotides = {"A", "T", "G", "C", "U"}
    for ch in seq_upper:
        if ch not in valid_nucleotides:
            return False
    if "T" in seq_upper and "U" in seq_upper:
        return False
    return True


def transcribe_sequence(seq: str) -> str:
    """Транскрибирует ДНК последовательность в РНК."""
    if not isinstance(seq, str):
        raise TypeError("transcribe_sequence expects a string")
    return seq.replace("T", "U").replace("t", "u")


def reverse_sequence(seq):
    """Возвращает перевернутую последовательность."""
    if not isinstance(seq, str):
        raise TypeError("reverse_sequence expects a string")
    return seq[::-1]


def complement_sequence(seq):
    """
    Возвращает комплементарную последовательность, используя
    логику с проверкой на РНК и словарями.
    """
    if not isinstance(seq, str):
        raise TypeError("complement_sequence expects a string")
    is_rna = "u" in seq.lower()
    if is_rna:
        complement_map = {
            "A": "U",
            "U": "A",
            "G": "C",
            "C": "G",
            "a": "u",
            "u": "a",
            "g": "c",
            "c": "g",
        }
    else:
        complement_map = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "a": "t",
            "t": "a",
            "g": "c",
            "c": "g",
        }

    return "".join(complement_map.get(ch, ch) for ch in seq)


def reverse_complement_sequence(seq: str) -> str:
    """Возвращает обратную комплементарную последовательность."""
    if not isinstance(seq, str):
        raise TypeError(
              "reverse_complement_sequence expects a string"
              )
    return reverse_sequence(complement_sequence(seq))


def get_command_handler(command: str) -> Tuple[str, Callable[[str], Union[str, bool]]]:
    """
    Возвращает тип команды и
    соответствующую ей функцию-обработчик.

    Args:
        command (str): Название команды.

    Returns:
        Tuple[str, Callable]: Кортеж из
        типа команды ("check" или "action") и функции.

    Raises:
        ValueError: Если команда неизвестна.
    """
    if command == "is_nucleic_acid":
        return ("check", is_nucleic_acid)
    commands = {
        "transcribe": transcribe_sequence,
        "reverse": reverse_sequence,
        "complement": complement_sequence,
        "reverse_complement": reverse_complement_sequence,
    }
    if command not in commands:
        raise ValueError("Unknown command: {!r}".format(command))
    return ("action", commands[command])


def apply_actions(
    sequences: Sequence[str], mode: str, func: Callable
) -> Union[List[Union[str, bool, None]], Union[str, bool, None]]:
    """
    Применяет функцию к списку последовательностей
    в зависимости от режима.

    Args:
        sequences (Sequence[str]):
            Кортеж или список последовательностей.
        mode (str):
            Режим работы ("check" или "action").
        func (Callable):
            Функция для применения.

    Returns:
        Union[List, str, bool, None]:
            Список результатов или один результат,
            если на вход подана одна последовательность.
    """
    results = []
    if mode == "check":
        for sequence in sequences:
            results.append(func(sequence))
    else:
        for sequence in sequences:
            if not is_nucleic_acid(sequence):
                results.append(None)
            else:
                results.append(func(sequence))
    return results[0] if len(results) == 1 else results


def run_dna_rna_tools(
    *sequences: str, command: str
) -> Union[List[Union[str, bool, None]], Union[str, bool, None]]:
    """
    Главная функция для работы с
    ДНК/РНК последовательностями.

    Принимает любое количество
    последовательностей и одну команду.

    Args:
        *sequences (str):
            Одна или несколько
            последовательностей для обработки.
        command (str):
            Команда для выполнения
            (например, "transcribe").

    Returns:
        Результат(ы) выполнения команды.
        Может быть одним значением или списком.
    """
    if not sequences:
        raise ValueError("At least one sequence is required")

    mode, handler = get_command_handler(command)
    return apply_actions(sequences, mode, handler)
