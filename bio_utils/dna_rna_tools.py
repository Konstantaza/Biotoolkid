# bio_utils/dna_rna_tools.py
from typing import Callable, List, Tuple, Union, Sequence


def is_nucleic_acid(seq: str) -> bool:
    """
    Проверяет, является ли строка нуклеиновой кислотой (ДНК или РНК).

    Args:
        seq (str): Входная последовательность.

    Returns:
        bool: True, если последовательность валидна, иначе False.
    """
    if not isinstance(seq, str) or not seq:
        return False
    # Проверка через множества
    valid_nucleotides = set("ATGCU")
    seq_upper = seq.upper()
    return set(seq_upper).issubset(valid_nucleotides) and not ("T" in seq_upper and "U" in seq_upper)


def transcribe_sequence(seq: str) -> str:
    """Транскрибирует ДНК последовательность в РНК."""
    return seq.upper().replace("T", "U")


def reverse_sequence(seq: str) -> str:
    """Возвращает перевернутую последовательность."""
    return seq.upper()[::-1]


def complement_sequence(seq: str) -> str:
    """Возвращает комплементарную последовательность."""
    seq_upper = seq.upper()
    complement_map = str.maketrans("ATGCU", "TACGA")
    return seq_upper.translate(complement_map)


def reverse_complement_sequence(seq: str) -> str:
    """Возвращает обратную комплементарную последовательность."""
    return complement_sequence(reverse_sequence(seq))


def get_command_handler(command: str) -> Tuple[str, Callable[[str], Union[str, bool]]]:
    """
    Возвращает тип команды и соответствующую ей функцию-обработчик.

    Args:
        command (str): Название команды.

    Returns:
        Tuple[str, Callable]: Кортеж из типа команды ("check" или "action") и функции.
    """
    if command == "is_nucleic_acid":
        return "check", is_nucleic_acid

    commands = {
        "transcribe": transcribe_sequence,
        "reverse": reverse_sequence,
        "complement": complement_sequence,
        "reverse_complement": reverse_complement_sequence,
    }

    if command in commands:
        return "action", commands[command]

    raise ValueError(f"Неизвестная команда: '{command}'")


def apply_actions(
    sequences: Sequence[str], mode: str, func: Callable
) -> Union[List[Union[str, bool, None]], Union[str, bool, None]]:
    """
    Применяет функцию к списку последовательностей в зависимости от режима.

    Args:
        sequences (Sequence[str]): Кортеж или список последовательностей.
        mode (str): Режим работы ("check" или "action").
        func (Callable): Функция для применения.

    Returns:
        Union[List, str, bool, None]: Список или одиночный результат.
    """
    results = []
    for seq in sequences:
        if mode == "check":
            results.append(func(seq))
        elif mode == "action":
            if not is_nucleic_acid(seq):
                results.append(None)
            else:
                results.append(func(seq))

    return results[0] if len(results) == 1 else results
