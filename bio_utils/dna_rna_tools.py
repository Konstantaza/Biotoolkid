def is_nucleic_acid(seq) -> bool:

    if not isinstance(seq, str):
        return False

    if seq == "":
        return False

    seq_upper = seq.upper()
    valid_nucleotides = {"A", "T", "G", "C", "U"}

    for ch in seq_upper:
        if ch not in valid_nucleotides:
            return False

    if "T" in seq_upper and "U" in seq_upper:
        return False

    return True


def transcribe_sequence(seq):

    if not isinstance(seq, str):
        raise TypeError("transcribe_sequence expects a string")

    return seq.replace("T", "U").replace("t", "u")


def reverse_sequence(seq):

    if not isinstance(seq, str):
        raise TypeError("reverse_sequence expects a string")

    return seq[::-1]


def complement_sequence(seq):

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


def reverse_complement_sequence(seq):

    if not isinstance(seq, str):
        raise TypeError("reverse_complement_sequence expects a string")

    return reverse_sequence(complement_sequence(seq))


def get_command_handler(command):

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


def apply_actions(sequences, mode, func):

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


def run_dna_rna_tools(*args):

    if len(args) < 2:
        raise ValueError("At least one sequence and one command are required")

    sequences = args[:-1]
    command = args[-1]

    mode, handler = get_command_handler(command)

    return apply_actions(sequences, mode, handler)
