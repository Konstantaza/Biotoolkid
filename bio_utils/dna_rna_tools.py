def is_nucleic_acid(seq: str) -> bool:
    """
    Checks if the sequence contains only valid nucleotides (A, T, G, C, U)
    and does not mix DNA (T) and RNA (U) specific nucleotides.
    """
    if not isinstance(seq, str) or not seq:
        return False

    seq_upper = seq.upper()
    valid_nucleotides = {"A", "T", "G", "C", "U"}

    unique_chars = set(seq_upper)
    if not unique_chars.issubset(valid_nucleotides):
        return False

    if "T" in unique_chars and "U" in unique_chars:
        return False

    return True


def transcribe_sequence(seq: str) -> str:
    """Transcribes DNA sequence to RNA."""
    return seq.replace("T", "U").replace("t", "u")


def reverse_sequence(seq: str) -> str:
    """Reverses the sequence."""
    return seq[::-1]


def complement_sequence(seq: str) -> str:
    """Returns the complementary sequence."""
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
    """Returns the reverse complementary sequence."""
    return reverse_sequence(complement_sequence(seq))


def run_dna_rna_tools(*args):
    """
    Main entry point for DNA/RNA processing tools.
    Accepts sequences and a command as the last argument.
    """
    if len(args) < 2:
        raise ValueError("At least one sequence and one command are required")

    sequences = args[:-1]
    command = args[-1]

    valid_commands = {
        "transcribe",
        "reverse",
        "complement",
        "reverse_complement",
        "is_nucleic_acid",
    }
    if command not in valid_commands:
        raise ValueError(f"Unknown command: {command}")

    results = []

    for seq in sequences:
        if not is_nucleic_acid(seq):
            if command == "is_nucleic_acid":
                results.append(False)
            else:
                print(f"Warning: Sequence '{seq}' is invalid. Skipping.")
                results.append(seq)
        else:
            if command == "is_nucleic_acid":
                results.append(True)
            elif command == "transcribe":
                results.append(transcribe_sequence(seq))
            elif command == "reverse":
                results.append(reverse_sequence(seq))
            elif command == "complement":
                results.append(complement_sequence(seq))
            elif command == "reverse_complement":
                results.append(reverse_complement_sequence(seq))

    if len(results) == 1:
        return results[0]
    return results
