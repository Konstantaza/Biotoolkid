import re
from pathlib import Path
from typing import Optional, Generator, Tuple

# Function 1: FASTA Processing
def _fasta_parser(input_fasta: str) -> Generator[Tuple[str, str], None, None]:
    """
    Helper generator to read FASTA files.
    Yields one record (header, sequence) at a time.
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
    Converts a multiline FASTA file to a oneline FASTA file.
    Saves the result in the 'processing_results' directory.
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
    print(f"Oneline FASTA saved to: {final_output_path}")


# Function 2: BLAST Processing
def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Extracts descriptions of the top hits from BLAST text output.
    Sorts the results alphabetically and saves to a file.
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
    print(f"BLAST parsing results saved to: {final_output_path}")


if __name__ == "__main__":
    print("--- Testing FASTA Converter ---")
    print("\n--- Testing BLAST Parser ---")