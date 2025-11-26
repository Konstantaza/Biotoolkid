# main.py
from pathlib import Path
from typing import Tuple, Union, List

# Import FASTQ processing tools
from bio_utils.fastq_tools import (
    is_passing_gc_filter,
    is_passing_length_filter,
    is_passing_quality_filter,
    read_fastq,
    write_fastq,
)

# Import the main DNA/RNA processing function
from bio_utils.dna_rna_tools import run_dna_rna_tools


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[int, float, Tuple[Union[int, float], Union[int, float]]] = (0, 100),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """
    Filters FASTQ file by GC content, length, and quality, then saves the result.
    """
    input_path = Path(input_fastq)
    sequences = read_fastq(input_path)
    filtered_seqs = {}
    
    for name, (seq, quality) in sequences.items():
        if (
            is_passing_gc_filter(seq, gc_bounds)
            and is_passing_length_filter(seq, length_bounds)
            and is_passing_quality_filter(quality, quality_threshold)
        ):
            filtered_seqs[name] = (seq, quality)

    write_fastq(filtered_seqs, output_fastq)


if __name__ == "__main__":
    print("--- Test 1: DNA/RNA Tools ---")
    try:
        result = run_dna_rna_tools("ATG", "transcribe")
        print(f"ATG -> transcribe -> {result}")
    except Exception as e:
        print(f"Error testing DNA/RNA tools: {e}")

    print("\n--- Test 2: FASTQ Filtering ---")
    print("Main script finished successfully.")