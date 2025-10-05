from bio_utils.dna_rna_tools import run_dna_rna_tools
from bio_utils.fastq_tools import filter_fastq
from bio_utils.example_data import FASTQ_SEQS


if __name__ == "__main__":
    # Тестируем фильтрацию FASTQ
    print("--->Testing FASTQ Filter")

    # Случай 1: без фильтров (должны вернуться почти все последовательности)
    print("\nOriginal data:", FASTQ_SEQS)
    filtered_default = filter_fastq(FASTQ_SEQS)
    print("Filtered (default):", filtered_default)

    # Случай 2: фильтруем по GC > 40-60%, длине 10-20 и качеству > 20
    filtered_custom = filter_fastq(
        seqs=FASTQ_SEQS,
        gc_bounds=(40, 60),
        length_bounds=(10, 20),
        quality_threshold=20,
    )
    print("\n--->Filtered (GC 40-60, len 10-20, qual>20):", filtered_custom)

    # Тестируем инструменты ДНК/РНК
    print("\n--->Testing DNA/RNA Tools")
    # Вызываем функцию с именованным аргументом command
    dna_results = run_dna_rna_tools("GATTACA",
                                    "AGCTTGAAACGT",
                                    command="transcribe")
    print(dna_results)
