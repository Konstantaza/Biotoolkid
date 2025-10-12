# BioTooLKid 🧬

A Python project with a collection of bioinformatics utilities. It includes a file-based FASTQ filter and a separate script for processing various bioinformatics file formats like FASTA and BLAST outputs. This project was created as part of a Python programming course.

---

## ## Features

### ### Main Toolkit (`main.py`)

-   **FASTQ Filtering**: Reads a FASTQ file, filters sequences based on GC content, length, and quality, and saves the results to a new file in a `filtered` directory.
-   **DNA/RNA Tools**: A module with functions for basic sequence manipulation (transcription, reverse, complement, etc.).

### ### Bio-Files Processor (`bio_files_processor.py`)

-   **FASTA Converter**: Converts multi-line FASTA files into a single-line format for easier parsing.
-   **BLAST Parser**: Extracts and sorts the top hit descriptions for each query from a BLAST text output file.

---

## ## Installation

To use it, please clone the repository:

```bash
git clone [https://github.com/Konstantaza/Biotoolkid.git](https://github.com/Konstantaza/Biotoolkid.git)
cd Biotoolkid
```

---

## ## Usage Example

The main entry point is main.py. It contains examples of how to use the implemented functions on sample data.

To run the demonstration, execute:

```bash
python3 main.py
```

This command will run an analysis on the sample data from bio_utils/example_data.py and print the results to the console.

### Running the Bio-Files Processor
The bio_files_processor.py script demonstrates FASTA and BLAST file processing

```bash
python3 bio_files_processor.py
```
This command will process sample FASTA and BLAST files and save the results in the processing_results/ directory.
---

## ## Author

-   **Name**: Konstantin Yamshchikov
-   **GitHub**: [@Konstantaza](https://github.com/Konstantaza)
