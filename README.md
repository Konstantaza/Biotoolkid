# Biotoolkid

# BioTooLKid 🧬

A small Python package for basic bioinformatics tasks, including FASTQ filtering and DNA/RNA sequence processing. This project was created as part of a Python programming training course.

---

## ## Main features

-   **FASTQ Filtering**:
    -   Filter by GC content bounds.
    -   Filter by sequence length bounds.
    -   Filter by minimum average quality score (Phred33).
-   **DNA/RNA Tools**:
    -   `transcribe`: Transcribes DNA into RNA (`T` → `U`).
    -   `reverse`: Reverses a sequence.
    -   `complement`: Generates the complementary strand.
    -   `reverse_complement`: Generates the reverse complementary strand.

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

---

## ## Author

-   **Name**: Konstantin Yamshchikov
-   **GitHub**: [@Konstantaza](https://github.com/Konstantaza)
