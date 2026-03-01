# BioTooLKid 🧬

A Python-based bioinformatics toolkit redesigned with OOP principles and Biopython. This project provides an intuitive interface for manipulating biological sequences and filtering FASTQ files.

---

## ## Features

### ### Biological Sequence Classes (OOP)
A robust hierarchy of classes for working with biological sequences:
* `BiologicalSequence`: The abstract base class providing common interfaces (length, slicing, alphabet validation).
* `NucleicAcidSequence`: Intermediate class with polymorphic methods for `complement`, `reverse`, and `reverse_complement`.
* `DNASequence` & `RNASequence`: Concrete classes with specific alphabets and translation maps. DNA includes a `transcribe()` method returning an RNA object.
* `AminoAcidSequence`: A class for protein sequences, featuring custom methods like `count_hydrophobic()`.

### ### FASTQ Filtering (Biopython)
The `filter_fastq` function utilizes `Bio.SeqIO` and `Bio.SeqUtils` to efficiently parse and filter FASTQ files based on:
* **Length bounds**
* **GC-content bounds**
* **Average Phred quality (Q-score)**

### ### (Addition) Bio-Files Processor (`bio_files_processor.py`)

-   **FASTA Converter**: Converts multi-line FASTA files into a single-line format for easier parsing.
-   **BLAST Parser**: Extracts and sorts the top hit descriptions for each query from a BLAST text output file.

---

## ## Installation

To use it, please clone the repository and install the required dependencies:

```bash
git clone [https://github.com/Konstantaza/Biotoolkid.git](https://github.com/Konstantaza/Biotoolkid.git)
cd Biotoolkid
pip install -r requirements.txt
```

---

## ## Usage Example

The main entry point is main.py. It contains examples of how to use the implemented functions on sample data.

To run the demonstration, execute:

```bash
python3 main.py

Test OOP Biological Sequences
Original DNA: DNASequence: ATGC
DNA length: 4
Is alphabet valid? True
Complement: DNASequence: TACG
Reverse complement: DNASequence: GCAT
Transcribed to RNA: RNASequence: AUGC
RNA complement: RNASequence: UACG
Sliced DNA [1:3]: DNASequence: TG (Type: DNASequence)
Protein: AminoAcidSequence: MAVW
Hydrophobic count: 4
```


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
