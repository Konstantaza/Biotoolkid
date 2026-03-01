# main.py
from pathlib import Path
from typing import Tuple, Union, List
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    def __init__(self, sequence: str):
        self.sequence = sequence


    @property
    @abstractmethod
    def alphabet(self) -> set:
        pass


    def check_alphabet(self) -> bool:
        """
        Проверяет алфавит. Работает для всех наследников благодаря self.alphabet.
        """
        return set(self.sequence.upper()).issubset(self.alphabet)
    

    def __len__(self) -> int:
        """Returns the length of the biological sequence."""
        return len(self.sequence)
    

    def __getitem__(self, index):
        """
        Allows indexing and slicing.
        Returns a new object of the same class if sliced, or a string character if indexed.
        """
        result = self.sequence[index]
        
        if isinstance(index, slice):
            return type(self)(result)
            
        return result
    

    def __str__(self) -> str:
        """
        Provides a formatted string representation of the sequence object.
        """
        return f"{self.__class__.__name__}: {self.sequence}"
    

class NucleicAcidSequence(BiologicalSequence):
    """
    Intermediate base class for nucleic acids (DNA and RNA).
    Provides common methods like complement, reverse, and reverse_complement.
    """

    @property
    @abstractmethod
    def complement_map(self) -> dict:
        """
        Abstract property. Forces child classes to provide their own
        dictionary for nucleotide complementation.
        """
        raise NotImplementedError("Subclasses must implement 'complement_map'")


    def complement(self):
        """
        Returns the complementary sequence as a new object of the exact same type.
        """
        comp_seq = "".join(self.complement_map.get(ch, ch) for ch in self.sequence)
        
        return type(self)(comp_seq)
    

    def reverse(self):
        """
        Returns the reversed sequence as a new object of the exact same type.
        """
        reversed_seq = self.sequence[::-1]
        return type(self)(reversed_seq)


    def reverse_complement(self):
        """
        Returns the reverse complement of the sequence.
        Demonstrates method chaining.
        """
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """
    Represents a DNA sequence.
    Inherits all sequence operations and complement logic from NucleicAcidSequence.
    """

    @property
    def alphabet(self) -> set:
        """Returns the valid DNA alphabet."""
        return {"A", "T", "G", "C"}


    @property
    def complement_map(self) -> dict:
        """Returns the complement mapping specific to DNA."""
        return {
            "A": "T", "T": "A", "G": "C", "C": "G",
            "a": "t", "t": "a", "g": "c", "c": "g"
        }


    def transcribe(self):
        """
        Transcribes the DNA sequence into an RNA sequence.
        Returns a new RNASequence object.
        """
        transcribed_seq = self.sequence.replace("T", "U").replace("t", "u")
        return RNASequence(transcribed_seq)
    

class RNASequence(NucleicAcidSequence):
    """
    Represents an RNA sequence.
    """

    @property
    def alphabet(self) -> set:
        """Returns the valid RNA alphabet."""
        return {"A", "U", "G", "C"}
    

    @property
    def complement_map(self) -> dict:
        """Returns the complement mapping specific to RNA."""
        return {
            "A": "U", "U": "A", "G": "C", "C": "G",
            "a": "u", "u": "a", "g": "c", "c": "g"
        }


class AminoAcidSequence(BiologicalSequence):
    """
    Represents an amino acid (protein) sequence.
    Inherits from BiologicalSequence, bypassing nucleic acid specific methods.
    """

    @property
    def alphabet(self) -> set:
        """Returns the valid standard amino acid alphabet."""
        return set("ACDEFGHIKLMNPQRSTVWY")
    

    def count_hydrophobic(self) -> int:
        """
        Calculates the number of hydrophobic amino acids in the sequence.
        """
        hydrophobic_aa = {"A", "V", "I", "L", "M", "F", "Y", "W"}
        return sum(1 for aa in self.sequence.upper() if aa in hydrophobic_aa)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[int, float, Tuple[Union[int, float], Union[int, float]]] = (0, 100),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """
    Filters a FASTQ file based on GC content, sequence length, and average quality.
    Uses Biopython's SeqIO and SeqRecord objects for parsing and writing.
    Saves the filtered records to a 'filtered' directory.
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    output_dir = Path.cwd() / "filtered"
    output_dir.mkdir(exist_ok=True)
    final_path = output_dir / output_fastq

    if final_path.exists():
        print(f"Warning: File {final_path} already exists. Write operation aborted.")
        return

    passed_records = []
    
    for record in SeqIO.parse(input_fastq, "fastq"):
        seq_len = len(record)
        if not (length_bounds[0] <= seq_len <= length_bounds[1]):
            continue  
        gc_content = gc_fraction(record.seq) * 100
        if not (gc_bounds[0] <= gc_content <= gc_bounds[1]):
            continue
        phred_scores = record.letter_annotations["phred_quality"]
        average_quality = sum(phred_scores) / seq_len if seq_len > 0 else 0
        if average_quality < quality_threshold:
            continue
        passed_records.append(record)

    SeqIO.write(passed_records, final_path, "fastq")


if __name__ == "__main__":
    print("Test OOP Biological Sequences")
    
    # DNA
    dna = DNASequence("ATGC")
    print(f"Original DNA: {dna}")  
    print(f"DNA length: {len(dna)}") 
    print(f"Is alphabet valid? {dna.check_alphabet()}")
    

    print(f"Complement: {dna.complement()}")
    print(f"Reverse complement: {dna.reverse_complement()}")
    

    rna = dna.transcribe()
    print(f"Transcribed to RNA: {rna}")
    print(f"RNA complement: {rna.complement()}")
    

    sliced_dna = dna[1:3]
    print(f"Sliced DNA [1:3]: {sliced_dna} (Type: {type(sliced_dna).__name__})")
    
    # Protein
    protein = AminoAcidSequence("MAVW")
    print(f"Protein: {protein}")
    print(f"Hydrophobic count: {protein.count_hydrophobic()}")