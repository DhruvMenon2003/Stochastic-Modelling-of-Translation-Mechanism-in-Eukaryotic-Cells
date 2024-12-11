from Bio import *  # Import all BioPython functionalities
from Bio.Seq import Seq
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio import SeqIO, Entrez
import numpy as np
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import os

class CodonUsageIndex:
    def __init__(self):
        self.codon_counts = {}
        self.rscu_values = {}

    def generate_index(self, sequence):
        total_codons = 0
        codon_table = unambiguous_dna_by_id[1].forward_table

        # Initialize codon counts
        for codon in codon_table.keys():
            self.codon_counts[codon] = 0

        # Count codons
        codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
        for codon in codons:
            if codon in self.codon_counts:
                self.codon_counts[codon] += 1
            total_codons += len(codons)

        # Calculate RSCU values
        for aa in set(codon_table.values()):
            synonymous_codons = [codon for codon, aa_mapped in codon_table.items() if aa_mapped == aa]
            aa_total = sum(self.codon_counts[codon] for codon in synonymous_codons)
            n_codons = len(synonymous_codons)
            for codon in synonymous_codons:
                expected = aa_total / n_codons if aa_total > 0 else 0
                self.rscu_values[codon] = self.codon_counts[codon] / expected if expected > 0 else 0

        print("Codon Usage Index generated successfully.")

    def display_rscu_values(self):
        for codon, value in self.rscu_values.items():
            print(f"{codon}: {value:.2f}")

    def classify_protein(self, disorder_scores, rare_threshold=0.3):
        disorder_median = np.median(disorder_scores)
        rare_codon_proportion = sum(1 for value in self.rscu_values.values() if value < rare_threshold) / len(self.rscu_values)

        if disorder_median > 0.5 and rare_codon_proportion > 0.2:
            return "Class 1"  # Disordered and co-translationally folded
        elif disorder_median <= 0.5 and rare_codon_proportion > 0.2:
            return "Class 2"  # Ordered and co-translationally folded
        elif disorder_median > 0.5 and rare_codon_proportion <= 0.2:
            return "Class 3"  # Disordered and not co-translationally folded
        else:
            return "Class 4"  # Ordered and not co-translationally folded

# Fetch protein sequences from UniProt using Entrez
Entrez.email = "dhruv.menon@plaksha.edu.in"  # Replace with your email

def fetch_sequences_from_uniprot(term, max_results=5):
    search_handle = Entrez.esearch(db="protein", term=term, retmax=max_results)
    search_results = Entrez.read(search_handle)
    search_handle.close()
    protein_ids = search_results.get("IdList", [])
    fetch_handle = Entrez.efetch(db="protein", id=protein_ids, rettype="fasta", retmode="text")
    records = list(SeqIO.parse(fetch_handle, "fasta"))
    fetch_handle.close()
    return records

# Example usage
if __name__ == "__main__":
    # Fetch sequences directly from UniProt
    proteins = fetch_sequences_from_uniprot("co-translational folding", max_results=100)
    codon_usage_index = CodonUsageIndex()
    for protein in proteins:
        sequence = str(protein.seq)
        codon_usage_index.generate_index(sequence)
        codon_usage_index.display_rscu_values()
        disorder_scores = np.random.uniform(0.2, 0.7, len(protein.seq))  # Simulated disorder scores
        classification = codon_usage_index.classify_protein(disorder_scores)
        print(f"Protein {protein.id} classified as: {classification}")
