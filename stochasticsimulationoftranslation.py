# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 11:49:41 2024

@author: lenovo
"""
import numpy as np
import random

# Normalize codon probabilities
def normalize_codon_probs(codons):
    total = sum(codons.values())
    return {k: v / total for k, v in codons.items()}

# Draw a random codon sequence
def generate_sequence(codons, length):
    codon_keys = list(codons.keys())
    codon_probs = list(codons.values())
    return random.choices(codon_keys, weights=codon_probs, k=length)

# Check outcomes
def analyze_sequence(sequence, start_codon, stop_codons, rare_codons, rare_threshold):
    # Check for start and stop codons
    has_start = start_codon in sequence
    has_stop = any(stop in sequence for stop in stop_codons)
    
    if not has_start or not has_stop:
        return "No Protein"
    
    # Count rare codons
    rare_count = sum(1 for codon in sequence if codon in rare_codons)
    if rare_count >= rare_threshold:
        return "Misfolded Protein"
    
    return "Properly Folded Protein"

# Simulation
def simulate(codons, num_sequences, seq_length, rare_threshold):
    codons = normalize_codon_probs(codons)
    
    # Define key codons
    start_codon = "AUG"
    stop_codons = {"UAA", "UAG", "UGA"}
    rare_codons = {k for k, v in codons.items() if v < 0.01}  # Rare codons
    
    outcomes = {"Properly Folded Protein": 0, "No Protein": 0, "Misfolded Protein": 0}
    
    # Run simulations
    for _ in range(num_sequences):
        sequence = generate_sequence(codons, seq_length)
        outcome = analyze_sequence(sequence, start_codon, stop_codons, rare_codons, rare_threshold)
        outcomes[outcome] += 1
    
    # Normalize outcomes to probabilities
    total = sum(outcomes.values())
    outcomes = {k: v / total for k, v in outcomes.items()}
    return outcomes

# Example Input
hmmcodons = {
    'AAA': 0.008655332302936622, 'AAU': 0.008655332302936622, 'AAG': 0.008655332302936622, 'AAC': 0.023080886141164325,
    'AUA': 0.008655332302936622, 'AUU': 0.008655332302936622, 'AUG': 0.010819165378670777, 'AUC': 0.023080886141164325,
    'AGA': 0.008655332302936622, 'AGU': 0.008655332302936622, 'AGG': 0.008655332302936622, 'AGC': 0.008655332302936622,
    'ACA': 0.008655332302936622, 'ACU': 0.008655332302936622, 'ACG': 0.008655332302936622, 'ACC': 0.023080886141164325,
    'UAA': 0.09891808346213285, 'UAU': 0.008655332302936622, 'UAG': 0.09891808346213285, 'UAC': 0.008655332302936622,
    'UUA': 0.008655332302936622, 'UUU': 0.023080886141164325, 'UUG': 0.008655332302936622, 'UUC': 0.008655332302936622,
    'UGA': 0.09891808346213285, 'UGU': 0.008655332302936622, 'UGG': 0.023080886141164325, 'UGC': 0.008655332302936622,
    'UCA': 0.008655332302936622, 'UCU': 0.008655332302936622, 'UCG': 0.008655332302936622, 'UCC': 0.008655332302936622,
    'GAA': 0.008655332302936622, 'GAU': 0.023080886141164325, 'GAG': 0.008655332302936622, 'GAC': 0.008655332302936622,
    'GUA': 0.008655332302936622, 'GUU': 0.008655332302936622, 'GUG': 0.008655332302936622, 'GUC': 0.008655332302936622,
    'GGA': 0.008655332302936622, 'GGU': 0.023080886141164325, 'GGG': 0.008655332302936622, 'GGC': 0.023080886141164325,
    'GCA': 0.008655332302936622, 'GCU': 0.008655332302936622, 'GCG': 0.008655332302936622, 'GCC': 0.008655332302936622,
    'CAA': 0.023080886141164325, 'CAU': 0.008655332302936622, 'CAG': 0.023080886141164325, 'CAC': 0.008655332302936622,
    'CUA': 0.008655332302936622, 'CUU': 0.008655332302936622, 'CUG': 0.023080886141164325, 'CUC': 0.008655332302936622,
    'CGA': 0.008655332302936622, 'CGU': 0.008655332302936622, 'CGG': 0.008655332302936622, 'CGC': 0.008655332302936622,
    'CCA': 0.008655332302936622, 'CCU': 0.023080886141164325, 'CCG': 0.008655332302936622, 'CCC': 0.008655332302936622,
}

# Run Simulation
num_sequences = 10000
sequence_length = 550  # Average sequence length
rare_threshold = 0.45*sequence_length  

outcomes = simulate(hmmcodons, num_sequences, sequence_length, rare_threshold)
print(outcomes)

