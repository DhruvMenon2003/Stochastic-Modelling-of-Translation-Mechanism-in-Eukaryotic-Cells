# Stochastic-Modelling-of-Translation-Mechanism-in-Eukaryotic-Cells

# Project Overview

This project aims to simulate the translation process in eukaryotic cells, focusing on how synonymous codon usage influences protein folding. By modeling the stochastic nature of translation, we seek to predict protein folding outcomes based solely on spliced mRNA sequences.

# Background

The genetic code's degeneracy allows multiple codons to encode the same amino acid. However, synonymous codons are not utilized uniformly; this phenomenon is known as codon usage bias. Rare codons can slow translation, affecting co-translational protein folding and potentially leading to misfolded proteins.[1]

# Notably, proteins with higher intrinsic disorder tendencies are more susceptible to the effects of codon usage bias. The IUPred2A tool can predict such intrinsically disordered regions from amino acid sequences[2]
# Objectives

Simulate Translation: Model the translation process, incorporating codon usage probabilities to generate mRNA sequences.
Analyze Folding Outcomes: Assess the impact of rare codon frequency on protein folding, categorizing outcomes as properly folded, misfolded, or no protein produced.
Optimize Parameters: Determine optimal sequence lengths and rare codon thresholds to predict folding outcomes accurately.
# Methodology

# Codon Probability Matrix: Utilize empirical codon probabilities to reflect natural codon usage. These probabilities are normalized to sum to one, ensuring they represent valid probabilities.[3]

# Sequence Generation: Generate random mRNA sequences based on the codon probability matrix. This stochastic approach mirrors the natural variability in mRNA sequences.

# Outcome Analysis: Evaluate each sequence for:

Presence of start (AUG) and stop codons (UAA, UAG, UGA).
Frequency of rare codons (defined as those with probabilities below a specified threshold).
Classify outcomes:
Properly Folded Protein: Sequence contains start and stop codons with rare codon frequency below the threshold.
Misfolded Protein: Sequence contains start and stop codons with rare codon frequency above the threshold.
No Protein: Sequence lacks either start or stop codons.
# Parameter Optimization: Conduct simulations across various sequence lengths and rare codon thresholds. Generate a 2D contour plot to visualize the relationship between these parameters and folding outcomes, aiding in identifying optimal conditions.

# Considerations

Simplifications: The model does not account for factors like temperature, pH, or chaperone proteins, which also influence protein folding.
Future Enhancements: Incorporating additional biological factors and refining the codon probability matrix with more comprehensive data could improve model accuracy.
# References
[1]Liu, Y. A code within the genetic code: codon usage regulates co-translational protein folding. Cell Commun Signal 18, 145 (2020). https://doi.org/10.1186/s12964-020-00642-6
#
[2]https://iupred2a.elte.hu/
#
[3]Yoon BJ. Hidden Markov Models and their Applications in Biological Sequence Analysis. Curr Genomics. 2009 Sep;10(6):402-15. doi: 10.2174/138920209789177575. PMID: 20190955; PMCID: PMC2766791.
