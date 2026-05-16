---
layout: single
title: "Methods"
permalink: /methods/
---

Thermostable and non-thermostable laccases were collected from UniProt, 
NCBI and BRENDA using the filters these databases provide. 
The protein sequences were downloaded in FASTA format using a bash 
script.

To reduce sequence redundancy within the groups, the sequences were 
clustered using CD-HIT into smaller representative groups.

Multiple sequence alignments of thermostable and non-thermostable 
laccases were then constructed using MAFFT. To improve their 
quality and downstream usability, the alignments were trimmed with trimAl.

The solvent accessibility of each residue in the B. subtilis CotA structure (PDB: 1GSK) was calculated using FreeSASA. Residues were classified as buried or solvent-accessible based on a 20% relative solvent accessibility threshold and exported from ChimeraX.

> ##### ***What is solvent accessibility?***
> Solvent accessibility measures how exposed a given amino acid is to the surrounding water molecules. It is expressed as a percentage of the maximum possible exposure for that amino acid — a value close to 0% means the residue is deeply buried inside the protein, while a value close to 100% means it is fully exposed on the surface. This distinction is important because buried and surface residues play very different roles in protein structure and stability.

Position Specific Scoring Matrices (PSSMs) were generated from the trimmed alignments using goalign. Using CotA as a reference, the PSSM rows corresponding to buried and solvent-accessible residues were extracted separately for both datasets using a custom Python script.

Three complementary analyses were then applied to compare the two groups. Amino acid group enrichment was assessed using Chi-squared tests, physicochemical properties were compared using Student's t-test and Mann-Whitney U tests, and the Kullback-Leibler divergence was used to identify the positions that differ most between thermostable and non-thermostable laccases.

The most discriminating positions were mapped onto the 1GSK structure using ChimeraX and their structural context was examined to interpret their potential role in thermostability.
