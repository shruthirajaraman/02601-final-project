# Quantifying the Phylogeny of the Presence of Antibiotic Resistance in Bacterial Antibiotic Resistance Genes (ARGs)
Isabelle, Shruthi, and Charlotte's 02613 Programming for Scientists Final Project.

Our project aims to identify ARG regions in bacterial genomes and create a phylogenetic tree from the identified ARG regions. We determine the preservation of the ARG-coded proteins throughout the tree and estimate the likelihood that each ARG region displays antibiotic resistance.

This repository contains source code and example data to run the RShiny UI for our project.

PREPARING YOUR MACHINE
----
1. Download [RStudio](https://posit.co/download/rstudio-desktop/) on your machine. Utilize version 2024.04.2 to ensure compatability with key R packages.
2. If desired, download [Visual Studio Code](https://code.visualstudio.com/download) (or your preferred code editor) to view/edit source code.
3. Download [Go language](https://go.dev/dl/) on your machine to ensure the source code can be run successfully.
4. Download source code (Group Project Rshiny App.R, main.go) and example data (genome_sequences_small_dataset.fasta, ciproflocaxin_gene_sequence.fasta) from this Github repository. Make sure that the source code is contained in a folder titled P4S-term-project. 

RUNNING THE PROGRAM
----
1. Open Group Project RShiny App.R (contained in the P4S-term-project folder you created above) in RStudio. 
2. Set working directory to source file location (Session > Set Working Directory > To Source File Location).
3. Click Run App.
4. Upload a file of genome sequences at the "Upload your file of genomes" prompt. This file must be in the FASTA format. It can either be a .txt or .fasta file. If utilizing the example data, upload genome_sequences_small_dataset.fasta. 
Note: file limit set to 2 GB. The program will run slowly (i.e. take hours) if your file is towards the upper limit.
5. Upload a file of an ARG sequence at the "Upload your antibiotic resistance gene sequence file" prompt. This file must be in the FASTA format. It can either be a .txt or .fasta file. If utilizing the example data, upload ciprofloxacin_gene_sequence.fasta. 
6. Set your desired percent identity threshold. Percent identity is the nucleotide similarity between two regions. Only genomes with a gene region with a nucleotide similarity at or above this threshold, compared to the input gene sequence, will be included in the calculations.
7. Set your desired mutation percentage threshold. Mutation percentage is the number of nucleotides in an identified gene region that do not match the input gene sequence. Genomes with a mutation percentage at or above this threshold will be considered more likely to have antibiotic resistance, and vice versa.
8. Click Run Genome Processing. Your trees and data tables will appear on their respective tabs when the program is complete.

OUTPUTS
----
- Tree and table visualizations in RShiny.
- Newick (NHX) representations of nucleotide and protein sequence trees created from the input genomes. Only genomes that have a percent identity at or above the specified threshold, compared to the input gene sequence, will be included in the trees.
- FASTA files of the nucleotide and protein sequences of the identified ARG regions in the genomes. Only genomes that have a percent identity at or above the specified threshold, compared to the input gene sequence, will have their sequences included in these files.
- CSV files with each genome's protein preservation data and predicted likelihood of antibiotic resistance. Protein preservation is calculated by the number of non-mutated nucleotides divided by the length of the protein, comparing to the input gene sequence. For antibiotic resistance likelihood prediction, the protein preservation value is compared to the input preservation threshold; if the protein preservation value is less than this value, we consider this genome to be more likely to code for antibiotic resistance (represented by a 1). Only genomes that have a percent identity at or above the specified threshold, compared to the input gene sequence, will have their sequences included in these files.

DEMO
----
You can find a recorded demonstration of this program ![here](GroupProjectDemo.mp4). You'll need to click the "View Raw" hyperlink on that page to download the demo video and watch it.

COMMON GOTCHAS
----
- If re-running the application, make sure to remove any previous output files created by the program in your working directory folder. This will ensure that your new run will work off the new output files.
- If no trees are created and you get a "panic: runtime error: makeslice: len out of range" error in the RStudio terminal, this is most likely because there aren't any genomes in your file with a percent identity at or above your set threshold. Please decrease the percent identity threshold until a tree is successfully created.
- If you have a very large genome dataset (thousands of genomes), the program will likely take a long time to run. Patience is key.

CONTRIBUTORS
----
- Isabelle D'Amico
- Lalithashruthi Rajaraman
- Charlotte Barretto

TESTS
----
Test files are included in the tests folder on this Github repository. To run them, please open main.go and functions_test.go in Visual Studio. First compile the programs with "go build" in the terminal. Then, run "go test" in the terminal.