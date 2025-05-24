package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"
)

type DistanceMatrix [][]float64

// ARG is a datatype that holds the information for an antibiotic resistance gene region identified in a genome.
type ARG struct {
	Nucleotides   string
	Protein       string
	GenomeStrain  string
	GenomeName    string
	IndexInGenome [2]int
	Mutations     int
	NumMutations  int
	ArgLabel      string
}

// ARG is a datatype that holds the information for the protein of a antibiotic resistance gene region identified in a genome.
type Protein struct {
	AminoAcids          string
	ARG                 *ARG
	Polarity            string
	Mutations           int
	NumNotableMutations int
}

// ARGNode is a datatype that holds the information for a ARG's respective node in the nucleotide ARG tree.
type ARGNode struct {
	ARGpointer   *ARG
	Protein      *Protein
	Preservation float64
	Likelihood   float64
	GenomeStrain string
	Descendents1 *ARGNode
	Descendents2 *ARGNode
}

// ProteinNode is a datatype that holds the information for a Protein's respective node in the protein ARG tree.
type ProteinNode struct {
	Protein      *Protein
	Label        string
	Preservation float64
	Descendants1 *ProteinNode
	Descendants2 *ProteinNode
}

// ARGTree is a datatype that holds the information for the nucleotide ARG tree, made up of ARGNodes.
type ARGTree []*ARGNode

// ProteinTree is a datatype that holds the information for the protein ARG tree, made up of ProteinNodes.
type ProteinTree []*ProteinNode

func main() {
	//Check that we have the right amount of arguments
	if len(os.Args) < 5 {
		panic("Not enough arguments provided.")
	}

	//Save off the argument values
	genomeFilePath := os.Args[1]
	geneFilePath := os.Args[2]
	percentIdentity, err := strconv.ParseFloat(os.Args[3], 64)
	if err != nil {
		panic(err)
	}
	mutationThreshold, err := strconv.ParseFloat(os.Args[4], 64)
	if err != nil {
		panic(err)
	}

	//create our slices of genome sequences and names
	genomes, genomeNames := ReadGenomesFromFile(genomeFilePath)
	fmt.Println("Read in genomes.")

	//save our searchGene
	searchGene, searchGeneName := ReadGeneFromFile(geneFilePath)
	fmt.Println("Read in gene.")

	//run IdentifyARGsForGenomes to pull the best ARG sequence in each genome and save the extended ARG regions to a slice of strings, and the indices to remove from our genome slices
	extendedARGRegions, _, removeIndices := IdentifyARGsForGenomes(genomes, searchGene, percentIdentity)
	fmt.Println("Identified ARG regions in genomes.")

	//run RemoveGenomes to remove any genomes that don't have the ARG region from our slice
	genomes, genomeNames = RemoveGenomes(genomes, genomeNames, removeIndices)

	//run TranslateProteinSequences on slice of extended ARG regions to create a slice of proteins
	proteins := TranslateProteinSequences(extendedARGRegions)
	fmt.Println("Translated ARG regions to proteins.")

	//run NucleotideDistanceMatrix and ProteinDistanceMatrix on the slice of extended ARG sequences and slice of proteins
	nucleotideDistMtx := NucleotideDistanceMatrix(extendedARGRegions)
	proteinDistMtx := ProteinDistanceMatrix(proteins)
	fmt.Println("Created distance matrices.")

	//run CreateNucleotideARGTree on nucleotide distance matrix
	argTree := CreateNucleotideARGTree(nucleotideDistMtx, genomes, genomeNames, searchGene, searchGeneName, percentIdentity)
	fmt.Println("Created nucleotide tree.")

	//run CreateProteinTree on protein distance matrix
	proteinTree := CreateProteinTree(proteinDistMtx, proteins, genomes, genomeNames, argTree, percentIdentity)
	fmt.Println("Created protein tree.")

	//run CalculatePreservationOfNucleotides on nucleotide tree
	argTree.CalculatePreservationOfNucleotides(searchGene)

	//run CalculatePreservationOfProteins on protein tree
	proteinTree.CalculatePreservationOfProteins(searchGene)
	fmt.Println("Calculated preservation.")

	//run CalculateLikelihood on preservation calculations
	argTree.CalculateLikelihood(mutationThreshold)
	fmt.Println("Calculated likelihood.")

	//summarize data with Newick() on both Protein and nucleotide tree
	WriteNewickToFile(("(" + argTree[len(argTree)-1].NewickNucleotide() + ")" + ";"), "Nucleotide")
	WriteNewickToFile(("(" + proteinTree[len(proteinTree)-1].NewickProtein() + ")" + ";"), "Protein")
	fmt.Println("Wrote trees to file.")

	//Create our FASTA output file for nucleotides
	nucleotideFile, err := os.Create("NucleotideFastaOutput.txt")
	if err != nil {
		panic(err)
	}

	// Ensure the file is closed when we're done
	defer nucleotideFile.Close()
	argTree[len(argTree)-1].WritetoNucleotideFasta(nucleotideFile)
	fmt.Println("Wrote nucleotides to file.")

	//Create our FASTA output file for proteins
	proteinFile, err2 := os.Create("ProteinFastaOutput.txt")
	if err2 != nil {
		panic(err2)
	}

	// Ensure the file is closed when we're done
	defer proteinFile.Close()
	proteinTree[len(proteinTree)-1].WritetoProteinFasta(proteinFile)
	fmt.Println("Wrote proteins to file.")

	//Create our likelihood and preservation csv files.
	argTree[len(argTree)-1].WriteTreeToCSV("NucleotideLikelihoodTable.csv")
	proteinTree[len(proteinTree)-1].WriteTreeToCSV("ProteinPreservationTable.csv")
	fmt.Println("Wrote preservation and likelihood to file.")
}

// Source for Codon Table:
// https://www.khanacademy.org/science/ap-biology/gene-expression-and-regulation/translation/a/the-genetic-code-discovery-and-properties

// TranslateProteinSequences (Shruthi) (edited by Isabelle)
// INPUT: a slice of strings dnaFragments consisting of DNA sequences corresponding to identified ARGs in the original genome
// OUTPUT: a slice of strings consisting of the amino acid sequences that the ARG sequences code for
func TranslateProteinSequences(dnaFragments []string) []string {
	proteinSeq := make([]string, 0)

	for _, fragment := range dnaFragments {
		var protein strings.Builder
		// Transcribe the DNA into mRNA
		rnaFragment := strings.Replace(fragment, "T", "U", -1)
		for i := 0; i < len(rnaFragment); i += 3 {
			// Translate each codon of length 3 into its corresponding amino acid to create the protien sequence
			if i+3 <= len(rnaFragment) {
				aminoAcid := RNAtoAminoAcid(rnaFragment[i : i+3])
				protein.WriteString(aminoAcid)
			}
		}
		//append each protein sequence to the slice
		proteinSeq = append(proteinSeq, protein.String())
	}

	return proteinSeq
}

// RNAtoAminoAcid (Shruthi) (edited by Isabelle)
// INPUT: a string of length three representing a codon in a mRNA sequence
// OUTPUT: a byte representing the one-letter code for the amino acid that a codon translates to
func RNAtoAminoAcid(codon string) string {
	aminoAcid := ""

	switch codon {
	case "UUU", "UUC":
		aminoAcid = "F" // Phenylalanine
		// polarity = "nonpolar"
	case "UUA", "UUG", "CUU", "CUC", "CUA", "CUG":
		aminoAcid = "L" // Leucine
		// polarity = "nonpolar"
	case "AUU", "AUC", "AUA":
		aminoAcid = "I" // Isoleucine
		// polarity = "nonpolar"
	case "AUG":
		aminoAcid = "M" // Methionine
		// polarity = "nonpolar"
	case "GUU", "GUC", "GUA", "GUG":
		aminoAcid = "V" // Valine
		// polarity = "nonpolar"
	case "UCU", "UCC", "UCA", "UCG", "AGU", "AGC":
		aminoAcid = "S" // Serine
		// polarity = "polar"
	case "CCU", "CCC", "CCA", "CCG":
		aminoAcid = "P" // Proline
		// polarity = "nonpolar"
	case "ACU", "ACC", "ACA", "ACG":
		aminoAcid = "T" // Threonine
		// polarity = "polar"
	case "GCU", "GCC", "GCA", "GCG":
		aminoAcid = "A" // Alanine
		// polarity = "nonpolar"
	case "UAU", "UAC":
		aminoAcid = "Y" // Tyrosine
		// polarity = "polar"
	case "UAA", "UAG", "UGA":
		aminoAcid = "*" // Stop codon
		// polarity = "nonpolar"
	case "CAU", "CAC":
		aminoAcid = "H" // Histidine
		// polarity = "basic"
	case "CAA", "CAG":
		aminoAcid = "Q" // Glutamine
		// polarity = "polar"
	case "AAU", "AAC":
		aminoAcid = "N" // Asparagine
		// polarity = "polar"
	case "AAA", "AAG":
		aminoAcid = "K" // Lysine
		// polarity = "basic"
	case "GAU", "GAC":
		aminoAcid = "D" // Aspartic Acid
		// polarity = "acidic"
	case "GAA", "GAG":
		aminoAcid = "E" // Glutamic Acid
		// polarity = "acidic"
	case "UGU", "UGC":
		aminoAcid = "C" // Cysteine
		// polarity = "nonpolar"
	case "UGG":
		aminoAcid = "W" // Tryptophan
		// polarity = "nonpolar"
	case "CGU", "CGC", "CGA", "CGG", "AGA", "AGG":
		aminoAcid = "R" // Arginine
		// polarity = "basic"
	case "GGU", "GGC", "GGA", "GGG":
		aminoAcid = "G" // Glycine
		// polarity = "nonpolar"
	}

	return aminoAcid
}

// ProteinDistanceMatrix (Shruthi), edited - Charlotte 12/3/24
// INPUT: a slice of AminoAcid objects representing the protein sequences of ARGs
// OUTPUT: a 2D array matrix consisting of float64 values representing the dissimilarity between protein sequences
func ProteinDistanceMatrix(proteinSeq []string) DistanceMatrix {
	// Initialize the distance matrix
	numSamples := len(proteinSeq)
	mtx := InitializeSquareMatrix(numSamples)

	// Populate the distance matrix with the percent protein similarity
	for i := range proteinSeq {
		for j := i + 1; j < len(proteinSeq); j++ {
			// When comparing a sequence to itself, the distance should be 0
			if i == j-1 {
				mtx[i][i] = 0
			}
			// Calculate the percent similarity
			percentSimilarity := PercentProteinIdentity(proteinSeq[i], proteinSeq[j])
			// Since the matrix is symmetrical across the diagonal, mtx[i][j] = mtx[j][i]
			mtx[i][j] = percentSimilarity
			mtx[j][i] = percentSimilarity
		}
	}

	return mtx
}

// created - idamico 11/11/24 (edited Charlotte)

// IdentifyARGsForGenomes takes as input a slice of genomes, a searchGene string, and a float percentIdentity.
// It returns a slice of extended ARG regions, a slice of the gene indices corresponding to the original ARG region in each genome, and a slice of indices for genomes to remove (they didn't have any ARGs)
func IdentifyARGsForGenomes(genomes []string, searchGene string, percentIdentity float64) ([]string, [][2]int, []int) {
	//create our slices to save the extended ARG regions and ARG indices
	ARGregions := make([]string, 0, len(genomes))
	geneIndicesSlice := make([][2]int, 0, len(genomes))
	genomesToRemove := make([]int, 0)

	//range over all the genomes and identify the best ARG region
	for i, val := range genomes {
		geneIndices := IdentifyARGs(val, searchGene, percentIdentity)

		//if we don't find any ARGs, let's make a note of the index in the genome slice so we can remove it later
		if geneIndices[0][0] == -1 {
			genomesToRemove = append(genomesToRemove, i)
			continue
		}

		//if we do find ARGs, we'll save their indices and extended regions
		geneIndicesSlice = append(geneIndicesSlice, geneIndices...)
		ARGregions = append(ARGregions, SaveArgRegions(val, geneIndices)...)

	}

	//return all our slices
	return ARGregions, geneIndicesSlice, genomesToRemove
}

// created - idamico 11/11/24
// RemoveGenomes takes as input a slice of genome strings, a slice of genome names, and a []int slice of indices.
// It returns an updated slice of genome strings amd updated slice of genome name strings, minus the indices specified.
func RemoveGenomes(genomes, genomenames []string, indices []int) ([]string, []string) {
	//reverse sort our indices to remove so we start with the largest first
	//this is to prevent any index out of bound errors
	sort.Sort(sort.Reverse(sort.IntSlice(indices)))

	//range over the indices and remove them from our genome slices
	for _, val := range indices {
		//if the index is at the end, need special handling
		if val == len(genomes)-1 {
			genomes = genomes[:val]
			genomenames = genomenames[:val]
		} else { //otherwise, remove as normal
			genomes = append(genomes[:val], genomes[val+1:]...)
			genomenames = append(genomenames[:val], genomenames[val+1:]...)
		}
	}

	//return our slices
	return genomes, genomenames
}

// created - idamico 10/25/24, edited - Charlotte 11/1/24
// IdentifyARGs takes as input a genome string, a searchGene string, and a float percentIdentity.
// It returns a [][2] matrix corresponding to the start/end indices of the region in the genome that has the highest nucleotide percentIdentity with the searchGene.
// If there are no regions in the genome that match the searchGene, we will return [[-1,-1]]
func IdentifyARGs(genome, searchGene string, percentIdentity float64) [][2]int {
	//initialize our variables and slices
	genomeLength := len(genome)
	geneLength := len(searchGene)
	var geneIndices [][2]int
	var tempArray [2]int

	//set our max to 0 for noww
	max := 0.0

	//search over the length of the genome
	for i := 10; i <= genomeLength-geneLength-10; i++ {
		//in windows of the length of searchGene, evaluate the nucleotide percent identity against searchGene
		percentNdId := PercentNucleotideIdentity(searchGene, genome[i:i+geneLength])

		//if the nucleotide percent identity is at or above the given value, AND it's the maximum percent identity so far, we save the start/end indices to the result array
		if (percentNdId >= percentIdentity && percentNdId > max) || (percentNdId >= percentIdentity && max == 0.0) {
			max = percentNdId
			tempArray = [2]int{i, i + geneLength - 1}
		}

	}

	//if we did not find a sequence that matched, we update the indices to -1,-1
	if max == 0.0 {
		tempArray = [2]int{-1, -1}
	}

	//append the highest percent identity ARG region to the return array and return
	geneIndices = append(geneIndices, tempArray)
	return geneIndices
}

// created - idamico 10/25/24; edited - Shruthi 12/1/24
// PercentNucleotideIdentity takes as input two gene strings.
// It returns the percentage (float) of nucleotides that are the same in both gene sequences.
func PercentNucleotideIdentity(gene1, gene2 string) float64 {
	//initialize our variables
	countSame := 0
	geneLength := 0

	//if gene 2 is longer, we'll only check the length of gene 1
	if len(gene1) <= len(gene2) {
		geneLength = len(gene1)

		//range over the genes and count the number of the same nucleotides
		for i := range gene1 {
			if gene1[i] == gene2[i] {
				countSame++
			}
		}
	} else { //if gene 1 is longer, we'll only check the length of gene 2
		geneLength = len(gene2)

		//range over the genes and count the number of the same nucleotides
		for i := range gene2 {
			if gene1[i] == gene2[i] {
				countSame++
			}
		}
	}

	//return the percentage of nucleotides that are the same in both genes
	return float64(countSame) / float64(geneLength)
}

// created - idamico 10/29/24; edited - Shruthi 12/1/24
// PercentProteinIdentity takes as input two protein strings
// It returns the percentage (float) of amino acids that are the same in both protein sequences
func PercentProteinIdentity(protein1, protein2 string) float64 {
	//initialize our variables
	countSame := 0
	proteinLength := 0

	// if protein 2 is longer, we'll only check the length of protein 1
	if len(protein1) <= len(protein2) {
		proteinLength = len(protein1)

		// Range over the proteins and count the number of same amino acids
		for i := range protein1 {
			if protein1[i] == protein2[i] {
				countSame++
			}
		}
	} else { //if protein 1 is longer, we'll only check the length of protein 2
		proteinLength = len(protein2)

		// Range over the proteins and count the number of same amino acids
		for i := range protein2 {
			if protein2[i] == protein1[i] {
				countSame++
			}
		}
	}

	// Return the percentage of amino acids that are the same in both proteins
	return float64(countSame) / float64(proteinLength)
}

// created - idamico 10/29/24, edited -Charlotte 12/3/24
// NucleotideDistanceMatrix takes as input a slice of nucleotide strings
// It returns a distance matrix for the nucleotide strings, taking the extended region into account
func NucleotideDistanceMatrix(nucleotides []string) DistanceMatrix {
	//initialize distance matrix
	numSamples := len(nucleotides)
	mtx := InitializeSquareMatrix(numSamples)

	//range through distance matrix and set values based on percent nucleotide identity
	for i := range nucleotides {
		for j := i + 1; j < len(mtx[i]); j++ {
			mtx[i][j] = PercentNucleotideIdentity(nucleotides[i], nucleotides[j])

			//set the value that is symmetric across the main diagonal
			mtx[j][i] = mtx[i][j]
		}
	}

	return mtx
}

// created - idamico 10/29/24, edited - Charlotte 12/3/24
// InitializeSquareMatrix takes as input an integer corresponding to the number of rows and columns.
// It returns a matrix with the specified number of row and columns.
func InitializeSquareMatrix(i int) DistanceMatrix {
	//create matrix rows
	mtx := make([][]float64, i)

	//create matrix columns
	for row := range mtx {
		mtx[row] = make([]float64, i)
	}

	return mtx
}

// Charlotte Barretto 11/1/24 (edited by Isabelle)
// extend the region of the arg as a string
// Input: Genome as a string, and the indicies of the found ARG region
func SaveArgRegions(genome string, args [][2]int) []string {
	extended := make([]string, 0)
	//we want to extend the region to account for any frame shifts that may have occured between genome strains
	//the percent identity will never be 100% because of this, except for those strains that are highly conserved in the ARG region
	for indicePair := range args {
		if args[indicePair][0]-10 <= 0 {
			args[indicePair][0] = 0
		} else {
			args[indicePair][0] -= 10
		}

		if args[indicePair][1]+10 >= len(genome)-1 {
			args[indicePair][1] = len(genome) - 1
		} else {
			args[indicePair][1] += 10
		}

		extended = append(extended, genome[args[indicePair][0]:args[indicePair][1]])
	}

	return extended
}

// Charlotte Barretto 11/10/24 (edited by Isabelle)
// Input: distance matrix, ammino acid sequences, genomeStrains, genomeNames, ARGTree, percentIdentity
// Output: ProteinTree updated using the clusters
func CreateProteinTree(mtx DistanceMatrix, AASequences []string, GenomeStrains, GenomeNames []string, argTree ARGTree, percentIdentity float64) ProteinTree {
	CheckSquareMatrixAndSameNumSpecies(mtx, AASequences)
	//the last node will have the greatest percent identity to the search string
	proteinTree := InitializeProteinTree(AASequences, GenomeStrains, GenomeNames, argTree)

	numLeaves := len(mtx)
	clusters := proteinTree.InitializeClustersForProtein()
	//traverse internal nodes of ARGTree
	for i := numLeaves; i < 2*numLeaves-1; i++ {
		//find the minimum value in the identity matrix and return the associated row and col
		row, col, _ := MinIdentityDifference(mtx)

		// set the decendents of each node as the row and col of the distance mtx
		proteinTree[i].Descendants1 = clusters[row]
		proteinTree[i].Descendants2 = clusters[col]

		//size of clusters
		clusters1 := CountLeavesProtein(proteinTree[i].Descendants1)
		clusters2 := CountLeavesProtein(proteinTree[i].Descendants2)

		mtx = AddRowCol(row, col, clusters1, clusters2, mtx)
		mtx = DeleteRowCol(mtx, row, col)

		clusters = append(clusters, proteinTree[i])
		clusters = DeleteClustersProtein(clusters, row, col)

	}
	return proteinTree
}

// Charlotte Barretto 11/10/24 (edited by Isabelle)
// Input: ammino acids as array of strings AASequences []string, the array of GenomeStrains, GenomeNames []string,
//
//	and the args ARGTree to use to find the matching GenomeNames between the ARGNodes and ProteinNodes
//
// Output: ProteinTree skeleton where only the leaves are set
func InitializeProteinTree(AASequences []string, GenomeStrains, GenomeNames []string, args ARGTree) ProteinTree {
	//check that there are the same number of AAsequences as GenomeStrains
	if len(AASequences) != len(GenomeStrains) {
		panic("The number of amino acid sequences are not the same as the number of genome strains.")
	}
	numLeaves := len(AASequences)
	//double check that there are the same number of AAsequences as GenomeNames
	if numLeaves != len(GenomeNames) {
		panic("Number of genomes does not match number of names. Incorrect input format")
	}
	tree := make([]*ProteinNode, 2*numLeaves-1)

	for i := range tree {

		if i < numLeaves {
			//get the corresponding ARG strain to linke the nodes of the two trees
			correspondingARG := FindARGInTree(GenomeStrains[i], args)
			tree[i] = &ProteinNode{}
			protein := &Protein{AminoAcids: AASequences[i], ARG: correspondingARG}
			correspondingARG.GenomeName = GenomeNames[i]
			tree[i].Protein = protein
		} else {
			//there will be an ancestor whenever i > than the num leaves
			tree[i] = &ProteinNode{}
			tree[i].Label = "Ancestor Protein: " + strconv.Itoa(i)

		}

	}
	return tree
}

// Charlotte Barretto 11/10/24
// Used during the inital creation of the protein nodes to set the correspondingARG feild in order to link the two binary trees
// Input: strain of the Genome at the current node of the Protein tree
// Output: A pointer to an ARGNode that either exists in the ARG tree currently or
//
//	indicates the potential for an Ancestor Strain of the ammino acid sequence because differing nucleotide sequences
//	may or may not result in a change in ammino acid sequence and thus may have a common ancestor with another more closely related strain
func FindARGInTree(strain string, args ARGTree) *ARG {
	//iterate through the ARG tee to find the GenomeStrain in the ARG tree
	for i := range args {
		if args[i].ARGpointer.GenomeStrain == strain {
			return args[i].ARGpointer
		}

	}

	//if there is no Label that matches the current strain, then set the lable to ancestor Strain to indicate potential ancestor
	return &ARG{ArgLabel: "Ancestor Strain"}
}

// Charlotte Barretto 11/10/24
// Input: pointer to ProteinNode, didn't make a method because the node itself isn't being updated in any of its feilds
// Output: num of leaves of a ProteinNode
func CountLeavesProtein(protein *ProteinNode) int {

	if protein.Descendants1 == nil && protein.Descendants2 == nil {
		return 1
	}
	//left treee traversal
	if protein.Descendants1 == nil {
		return CountLeavesProtein(protein.Descendants1)
	}
	//right tree traversal
	if protein.Descendants1 == nil {
		return CountLeavesProtein(protein.Descendants2)
	}
	return CountLeavesProtein(protein.Descendants1) + CountLeavesProtein(protein.Descendants2)
}

// Charlotte Barretto
// Charlotte Barretto 11/10/24
// Input: method for the ProteinTree already created
// Output: An array of pointers to ProteinNode objects that represent the clusters in the distance matrix and the leaves of the tree
func (proteinTree ProteinTree) InitializeClustersForProtein() []*ProteinNode {
	//want to get the num nodes n + 1/2 for the row and the col because the row + col nodes is 2n and a tree has 2n +1/ 2 total nodes
	numNodes := len(proteinTree)
	numLeaves := (numNodes + 1) / 2

	clusters := make([]*ProteinNode, numLeaves)

	for i := range clusters {
		clusters[i] = proteinTree[i] // copy pointer to clusters[i]
	}

	return clusters
}

// UPGMA tree construction
// tree of ARG relation, a child of an ARG is the minimum number of mutations away from the parent or ancestor node
// need to update the tree via preorder traversal
// will start with one arg and link that root to another tree?
// Charlotte Barretto (edited by Isabelle)
// Input: Passes *ARG.genome which will be stored in an array of strings
// Output: ARGTree that generalizes the distance between different strains' mutations
func CreateNucleotideARGTree(mtx DistanceMatrix, GenomeStrains, GenomeNames []string, searchString string, ARGName string, percentIdentity float64) ARGTree {
	CheckSquareMatrixAndSameNumSpecies(mtx, GenomeStrains)
	ARGTree := InitializeTree(GenomeStrains, GenomeNames, searchString, ARGName, percentIdentity)
	//the last node will have the greatest percent identity to the search string
	numLeaves := len(mtx)
	clusters := ARGTree.InitializeClustersForARG()

	//traverse internal nodes of ARGTree
	for i := numLeaves; i < 2*numLeaves-1; i++ {
		row, col, _ := MinIdentityDifference(mtx)
		//MinIdentityDifference

		ARGTree[i].Descendents1 = clusters[row]
		ARGTree[i].Descendents2 = clusters[col]

		//size of clusters
		clusters1 := CountLeavesARG(ARGTree[i].Descendents1)
		clusters2 := CountLeavesARG(ARGTree[i].Descendents2)

		mtx = AddRowCol(row, col, clusters1, clusters2, mtx)
		mtx = DeleteRowCol(mtx, row, col)

		clusters = append(clusters, ARGTree[i])
		clusters = DeleteClustersARG(clusters, row, col)

	}
	return ARGTree
}

// Charlotte Barretto
// Input: ARGNode
// Output: the number of leaves for that node
func CountLeavesARG(arg *ARGNode) int {

	if arg.Descendents1 == nil && arg.Descendents2 == nil {
		return 1
	}
	//left treee traversal

	if arg.Descendents1 == nil {
		return CountLeavesARG(arg.Descendents1)
	}
	//right tree traversal
	if arg.Descendents1 == nil {
		return CountLeavesARG(arg.Descendents2)
	}
	return CountLeavesARG(arg.Descendents1) + CountLeavesARG(arg.Descendents2)

}

// Charlotte Barretto 11/9/24
// Input: method for the ARGTree already created
// Output: An array of pointers to ARGNode objects that represent the clusters in the distance matrix and the leaves of the tree
func (ARGTree ARGTree) InitializeClustersForARG() []*ARGNode {
	//want to get the num nodes n + 1/2 for the row and the col because the row + col nodes is 2n and a tree has 2n +1/ 2 total nodes
	numNodes := len(ARGTree)
	numLeaves := (numNodes + 1) / 2

	clusters := make([]*ARGNode, numLeaves)

	for i := range clusters {
		clusters[i] = ARGTree[i] // copy pointer to clusters[i]
	}

	return clusters
}

// Charlotte Barretto 11/9/24 (edited by Isabelle)
// Input: current DistanceMatrix
// Output: the row and col for the percent identity that is the lowest
func MinIdentityDifference(mtx DistanceMatrix) (int, int, float64) {
	row := 0
	col := 1
	minIdentityDifference := mtx[row][col]

	//range over mtx and check to see if there are any values less than the first element
	// because row is set to 0 i := 1 and because we dont want to compare the same strains j := i + 1
	for i := 0; i < len(mtx); i++ {
		for j := i + 1; j < len(mtx[i]); j++ {

			if mtx[i][j] > minIdentityDifference {
				minIdentityDifference = mtx[i][j]
				row = i
				col = j
			}
		}

	}
	return row, col, minIdentityDifference
}

// Charlotte Barretto 11/9/24
// Input: DistanceMatrix and genomeStrains []string
// Output: None, checks to make sure the distance matrix was initialized correctly
func CheckSquareMatrixAndSameNumSpecies(mtx DistanceMatrix, genomeStrains []string) {
	numRows := len(mtx)
	// check that number of elements in each row is equal to numRows
	for r := 0; r < numRows; r++ {
		if len(mtx[r]) != numRows {
			//fmt.Println("Row", r, "of matrix has length", len(mtx[r]), "and matrix has", numRows, "rows.")
			panic("ERROR! ")
		}
	}
	if len(mtx) != len(genomeStrains) {
		panic("Error: Number of rows of matrix don't match number of species.")
	}
}

// Charlotte Barretto 11/9/24 (edited by Isabelle)
// Input: GenomeStrains is an array of genome sequences as strings, GenomeNames is a string denoting the genome from
//
//	the FASTA file, searchString is the nucleotide sequence of the genome being searched, percentIdentity is the similarity between sequences
//
// Output: ARGTree is an array of ARGNodes to represent the nucleotide sequences for each genome
func InitializeTree(GenomeStrains, GenomeNames []string, searchString, ARGName string, percentIdentity float64) ARGTree {
	//create the total number of leaves based on the number of genome strains
	////-> could also be done with the num GenomeNames because this is inputed from the same place
	numLeaves := len(GenomeStrains)

	if numLeaves != len(GenomeNames) {
		panic("Number of genomes does not match number of names. Incorrect input format")
	}

	tree := make([]*ARGNode, 2*numLeaves-1)

	//iterate over all the nodes in the tree, setting the feild of the ARGNode as we traverse initially
	for i := range tree {

		if i < numLeaves {
			ARGIndexInGenome := IdentifyARGs(GenomeStrains[i], searchString, percentIdentity)
			tree[i] = &ARGNode{}
			tree[i].ARGpointer = &ARG{GenomeStrain: GenomeStrains[i],
				Nucleotides: SaveArgRegions(GenomeStrains[i],
					ARGIndexInGenome)[0],
				IndexInGenome: IdentifyARGs(GenomeStrains[i], searchString, percentIdentity)[0],
				ArgLabel:      searchString,
				GenomeName:    GenomeNames[i]}
		} else {
			//create a "dummy node" for the tree output to prevent nil poitnter dereferencing bug
			tree[i] = &ARGNode{}
			tree[i].ARGpointer = &ARG{GenomeStrain: "Ancestor Strain: " + strconv.Itoa(i)}
		}

	}
	return tree
}

// Charlotte Barretto 11/9/24 (edited by Isabelle)
// adds a new row and column that is the arithmatic mean of the items in the two clusters passed
// Input: row, col int and the number of items in cluster1, cluster 2 int, and the current Distance Matrix
// Output: updated DistanceMatrix with added row and column
func AddRowCol(row, col, clusterSize1, clusterSize2 int, mtx DistanceMatrix) DistanceMatrix {
	numRows := len(mtx)
	newRow := make([]float64, numRows+1)

	for r := 0; r < len(newRow)-1; r++ {
		//if at the row and col of the new row, set values as average of cluster1[row][col] and cluster2[row][col]
		if r != row && r != col {
			//average the values for the two clusters for each col of the row being added
			newRow[r] = (float64(clusterSize1)*mtx[r][col] + float64(clusterSize2)*mtx[r][row]) / float64(clusterSize1+clusterSize2)
		}
	}

	mtx = append(mtx, newRow)

	//append values to the col correspoding to the row
	for c := 0; c < numRows; c++ {
		mtx[c] = append(mtx[c], newRow[c])
	}

	return mtx
}

// Charlotte Barretto 11/9/24 (edited by Isabelle)
// Input: an array of pointers to ARGNode clusters, and the row, col int of the cluster we want to delete
// Output: a pointer to the updated clusters
func DeleteClustersARG(clusters []*ARGNode, row, col int) []*ARGNode {
	//delete col from clusters
	clusters = append(clusters[:col], clusters[col+1:]...)

	// delete row from cluster
	clusters = append(clusters[:row], clusters[row+1:]...)

	return clusters
}

// Charlotte Barretto 11/9/24 (edited by Isabelle)
// Input: array of ProteinNode pointers named clusters and the row, col int to be deleted from the protein tree
// Output: updated protein clusters
func DeleteClustersProtein(clusters []*ProteinNode, row, col int) []*ProteinNode {
	//get rid of clusters col
	clusters = append(clusters[:col], clusters[col+1:]...)

	// next, get rid of clusters row
	clusters = append(clusters[:row], clusters[row+1:]...)

	return clusters
}

// Charlotte Barretto 11/9/24
// Input: DistanceMatrix, and the row, col int that is intended to be deleted
// Output: Updated DistanceMatrix with input row and colum deleted
func DeleteRowCol(mtx DistanceMatrix, row, col int) DistanceMatrix {
	//deleting the row
	mtx = append(mtx[:col], mtx[col+1:]...)
	mtx = append(mtx[:row], mtx[row+1:]...)

	//deleteing the colum values for the row being deleted
	for r := range mtx {
		mtx[r] = append(mtx[r][:col], mtx[r][col+1:]...)
		mtx[r] = append(mtx[r][:row], mtx[r][row+1:]...)
	}
	return mtx
}

// CalculatePreservationOfNucleotides (Shruthi)
// INPUT: nothing
// OUTPUT: nothing (it is a method called on an ARGTree in which each node contains a preservation parameter that represents the number of mutations that have been preserved across the wildtype ARG sequence and the ARG sequence associated with a given node)
func (nucleotideTree ARGTree) CalculatePreservationOfNucleotides(wildtypeSeq string) {
	for _, node := range nucleotideTree {
		node.Preservation = PercentNucleotideIdentity(node.ARGpointer.Nucleotides, wildtypeSeq)
	}
}

// CalculatePreservationOfProteins (Shruthi)
// INPUT: nothing
// OUTPUT: nothing (it is a method called on a ProteinTree object in which each node contains a preservation parameter that represents the number of mutations that have been preserved across the wildtype protein sequence and the protein sequence associated with a given node)
func (proteinTree ProteinTree) CalculatePreservationOfProteins(wildtypeSeq string) {
	// Create a slice containing the wildtype nucleotide sequence
	wildtypeSeqSlice := make([]string, 1)
	wildtypeSeqSlice = append(wildtypeSeqSlice, wildtypeSeq)

	// Create a string containing the translated wildtype nucleotide sequence
	wildtypeProteinSeq := TranslateProteinSequences(wildtypeSeqSlice)
	wildtypeProteinSeqString := strings.Join(wildtypeProteinSeq, "")

	// Set the preservation attribute for each node of proteinTree
	for _, node := range proteinTree {
		if node.Protein != nil {
			node.Preservation = PercentProteinIdentity(node.Protein.AminoAcids, wildtypeProteinSeqString)
		}
	}
}

// CalculateLikelihood (Shruthi)
// INPUT: a threshold value set by the user
// OUTPUT: nothing (it is a method that is called on an ARGTree object where each node contains a likelihood attribute that is set equal to 1, which represents increased likelihood of antibiotic resistance, or 0, which represents decreased likelihood of antibiotic resistance)
func (nucleotideTree ARGTree) CalculateLikelihood(threshold float64) {
	for _, node := range nucleotideTree {
		if node.Preservation < threshold {
			node.Likelihood = 1.0 // 1.0 represents increased likelihood of antibiotic resistance
		} else {
			node.Likelihood = 0.0 // 0.0 represents decreased likelihood of antibiotic resistance
		}
	}
}

// idamico
// NewickNucleotide is an ARGNode method, called on the root of a tree.
// It returns the Newick representation of the nucleotide tree.
func (n *ARGNode) NewickNucleotide() string {
	//check to make sure the Node exists
	if n == nil {
		return ""
	}

	//create our string
	var sb strings.Builder

	//if we have both descendents, send those both off recursively
	if n.Descendents1 != nil && n.Descendents2 != nil {
		sb.WriteString("(")
		sb.WriteString(n.Descendents1.NewickNucleotide())
		sb.WriteString(",")
		sb.WriteString(n.Descendents2.NewickNucleotide())
		sb.WriteString(")")
	} else if n.Descendents1 != nil { //if we only have the first one, send that off recursively
		sb.WriteString("(")
		sb.WriteString(n.Descendents1.NewickNucleotide())
		sb.WriteString(")")
	} else if n.Descendents2 != nil { //if we only have the second one, send that off recursively
		sb.WriteString("(")
		sb.WriteString(n.Descendents2.NewickNucleotide())
		sb.WriteString(")")
	}

	//write the name of this genome
	sb.WriteString(n.ARGpointer.GenomeName)

	//add the likelihood of antibiotic resistance and the name of this genome in NHX format
	sb.WriteString("[&&NHX:S=")
	sb.WriteString(strconv.FormatFloat(n.Likelihood, 'f', 3, 64))
	sb.WriteString(":E=")
	sb.WriteString(n.ARGpointer.GenomeName)
	sb.WriteString("]")

	//return our string
	return sb.String()
}

// idamico
// NewickNucleotide is an ARGNode method, called on the root of a tree.
// It returns the Newick representation of the nucleotide tree.
func (n *ProteinNode) NewickProtein() string {
	//check to make sure the Node exists
	if n == nil {
		return ""
	}

	//create our string
	var sb strings.Builder

	//if we have both descendents, send those both off recursively
	if n.Descendants1 != nil && n.Descendants2 != nil {
		sb.WriteString("(")
		sb.WriteString(n.Descendants1.NewickProtein())
		sb.WriteString(",")
		sb.WriteString(n.Descendants2.NewickProtein())
		sb.WriteString(")")
	} else if n.Descendants1 != nil { //if we only have the first descendent, send those off recursively
		sb.WriteString("(")
		sb.WriteString(n.Descendants1.NewickProtein())
		sb.WriteString(")")
	} else if n.Descendants2 != nil { //if we only have the second descendent, send those off recursively
		sb.WriteString("(")
		sb.WriteString(n.Descendants2.NewickProtein())
		sb.WriteString(")")
	}

	//write the name of this genome, if it has one
	if n.Protein != nil {
		sb.WriteString(n.Protein.ARG.GenomeName)
	} else {
		sb.WriteString("")
	}

	//add the preservation of amino acid sequence (if it is preserved) and the name of this genome in NHX format
	if n.Preservation > 0.05 {
		sb.WriteString("[&&NHX:S=PP:E=")
		sb.WriteString(n.Protein.ARG.GenomeName)
		sb.WriteString("]")
	} else {
		sb.WriteString("[&&NHX:E=")
		if n.Protein != nil {
			sb.WriteString(n.Protein.ARG.GenomeName)
		} else {
			sb.WriteString("")
		}
		sb.WriteString("]")
	}

	//return out string
	return sb.String()
}

// idamico
// WriteNewickToFile takes as input a Newick-formatted string and a string to differentiate nucleotide from protein.
// It returns no outputs but writes the string to a file.
func WriteNewickToFile(input, inputType string) {
	//initializing variables used in if blocks below
	var file *os.File
	var err error

	//create and name our files, depending on which type is passed in
	if inputType == "Nucleotide" {
		//Create our output file
		file, err = os.Create("NucleotideNewickOutput.txt")
		if err != nil {
			panic(err)
		}

	} else if inputType == "Protein" {
		//Create our output file
		file, err = os.Create("ProteinNewickOutput.txt")
		if err != nil {
			panic(err)
		}

	} else {
		panic("Was not expecting that type")
	}

	// Ensure the file is closed when we're done
	defer file.Close()

	// Create a buffered writer
	writer := bufio.NewWriter(file)

	// Write the FASTA string for this node to the buffered writer
	_, err = writer.WriteString(input)
	if err != nil {
		panic(err)
	}

	// Flush the buffered writer to ensure all data is written to the file
	err = writer.Flush()
	if err != nil {
		panic(err)
	}

}

// idamico
// WritetoNucleotideFasta is an ARGNode method, called on the root of a tree.
// It takes as input a file.
// It writes the nucleotide sequences for each node in the tree to a file, in FASTA format.
func (n *ARGNode) WritetoNucleotideFasta(file *os.File) {
	//check to make sure the Node exists
	if n == nil {
		panic("there's nothing here at this node.")
	}

	// Create a buffered writer
	writer := bufio.NewWriter(file)

	//if we have both descendents, send off the calls for this node's children recursively
	if n.Descendents1 != nil && n.Descendents2 != nil {
		n.Descendents1.WritetoNucleotideFasta(file)
		n.Descendents2.WritetoNucleotideFasta(file)

	} else if n.Descendents1 != nil { //if we have only the firsr descendent, send off the call for this node's child recursively
		n.Descendents1.WritetoNucleotideFasta(file)

	} else if n.Descendents2 != nil { //if we have only the second descendent, send off the call for this node's child recursively
		n.Descendents2.WritetoNucleotideFasta(file)

	} else { //if we get to the leaf, we write
		// Write the FASTA string for this node to the buffered writer
		_, err := writer.WriteString(">" + n.ARGpointer.GenomeName + "\n" + n.ARGpointer.Nucleotides + "\n")
		if err != nil {
			panic(err)
		}

	}

	// Flush the buffered writer to ensure all data is written to the file
	err := writer.Flush()
	if err != nil {
		panic(err)
	}
}

// idamico
// WritetoProteinFasta is an ProteinNode method, called on the root of a tree.
// It takes as input a file.
// It writes the amino acid sequences for each node in the tree to a file, in FASTA format.
func (n *ProteinNode) WritetoProteinFasta(file *os.File) {
	//check to make sure the Node exists
	if n == nil {
		panic("there's nothing here at this node.")
	}

	// Create a buffered writer
	writer := bufio.NewWriter(file)

	//if we have both second descendent, send off the calls for this node's children recursively
	if n.Descendants1 != nil && n.Descendants2 != nil {
		n.Descendants1.WritetoProteinFasta(file)
		n.Descendants2.WritetoProteinFasta(file)

	} else if n.Descendants1 != nil { //if we have only the first descendent, send off the call for this node's child recursively
		n.Descendants1.WritetoProteinFasta(file)

	} else if n.Descendants2 != nil { //if we have only the second descendent, send off the call for this node's child recursively
		n.Descendants2.WritetoProteinFasta(file)

	} else { //if we get to the leaf, we write
		// Write the FASTA string for this node to the buffered writer
		_, err := writer.WriteString(">" + n.Protein.ARG.GenomeName + "\n" + n.Protein.AminoAcids + "\n")
		if err != nil {
			panic(err)
		}

	}

	// Flush the buffered writer to ensure all data is written to the file
	err := writer.Flush()
	if err != nil {
		panic(err)
	}

}

// idamico
// ReadGenomesFromFile takes as input a string filepath for a FASTA-formatted file.
// It returns two slices of strings, one corresponding to genome nucleotide sequences, and one corresponding to genome names, read from the file.
func ReadGenomesFromFile(filepath string) ([]string, []string) {
	// Open the file
	file, err := os.Open(filepath)
	if err != nil {
		panic(err)
	}

	// Ensure the file is closed after we're done
	defer file.Close()

	//initialize variables
	var currentGenome strings.Builder
	genomes := make([]string, 0)
	genomeNames := make([]string, 0)
	count := 0

	// Create a scanner to read the file line by line
	scanner := bufio.NewScanner(file)

	// Set the maximum token size to 50MB because these genomes can get rather large
	const maxTokenSize = 50 * 1024 * 1024
	buf := make([]byte, 0, 5*1024*1024)
	scanner.Buffer(buf, maxTokenSize)

	//scan through our file
	for scanner.Scan() {
		line := scanner.Text()

		//if we get to a line with the >, we know it's a new genome
		if strings.HasPrefix(line, ">") {
			//save off our genome name and append to genomeName slice
			name := line[1:]
			nameLine := strings.Split(name, " ")
			genomeNames = append(genomeNames, nameLine[0])

			// If we have a current genome, add it to the genome slice
			if currentGenome.String() != "" {
				genomes = append(genomes, currentGenome.String())
				currentGenome.Reset()
			}

			//increment our counter and skip this line
			count++
			continue // Skip the line with '>'

		} else if count == 0 { //if we don't have a first line with >, we know it's not formatted correctly
			panic("Incorrectly formatted file. Please format file as a FASTA file.")
		} else { //otherwise, this line is part of the current genome - we add it to the string
			currentGenome.WriteString(line)

		}
	}

	//check for any errors with scanning, specifically
	if err := scanner.Err(); err != nil {
		fmt.Println("Error scanning:", err)
	}

	// Add the last genome if it exists
	if currentGenome.String() != "" {
		genomes = append(genomes, currentGenome.String())
	}

	return genomes, genomeNames
}

// idamico
// ReadGeneFromFile takes as input a string filepath for a FASTA-formatted file.
// It returns a string of gene sequence, read from the file.
func ReadGeneFromFile(filepath string) (string, string) {
	// Open the file
	file, err := os.Open(filepath)
	if err != nil {
		panic(err)
	}

	// Ensure the file is closed after we're done
	defer file.Close()

	//initialize variables
	var gene string

	// Create a scanner to read the file
	scanner := bufio.NewScanner(file)

	//read in the gene name from the file
	scanner.Scan()
	geneName := scanner.Text()

	//if the line doesn't have a > in front of it, we know it's incorrectly formatted
	//otherwise, we save off the gene name
	if strings.HasPrefix(geneName, ">") {
		geneName = geneName[1:]
		nameLine := strings.Split(geneName, " ")
		geneName = nameLine[0]
	} else {
		panic("Incorrectly formatted file. Please format file as a FASTA file.")
	}

	//read in the gene from the file
	scanner.Scan()
	gene = scanner.Text()

	//if we didn't pick up the gene, something went wrong
	if gene == "" {
		panic("Could not read in the gene from file.")
	}

	//return our gene information
	return gene, geneName
}

// idamico
// WriteLikelihoodtoCSV is an ARGNode method.
// It traverses the tree, starting at the root, and prints out each node's genome name and likelihood to a csv file.
func (n *ARGNode) WriteLikelihoodtoCSV(writer *csv.Writer) {
	//if this node is blank, we back out
	if n == nil {
		return
	}

	if n.Descendents1 == nil && n.Descendents2 == nil {
		// Write the current node's label and likelihood to the CSV file
		err := writer.Write([]string{n.ARGpointer.GenomeName, fmt.Sprintf("%f", n.Likelihood)})
		if err != nil {
			panic(err)
		}
	}

	// Recursively traverse the left and right children, if this node has any
	n.Descendents1.WriteLikelihoodtoCSV(writer)
	n.Descendents2.WriteLikelihoodtoCSV(writer)
}

// idamico
// WriteTreeToCSV is an ARGNode method.
// It writes the likelihood and genome name for each node in the tree to a csv file.
func (root *ARGNode) WriteTreeToCSV(filename string) {
	// Create a CSV file
	file, err := os.Create(filename)
	if err != nil {
		panic(err)
	}

	//defer file closing until we are done
	defer file.Close()

	//create a writer to write to CSV file
	writer := csv.NewWriter(file)

	//wait to flush until we are done
	defer writer.Flush()

	// Traverse the tree and write to the CSV file
	root.WriteLikelihoodtoCSV(writer)
}

// idamico
// WritePreservationtoCSV is an ARGNode method.
// It traverses the tree, starting at the root, and prints out each node's genome name and preservation to a csv file.
func (n *ProteinNode) WritePreservationtoCSV(writer *csv.Writer) {
	//if this node is blank, we back out
	if n == nil {
		return
	}

	//if this is a leaf node, we write the information
	if n.Descendants1 == nil && n.Descendants2 == nil {
		// Write the current node's label and preservation to the CSV file
		err := writer.Write([]string{n.Protein.ARG.GenomeName, fmt.Sprintf("%f", n.Preservation)})

		if err != nil {
			panic(err)
		}
	}

	// Recursively traverse the left and right children, if this node has any
	n.Descendants1.WritePreservationtoCSV(writer)
	n.Descendants2.WritePreservationtoCSV(writer)
}

// idamico
// WriteTreeToCSV is an ProteinNode method.
// It writes the preservation and genome name for each node in the tree to a csv file.
func (root *ProteinNode) WriteTreeToCSV(filename string) {
	// Create a CSV file
	file, err := os.Create(filename)

	if err != nil {
		panic(err)
	}

	//defer the file closing until we are done
	defer file.Close()

	//create a writer to write to CSV file
	writer := csv.NewWriter(file)

	//defer flushing the file until we are done
	defer writer.Flush()

	// Traverse the tree and write to the CSV file
	root.WritePreservationtoCSV(writer)
}
