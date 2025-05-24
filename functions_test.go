package main

import (
	"bufio"
	"fmt"
	"io/fs"
	"math"
	"os"
	"path/filepath"
	"reflect"
	"strconv"
	"strings"
	"testing"
)

// idamico
// PercentIdentityTest is a struct that holds the information for a test that takes two strings.
type PercentIdentityTest struct {
	str1, str2 string
	result     float64
}

// idamico
// DistanceMatrixTest is a struct that holds the information for a test that takes a slice of strings.
type DistanceMatrixTest struct {
	strings []string
	result  [][]float64
}

// idamico
// InitializeMatrixTest is a struct that holds the information for a test that takes an integer.
type InitializeMatrixTest struct {
	rowCols int
	result  [][]float64
}

// idamico
// IdentifyARGsTest is a struct that holds the information for a test that takes two strings and a float.
type IdentifyARGsTest struct {
	genome, searchGene string
	percentIdentity    float64
	result             [][2]int
}

// idamico
// IdentifyARGsForGenomesTest is a struct that holds the information for a test that takes a slice of strings, a string, and a float.
type IdentifyARGsForGenomesTest struct {
	genomes             []string
	searchGene          string
	percentIndentity    float64
	resultARGregions    []string
	resultGeneIndices   [][2]int
	resultGenomesRemove []int
}

// idamico
// RemoveGenomesTest is a struct that holds the information for a test that takes two slices of strings and a slice of integers.
type RemoveGenomesTest struct {
	genomes, names             []string
	indices                    []int
	resultGenomes, resultNames []string
}

// idamico
// NewickFormatTest is a struct that holds the information for a test that takes a tree.
type NewickFormatTest struct {
	proteinTreeRoot                 *ProteinNode
	nucleotideTreeRoot              *ARGNode
	proteinResult, nucleotideResult string
}

// TranslateProteinSequencesTest is a struct that holds the information for a test that takes a slice of strings (Shruthi)
type TranslateProteinSequencesTest struct {
	dnaFragments []string
	result       []string
}

// RNAToAminoAcidTest is a struct that holds the information for a test that takes a string (Shruthi)
type RNAtoAminoAcidTest struct {
	codon  string
	result string
}

// PreservationOfNucleotidesTest is a struct that holds the information for a test that takes an ARGTree (Shruthi)
type PreservationOfNucleotidesTest struct {
	wildtypeSeq          string
	nucleotideTree       ARGTree
	expectedPreservation []float64
}

// PreservationOfProteinsTest is a struct that holds the information for a test that takes a ProteinTree (Shruthi)
type PreservationOfProteinsTest struct {
	wildtypeSeq          string
	proteinTree          ProteinTree
	expectedPreservation []float64
}

// CalculateLikelihoodTest is a struct that holds the information for a test that takes an ARGTree and a float64 (Shruthi)
type CalculateLikelihoodTest struct {
	nucleotideTree     ARGTree
	threshold          float64
	expectedLikelihood []float64
}

// idamico
// ReadDirectory reads in a directory and returns a slice of fs.DirEntry objects containing file info for directory
func ReadDirectory(dir string) []fs.DirEntry {
	//read in all files in the given directory
	files, err := os.ReadDir(dir)
	if err != nil {
		panic(err)
	}
	return files
}

// idamico
// ReadIdentifyARGsFromFile reads in two strings and a float from a file.
func ReadIdentifyARGsFromFile(filename string) (string, string, float64) {
	// Open the file
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	// Create variables
	var genome, searchGene string
	var percentIdentity float64
	count := 1

	// Create a scanner to read the file line by line
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		// Read the line and split it into strings
		line := scanner.Text()

		if count == 1 {
			genome = line
		} else if count == 2 {
			searchGene = line
		} else if count == 3 {
			percentIdentity, err = strconv.ParseFloat(line, 64)
			if err != nil {
				panic(err)
			}
		}

		count++
	}

	return genome, searchGene, percentIdentity
}

// idamico
// ReadNucleotideTreeFile takes as input a filepath string. It reads the tree information from this file and outputs the root ARGnode from the tree.
func ReadNucleotideTreeFile(filepath string) *ARGNode {
	//Open the file
	file, err := os.Open(filepath)
	if err != nil {
		panic(err)
	}

	//Defer the file closing
	defer file.Close()

	//initialize variables
	relationships := make(map[int][]int) // Parent -> Children mapping
	nodes := make(map[int]*ARGNode)      // Node ID -> Node object
	argNodes := make(map[int]*ARG)       // ARGNode ID -> ARGNode object

	scanner := bufio.NewScanner(file)
	isRelationshipSection := true
	count := 0
	var rootNodeID int

	// Read the file line by line
	for scanner.Scan() {
		line := scanner.Text()

		// Handle the relationships section
		if isRelationshipSection {
			if strings.TrimSpace(line) == "" {
				// Transition to the attributes section once an empty line is found
				isRelationshipSection = false
				continue
			}

			// Parse the Node -> Node relationship
			parts := strings.Fields(line)
			if len(parts) == 2 {
				parentID, err := strconv.Atoi(parts[0])
				if err != nil {
					panic(err)
				}
				childID, err := strconv.Atoi(parts[1])
				if err != nil {
					panic(err)
				}
				relationships[parentID] = append(relationships[parentID], childID)
				if count == 0 {
					rootNodeID = parentID
				}
				count++
			}
		} else {
			// Parse the attributes section
			parts := strings.Fields(line)
			if len(parts) != 4 {
				panic("invalid line format")
			}

			nodeID, err := strconv.Atoi(parts[0])
			if err != nil {
				panic(err)
			}
			likelihood, err := strconv.ParseFloat(parts[1], 64)
			if err != nil {
				panic(err)
			}
			argNodeID, err := strconv.Atoi(parts[2])
			if err != nil {
				panic(err)
			}
			genomeName := parts[3]

			// Create the ARGNode if it doesn't exist
			if _, exists := argNodes[argNodeID]; !exists {
				argNodes[argNodeID] = &ARG{GenomeName: genomeName}
			}

			// Create the Node object
			node := ARGNode{
				Likelihood: likelihood,
				ARGpointer: argNodes[argNodeID],
			}
			nodes[nodeID] = &node
		}
	}

	if err := scanner.Err(); err != nil {
		panic(err)
	}

	// Now that we have all the nodes and relationships, link the Descendent1 and Descendent2 fields
	for parentID, children := range relationships {
		parent, ok := nodes[parentID]
		if !ok {
			panic("parent node not found")
		}

		if len(children) > 0 {
			childID := children[0]
			child, ok := nodes[childID]
			if !ok {
				panic("child node not found")
			}
			parent.Descendents1 = child
		}
		if len(children) > 1 {
			childID := children[1]
			child, ok := nodes[childID]
			if !ok {
				panic("child node not found")
			}
			parent.Descendents2 = child
		}
	}

	// Return the root node (typically the first node or the one without a parent)
	// Assuming the root is the node with no parent in the relationships map
	var rootNode *ARGNode
	rootNode = nodes[rootNodeID]

	return rootNode
}

// idamico
// ReadProteinTreeFile takes as input a filepath string. It reads the tree information from this file and outputs the root ProteinNode from the tree.
func ReadProteinTreeFile(filepath string) *ProteinNode {
	//Open the file
	file, err := os.Open(filepath)
	if err != nil {
		panic(err)
	}

	//Defer the file closing
	defer file.Close()

	//initialize variables
	relationships := make(map[int][]int) // Parent -> Children mapping
	nodes := make(map[int]*ProteinNode)  // Node ID -> Node object
	argNodes := make(map[int]*ARG)       // ARGNode ID -> ARGNode object

	scanner := bufio.NewScanner(file)
	isRelationshipSection := true
	count := 0
	var rootNodeID int

	// Read the file line by line
	for scanner.Scan() {
		line := scanner.Text()

		// Handle the relationships section
		if isRelationshipSection {
			if strings.TrimSpace(line) == "" {
				// Transition to the attributes section once an empty line is found
				isRelationshipSection = false
				continue
			}

			// Parse the Node -> Node relationship
			parts := strings.Fields(line)
			if len(parts) == 2 {
				parentID, err := strconv.Atoi(parts[0])
				if err != nil {
					panic(err)
				}
				childID, err := strconv.Atoi(parts[1])
				if err != nil {
					panic(err)
				}
				relationships[parentID] = append(relationships[parentID], childID)
				if count == 0 {
					rootNodeID = parentID
				}
				count++
			}
		} else {
			// Parse the attributes section
			parts := strings.Fields(line)
			if len(parts) != 4 {
				panic("invalid line format")
			}

			nodeID, err := strconv.Atoi(parts[0])
			if err != nil {
				panic(err)
			}
			preservation, err := strconv.ParseFloat(parts[1], 64)
			if err != nil {
				panic(err)
			}
			argNodeID, err := strconv.Atoi(parts[2])
			if err != nil {
				panic(err)
			}
			genomeName := parts[3]

			// Create the ARGNode if it doesn't exist
			if _, exists := argNodes[argNodeID]; !exists {
				argNodes[argNodeID] = &ARG{GenomeName: genomeName}
			}

			// Create the Node object
			node := ProteinNode{
				Preservation: preservation,
				Protein:      &Protein{ARG: argNodes[argNodeID]},
			}
			nodes[nodeID] = &node
		}
	}

	if err := scanner.Err(); err != nil {
		panic(err)
	}

	// Now that we have all the nodes and relationships, link the Descendent1 and Descendent2 fields
	for parentID, children := range relationships {
		parent, ok := nodes[parentID]
		if !ok {
			panic("parent node not found")
		}

		if len(children) > 0 {
			childID := children[0]
			child, ok := nodes[childID]
			if !ok {
				panic("child node not found")
			}
			parent.Descendants1 = child
		}
		if len(children) > 1 {
			childID := children[1]
			child, ok := nodes[childID]
			if !ok {
				panic("child node not found")
			}
			parent.Descendants2 = child
		}
	}

	// Return the root node (typically the first node or the one without a parent)
	// Assuming the root is the node with no parent in the relationships map
	var rootNode *ProteinNode
	rootNode = nodes[rootNodeID]

	return rootNode
}

// idamico
// ReadIdentifyARGsForGenomesFromFile reads in a slice of strings, a string, and a float from a file.
func ReadIdentifyARGsForGenomesFromFile(filename string) ([]string, string, float64) {
	// Open the file
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	// Create variables
	genomes := make([]string, 0)
	var searchGene string
	var currentGenome string
	var percentIdentity float64
	count := 1

	// Create a scanner to read the file line by line
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		// Read the line and split it into strings
		line := scanner.Text()

		if count == 1 {
			searchGene = line
		} else if count == 2 {
			percentIdentity, err = strconv.ParseFloat(line, 64)
			if err != nil {
				panic(err)
			}
		} else if count >= 3 {
			if strings.HasPrefix(line, ">") {
				// If we have a current genome, add it to the genome slice
				if currentGenome != "" {
					genomes = append(genomes, currentGenome)
					currentGenome = ""
				}
				continue // Skip the line with '>'
			}

			// Process the genome
			currentGenome = line
		}

		count++
	}

	// Add the last genome if it exists
	if currentGenome != "" {
		genomes = append(genomes, currentGenome)
	}

	return genomes, searchGene, percentIdentity
}

// idamico
// ReadRemoveGenomeSample reads in two slices of strings and a slice of integers from a file.
func ReadRemoveGenomeSample(filename string) ([]string, []string, []int) {
	// Open the file
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	// Create variables
	genomes := make([]string, 0)
	var currentGenome string
	names := make([]string, 0)
	indices := make([]int, 0)
	count := 1

	// Create a scanner to read the file line by line
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		// Read the line and split it into strings
		line := scanner.Text()

		if count == 1 { //grab our indices
			fields := strings.Fields(line)
			for _, val := range fields {
				value, err := strconv.Atoi(val)
				if err != nil {
					panic(err)
				}
				indices = append(indices, value)
			}
		} else if count == 2 { //grab our names
			fields := strings.Fields(line)
			for _, val := range fields {
				names = append(names, val)
			}
		} else if count >= 3 {
			if strings.HasPrefix(line, ">") {
				// If we have a current genome, add it to the genome slice
				if currentGenome != "" {
					genomes = append(genomes, currentGenome)
					currentGenome = ""
				}
				continue // Skip the line with '>'
			}

			// Process the genome
			currentGenome = line
		}

		count++
	}

	// Add the last genome if it exists
	if currentGenome != "" {
		genomes = append(genomes, currentGenome)
	}

	return genomes, names, indices
}

// idamico
// ReadFloatFromFile reads in a single float from a file
func ReadFloatFromFile(file string) float64 {
	//open the file
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	//create a new scanner
	scanner := bufio.NewScanner(f)

	//read in the line
	scanner.Scan()
	line := scanner.Text()

	//convert the line to an int using strconv
	value, err := strconv.ParseFloat(line, 64)
	if err != nil {
		panic(err)
	}

	return value
}

// idamico
// ReadStringFromFile reads in a single string from a file
func ReadStringFromFile(file string) string {
	//open the file
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	//create a new scanner
	scanner := bufio.NewScanner(f)

	//read in the line
	scanner.Scan()
	line := scanner.Text()

	return line
}

// idamico
// ReadIntFromFile reads in a single integer from a file
func ReadIntFromFile(file string) int {
	//open the file
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	//create a new scanner
	scanner := bufio.NewScanner(f)

	//read in the line
	scanner.Scan()
	line := scanner.Text()

	//convert the line to an int using strconv
	value, err := strconv.Atoi(line)
	if err != nil {
		panic(err)
	}

	return value
}

// idamico
// ReadStringSliceSample reads in a slice of strings.
func ReadStringSliceSample(filename string) []string {
	// Open the file
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	// Create a slice to hold the string values
	strings := make([]string, 0)

	// Create a scanner to read the file line by line
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		// Read the line and split it into strings
		line := scanner.Text()
		strings = append(strings, line)
	}

	// Check for errors during scanning
	if err := scanner.Err(); err != nil {
		panic(err)
	}

	return strings
}

// idamico
// ReadIntSliceSample reads in a slice of integers.
func ReadIntSliceSample(filename string) []int {
	// Open the file
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	// Create a slice to hold the integer values
	var ints []int

	// Create a scanner to read the file line by line
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		// Read the line and split it into strings
		line := scanner.Text()
		fields := strings.Fields(line)

		for _, val := range fields {
			//convert each value to an integer
			value, err := strconv.Atoi(val)

			if err != nil {
				panic(err)
			}

			//append the integer to the slice
			ints = append(ints, value)
		}
	}

	// Check for errors during scanning
	if err := scanner.Err(); err != nil {
		panic(err)
	}

	return ints
}

// idamico
// ReadMatrixSample reads in a matrix from a file.
func ReadMatrixSample(filename string) [][]float64 {
	// Open the file
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	// Create a matrix to hold the float64 values
	var matrix [][]float64

	// Create a scanner to read the file line by line
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		// Read the line and split it into strings
		line := scanner.Text()
		strValues := strings.Fields(line)

		// Create a slice to hold the float64 values for the current row
		var row []float64
		for _, str := range strValues {
			// Convert each string to float64 and append to the row
			value, err := strconv.ParseFloat(str, 64)
			if err != nil {
				panic(err) // Return error if conversion fails
			}
			row = append(row, value)
		}

		// Append the row to the matrix
		matrix = append(matrix, row)
	}

	// Check for errors during scanning
	if err := scanner.Err(); err != nil {
		panic(err)
	}

	return matrix
}

// idamico
// ReadMatrixSample reads in a matrix from a file.
func ReadIntMatrixSample(filename string) [][2]int {
	// Open the file
	file, err := os.Open(filename)
	if err != nil {
		panic(err)
	}
	defer file.Close()

	// Create a matrix to hold the float64 values
	var matrix [][2]int

	// Create a scanner to read the file line by line
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		// Read the line and split it into strings
		line := scanner.Text()
		strValues := strings.Fields(line)

		// Create a slice to hold the float64 values for the current row
		var row [2]int
		for i, str := range strValues {
			// Convert each string to float64 and append to the row
			value, err := strconv.Atoi(str)
			if err != nil {
				panic(err) // Return error if conversion fails
			}
			row[i] = value
		}

		// Append the row to the matrix
		matrix = append(matrix, row)
	}

	// Check for errors during scanning
	if err := scanner.Err(); err != nil {
		panic(err)
	}

	return matrix
}

// idamico
// ReadDistanceMatrixSample reads in a slice of strings from a file
func ReadDistanceMatrixSample(file string) []string {
	//open the file
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	//create a new scanner
	scanner := bufio.NewScanner(f)

	strings := make([]string, 0)

	//while scanner still has lines to read, read in
	for scanner.Scan() {
		//read in the line
		line := scanner.Text()

		//append each string to the slice
		strings = append(strings, line)
	}

	return strings
}

// idamico
// ReadPercentIdentitySample reads in two strings from a file
func ReadPercentIdentitySample(file string) (string, string) {
	//open the file
	f, err := os.Open(file)
	if err != nil {
		panic(err)
	}
	defer f.Close()

	//create a new scanner
	scanner := bufio.NewScanner(f)

	var str1, str2 string
	count := 1

	//while scanner still has lines to read, read in
	for scanner.Scan() {
		//read in the line
		line := scanner.Text()

		//read in first string
		if count == 1 {
			str1 = line

		} else if count == 2 { //read in second string
			str2 = line
		}

		count++
	}

	return str1, str2
}

// idamico
// ReadIdentifyARGsTests takes as input a directory and returns a slice of IdentifyARGsTest objects.
func ReadIdentifyARGsTests(directory string) []IdentifyARGsTest {
	//read in all tests from the directory and run them
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]IdentifyARGsTest, numFiles)
	for i, inputFile := range inputFiles {
		//read in the test's Boids
		tests[i].genome, tests[i].searchGene, tests[i].percentIdentity = ReadIdentifyARGsFromFile(directory + "input/" + inputFile.Name())
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match.")
	}

	for i, outputFile := range outputFiles {
		//read in the test's result
		tests[i].result = ReadIntMatrixSample(directory + "output/" + outputFile.Name())
	}

	return tests
}

// idamico
// ReadIdentifyARGsForGenomesTests takes as input a directory and returns a slice of IdentifyARGsForGenomesTest objects.
func ReadIdentifyARGsForGenomesTests(directory string) []IdentifyARGsForGenomesTest {
	//read in all tests from the directory and run them
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]IdentifyARGsForGenomesTest, numFiles)
	for i, inputFile := range inputFiles {
		//read in the test's Boids
		tests[i].genomes, tests[i].searchGene, tests[i].percentIndentity = ReadIdentifyARGsForGenomesFromFile(directory + "input/" + inputFile.Name())
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")
	if len(outputFiles)/3 != numFiles {
		panic("Error: number of input and output files do not match.")
	}

	for i, outputFile := range outputFiles {
		//read in the test's result
		if i%3 == 0 {
			tests[i/3].resultARGregions = ReadStringSliceSample(directory + "output/" + outputFile.Name())
		} else if i%3 == 1 {
			tests[i/3].resultGeneIndices = ReadIntMatrixSample(directory + "output/" + outputFile.Name())
		} else if i%3 == 2 {
			tests[i/3].resultGenomesRemove = ReadIntSliceSample(directory + "output/" + outputFile.Name())
		}
	}

	return tests
}

// idamico
// ReadRemoveGenomeTest takes as input a directory and returns a slice of RemoveGenomesTest objects
func ReadRemoveGenomeTest(directory string) []RemoveGenomesTest {
	//read in all tests from the directory and run them
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]RemoveGenomesTest, numFiles)
	for i, inputFile := range inputFiles {
		//read in the test's rows/columns
		tests[i].genomes, tests[i].names, tests[i].indices = ReadRemoveGenomeSample(directory + "input/" + inputFile.Name())
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")
	if len(outputFiles)/2 != numFiles {
		panic("Error: number of input and output files do not match.")
	}

	for i, outputFile := range outputFiles {
		//read in the test's result
		if i%2 == 0 {
			tests[i/2].resultGenomes = ReadStringSliceSample(directory + "output/" + outputFile.Name())
		} else if i%2 == 1 {
			tests[i/2].resultNames = ReadStringSliceSample(directory + "output/" + outputFile.Name())
		}

	}

	return tests
}

// idamico
// ReadInitializeMatrixTests takes as input a directory and returns a slice of InitializeMatrixTest objects
func ReadInitializeMatrixTests(directory string) []InitializeMatrixTest {
	//read in all tests from the directory and run them
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]InitializeMatrixTest, numFiles)
	for i, inputFile := range inputFiles {
		//read in the test's rows/columns
		tests[i].rowCols = ReadIntFromFile(directory + "input/" + inputFile.Name())
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match.")
	}

	for i, outputFile := range outputFiles {
		//read in the test's result
		tests[i].result = ReadMatrixSample(directory + "output/" + outputFile.Name())
	}

	return tests
}

// idamico
// ReadDistanceMatrixTests takes as input a directory and returns a slice of DistanceMatrixTest objects
func ReadDistanceMatrixTests(directory string) []DistanceMatrixTest {
	//read in all tests from the directory and run them
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]DistanceMatrixTest, numFiles)
	for i, inputFile := range inputFiles {
		//read in the test's rows/columns
		tests[i].strings = ReadDistanceMatrixSample(directory + "input/" + inputFile.Name())
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match.")
	}

	for i, outputFile := range outputFiles {
		//read in the test's result
		tests[i].result = ReadMatrixSample(directory + "output/" + outputFile.Name())
	}

	return tests
}

// idamico
// ReadNewickTests takes as input a directory and returns a slice of NewickFormatTests
func ReadNewickTests(directory string) []NewickFormatTest {
	//read in all tests from the directory and run them
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]NewickFormatTest, numFiles)
	for i, inputFile := range inputFiles {
		//read in the test's roots
		if i%2 == 0 {
			tests[i].nucleotideTreeRoot = ReadNucleotideTreeFile(directory + "input/" + inputFile.Name())
		} else if i%2 == 1 {
			tests[i].proteinTreeRoot = ReadProteinTreeFile(directory + "input/" + inputFile.Name())
		}
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match.")
	}

	for i, outputFile := range outputFiles {
		//read in the test's result
		if i%2 == 0 {
			tests[i].nucleotideResult = ReadStringFromFile(directory + "output/" + outputFile.Name())
		} else if i%2 == 1 {
			tests[i].proteinResult = ReadStringFromFile(directory + "output/" + outputFile.Name())
		}
	}

	return tests
}

// idamico
// ReadPercentIdentityTests takes as input a directory and returns a slice of PercentIdentityTest objects
func ReadPercentIdentityTests(directory string) []PercentIdentityTest {
	//read in all tests from the directory and run them
	inputFiles := ReadDirectory(directory + "/input")
	numFiles := len(inputFiles)

	tests := make([]PercentIdentityTest, numFiles)
	for i, inputFile := range inputFiles {
		//read in the test's Boids
		tests[i].str1, tests[i].str2 = ReadPercentIdentitySample(directory + "input/" + inputFile.Name())
	}

	//now, read output files
	outputFiles := ReadDirectory(directory + "/output")
	if len(outputFiles) != numFiles {
		panic("Error: number of input and output files do not match.")
	}

	for i, outputFile := range outputFiles {
		//read in the test's result
		tests[i].result = ReadFloatFromFile(directory + "output/" + outputFile.Name())
	}

	return tests
}

// idamico
// TestPercentIdentity tests the NucleotidePercentIdentity and ProteinPercentIdentity functions
func TestPercentIdentity(t *testing.T) {
	//read in all tests from the tests/PercentNucleotideIdentity directory and run them
	tests := ReadPercentIdentityTests("tests/PercentNucleotideIdentity/")
	for _, test := range tests {
		//run the test
		result := PercentNucleotideIdentity(test.str1, test.str2)
		//check the result
		if result != test.result {
			t.Errorf("PercentNucleotideIdentity(%v,%v) %v, want %v", test.str1, test.str2, result, test.result)
		}
	}

	//read in all tests from the tests/PercentNucleotideIdentity directory and run them
	tests2 := ReadPercentIdentityTests("tests/PercentProteinIdentity/")
	for _, test := range tests2 {
		//run the test
		result := PercentProteinIdentity(test.str1, test.str2)
		//check the result
		if result != test.result {
			t.Errorf("PercentProteinIdentity(%v,%v) %v, want %v", test.str1, test.str2, result, test.result)
		}
	}
}

// idamico
// TestDistanceMatrix tests the NucleotideDistanceMatrix, ProteinDistanceMatrix, and InitializeMatrix functions
func TestDistanceMatrix(t *testing.T) {
	//read in all tests from the tests/NucleotideDistanceMatrix directory and run them
	tests := ReadDistanceMatrixTests("tests/NucleotideDistanceMatrix/")
	for _, test := range tests {
		//run the test
		result := NucleotideDistanceMatrix(test.strings)
		//check the result
		if !CheckTwoMatrices(result, test.result) {
			t.Errorf("NucleotideDistanceMatrix(%v) %v, want %v", test.strings, result, test.result)
		}
	}

	//read in all tests from the tests/ProteinDistanceMatrix directory and run them
	tests2 := ReadDistanceMatrixTests("tests/ProteinDistanceMatrix/")
	for _, test := range tests2 {
		//run the test
		result := ProteinDistanceMatrix(test.strings)
		//check the result
		if !CheckTwoMatrices(result, test.result) {
			t.Errorf("ProteinDistanceMatrix(%v) %v, want %v", test.strings, result, test.result)
		}
	}

	//read in all tests from the tests/InitializeSquareMatrix directory and run them
	tests3 := ReadInitializeMatrixTests("tests/InitializeSquareMatrix/")
	for _, test := range tests3 {
		//run the test
		result := InitializeSquareMatrix(test.rowCols)
		//check the result
		if !CheckTwoMatrices(result, test.result) {
			t.Errorf("InitializeSquareMatrix(%v) %v, want %v", test.rowCols, result, test.result)
		}
	}
}

// idamico
// TestIdentifyARGS tests the IdentifyARGs and IdentifyARGsForGenomes functions.
func TestIdentifyARGS(t *testing.T) {
	//read in all tests from the tests/IdentifyARGs directory and run them
	tests := ReadIdentifyARGsTests("tests/IdentifyARGs/")
	for _, test := range tests {
		//run the test
		result := IdentifyARGs(test.genome, test.searchGene, test.percentIdentity)
		//check the result
		if !CheckTwoIntMatrices(result, test.result) {
			t.Errorf("IdentifyARGs(%v,%v,%v) %v, want %v", test.genome, test.searchGene, test.percentIdentity, result, test.result)
		}
	}

	//read in all tests from the tests/IdentifyARGsForGenomes directory and run them
	tests2 := ReadIdentifyARGsForGenomesTests("tests/IdentifyARGsForGenomes/")
	for _, test := range tests2 {
		//run the test
		resultARGregions, resultGeneIndices, resultGenomesRemove := IdentifyARGsForGenomes(test.genomes, test.searchGene, test.percentIndentity)
		//check the result
		if !CheckTwoStringMatrices(resultARGregions, test.resultARGregions) || !CheckTwoIntMatrices(resultGeneIndices, test.resultGeneIndices) || !CheckTwoIntSlices(resultGenomesRemove, test.resultGenomesRemove) {
			t.Errorf("IdentifyARGsForGenomes(%v,%v) %v %v %v, want %v %v %v", test.searchGene, test.percentIndentity, resultARGregions, resultGeneIndices, resultGenomesRemove, test.resultARGregions, test.resultGeneIndices, test.resultGenomesRemove)
		}
	}

	//read in all tests from the tests/IdentifyARGsForGenomes directory and run them
	tests3 := ReadRemoveGenomeTest("tests/RemoveGenomes/")
	for _, test := range tests3 {
		//run the test
		_, genomeNames := RemoveGenomes(test.genomes, test.names, test.indices)
		//check the result
		if !CheckTwoStringSlices(genomeNames, test.resultNames) { //||  {
			t.Errorf("RemoveGenomes(%v,%v) %v, want %v", test.names, test.indices, genomeNames, test.resultNames)
		}
	}
}

// idamico
// TestNewick tests the NewickNucleotide and NewickProtein functions.
func TestNewick(t *testing.T) {
	//read in all tests from the tests/Newick directory and run them
	tests := ReadNewickTests("tests/Newick/")
	for i, test := range tests {
		//run the test
		if i%2 == 0 {
			result := test.nucleotideTreeRoot.NewickNucleotide()
			//check the result
			if result != test.nucleotideResult {
				t.Errorf("%v.NewickNucleotide() %v, want %v", test.nucleotideTreeRoot, result, test.nucleotideResult)
			}
		} else if i%2 == 1 {
			result := test.proteinTreeRoot.NewickProtein()
			//check the result
			if result != test.proteinResult {
				t.Errorf("%v.NewickProtein() %v, want %v", test.proteinTreeRoot, result, test.proteinResult)
			}
		}
	}
}

// idamico
// CheckTwoStringSlices takes two slices as input
// It returns true if all the values are the same, false if there is a difference in the matrix values.
func CheckTwoStringSlices(mtx1, mtx2 []string) bool {
	//enforce the same dimensions for the two matrices
	if len(mtx1) != len(mtx2) {
		panic("The dimensions of these two matrices are not the same.")
	}

	//range over the rows and columns of the matrices and return false if any of them do not match
	for i := range mtx1 {
		if mtx1[i] != mtx2[i] {
			fmt.Println("discrepancy at", i)
			return false
		}
	}

	//if we get here, all the values are the same
	return true
}

// idamico
// CheckTwoMatrices takes two matrices as input
// It returns true if all the values are the same, false if there is a difference in the matrix values.
func CheckTwoMatrices(mtx1, mtx2 [][]float64) bool {
	//enforce the same dimensions for the two matrices
	if len(mtx1) != len(mtx2) || len(mtx1[0]) != len(mtx2[0]) {
		panic("The dimensions of these two matrices are not the same.")
	}

	//range over the rows and columns of the matrices and return false if any of them do not match
	for i := range mtx1 {
		for j := range mtx1[i] {
			if mtx1[i][j] != mtx2[i][j] {
				return false
			}
		}
	}

	//if we get here, all the values are the same
	return true
}

// idamico
// CheckTwoIntMatrices takes two matrices as input
// It returns true if all the values are the same, false if there is a difference in the matrix values.
func CheckTwoIntMatrices(mtx1, mtx2 [][2]int) bool {
	//enforce the same dimensions for the two matrices
	if len(mtx1) != len(mtx2) {
		panic("The dimensions of these two matrices are not the same.")
	}

	//range over the rows and columns of the matrices and return false if any of them do not match
	for i := range mtx1 {
		for j := range mtx1[i] {
			if mtx1[i][j] != mtx2[i][j] {
				return false
			}
		}
	}

	//if we get here, all the values are the same
	return true
}

// idamico
// CheckTwoIntSlices takes two matrices as input
// It returns true if all the values are the same, false if there is a difference in the matrix values.
func CheckTwoIntSlices(mtx1, mtx2 []int) bool {
	//enforce the same dimensions for the two matrices
	if len(mtx1) != len(mtx2) {
		panic("The dimensions of these two matrices are not the same.")
	}

	//range over the rows and columns of the matrices and return false if any of them do not match
	for i := range mtx1 {
		if mtx1[i] != mtx2[i] {
			return false
		}
	}

	//if we get here, all the values are the same
	return true
}

// idamico
// CheckTwoStringMatrices takes two matrices as input
// It returns true if all the values are the same, false if there is a difference in the matrix values.
func CheckTwoStringMatrices(mtx1, mtx2 []string) bool {
	//enforce the same dimensions for the two matrices
	if len(mtx1) != len(mtx2) {
		panic("The dimensions of these two matrices are not the same.")
	}

	//range over the rows and columns of the matrices and return false if any of them do not match
	for i := range mtx1 {
		for j := range mtx1[i] {
			if mtx1[i][j] != mtx2[i][j] {
				return false
			}
		}
	}

	//if we get here, all the values are the same
	return true
}

// Charlotte Barretto
func TestCreateNucleotideARGTree(t *testing.T) {
	inputDir := "Tests/ARGTree/Input"
	outputDir := "Tests/ARGTree/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			matrix, genomeStrains, searchString, argName, err := readARGInputFile(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			percentIdentity := 0.95

			// Create ARG tree
			resultTree := CreateNucleotideARGTree(matrix, genomeStrains, genomeStrains, searchString, argName, percentIdentity)

			// Validate resultTree
			if resultTree == nil {
				t.Fatal("CreateNucleotideARGTree returned nil")
			}

			// Read expected output
			expectedTree, err := readARGOutputFile(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Compare results
			if !compareARGTrees(resultTree, expectedTree) {
				t.Errorf("ARG tree result does not match expected output for file %s", inputFile)
				t.Errorf("Input Matrix:\n%v", formatMatrix(matrix))
				t.Errorf("Result Tree:\n%v", formatARGTree(resultTree))
				t.Errorf("Expected Tree:\n%v", formatARGTree(expectedTree))
			}
		})
	}
}

// Charlotte Barretto
func readARGInputFile(filename string) (DistanceMatrix, []string, string, string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, "", "", err
	}
	defer file.Close()

	var matrix DistanceMatrix
	var genomeStrains []string
	var searchString, argName string

	scanner := bufio.NewScanner(file)
	section := ""

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "//") {
			continue
		}

		if strings.HasPrefix(line, "#") {
			section = strings.TrimSpace(strings.TrimPrefix(line, "#"))
			continue
		}

		switch section {
		case "Distance Matrix":
			row := make([]float64, 0)
			for _, numStr := range strings.Fields(line) {
				num, err := strconv.ParseFloat(numStr, 64)
				if err != nil {
					return nil, nil, "", "", fmt.Errorf("invalid number in matrix: %s", numStr)
				}
				row = append(row, num)
			}
			matrix = append(matrix, row)
		case "Genome Strains":
			genomeStrains = append(genomeStrains, line)
		case "Search String":
			searchString = line
		case "ARG Name":
			argName = line
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, nil, "", "", err
	}

	// Validate input
	if len(matrix) == 0 || len(genomeStrains) == 0 || searchString == "" || argName == "" {
		return nil, nil, "", "", fmt.Errorf("missing required input data")
	}

	// Validate matrix dimensions
	for _, row := range matrix {
		if len(row) != len(matrix) {
			return nil, nil, "", "", fmt.Errorf("matrix is not square")
		}
	}

	// Validate number of genomes matches matrix size
	if len(genomeStrains) != len(matrix) {
		return nil, nil, "", "", fmt.Errorf("number of genomes does not match matrix size")
	}

	return matrix, genomeStrains, searchString, argName, nil
}

// Charlotte Barretto
func readARGOutputFile(filename string) (ARGTree, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var tree ARGTree
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := scanner.Text()
		if line == "" {
			continue
		}

		// Parse node information
		// Format: Node X: strainName sequence [indices] -> [descendants]
		node := &ARGNode{
			ARGpointer: &ARG{
				IndexInGenome: [2]int{-1, -1}, // Initialize with default values
			},
		}

		// Split the line into parts
		parts := strings.Split(line, ":")
		if len(parts) != 2 {
			continue
		}

		// Parse strain information
		strainInfo := strings.TrimSpace(parts[1])
		if strings.Contains(strainInfo, "->") {
			// This is an ancestor node
			ancestorParts := strings.Split(strainInfo, "->")
			node.GenomeStrain = strings.TrimSpace(ancestorParts[0])
			node.ARGpointer = nil // Ancestor nodes don't have ARG pointers

			// Parse descendants
			descendantStr := strings.Trim(ancestorParts[1], "[] ")
			descendants := strings.Split(descendantStr, ",")
			if len(descendants) == 2 {
				desc1, desc2 := strings.TrimSpace(descendants[0]), strings.TrimSpace(descendants[1])
				node.Descendents1 = &ARGNode{GenomeStrain: desc1}
				node.Descendents2 = &ARGNode{GenomeStrain: desc2}
			}
		} else {
			// This is a leaf node
			parts = strings.Fields(strainInfo)
			node.GenomeStrain = parts[0]
			if len(parts) > 1 {
				node.ARGpointer.GenomeStrain = parts[0]
				node.ARGpointer.Nucleotides = parts[1]
				if len(parts) > 2 {
					// Parse indices
					indicesStr := strings.Trim(parts[2], "[]")
					if indicesStr != "" {
						indices := strings.Split(indicesStr, ",")
						// Only store up to 2 indices since IndexInGenome is [2]int
						for i := 0; i < len(indices) && i < 2; i++ {
							num, err := strconv.Atoi(indices[i])
							if err != nil {
								continue
							}
							node.ARGpointer.IndexInGenome[i] = num
						}
					}
				}
			}
		}
		tree = append(tree, node)
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return tree, nil
}

// Charlotte Barretto
func compareARGTrees(tree1, tree2 ARGTree) bool {
	if len(tree1) != len(tree2) {
		return false
	}

	for i := range tree1 {
		if !compareARGNodes(tree1[i], tree2[i]) {
			return false
		}
	}
	return true
}

// Charlotte Barretto
func compareARGNodes(node1, node2 *ARGNode) bool {
	if node1.GenomeStrain != node2.GenomeStrain {
		return false
	}

	// Compare ARG pointers
	if (node1.ARGpointer == nil) != (node2.ARGpointer == nil) {
		return false
	}
	if node1.ARGpointer != nil {
		if !compareARGs(node1.ARGpointer, node2.ARGpointer) {
			return false
		}
	}

	// Compare descendants
	if !compareDescendants(node1.Descendents1, node2.Descendents1) {
		return false
	}
	if !compareDescendants(node1.Descendents2, node2.Descendents2) {
		return false
	}

	return true
}

// Charlotte Barretto
func compareDescendants(desc1, desc2 *ARGNode) bool {
	if (desc1 == nil) != (desc2 == nil) {
		return false
	}
	if desc1 != nil {
		return desc1.GenomeStrain == desc2.GenomeStrain
	}
	return true
}

// Charlotte Barretto
func compareARGs(arg1, arg2 *ARG) bool {
	if arg1.GenomeStrain != arg2.GenomeStrain {
		return false
	}
	if !reflect.DeepEqual(arg1.IndexInGenome, arg2.IndexInGenome) {
		return false
	}
	if arg1.Nucleotides != arg2.Nucleotides {
		return false
	}
	return true
}

// Charlotte Barretto
func formatMatrix(matrix DistanceMatrix) string {
	var sb strings.Builder
	for _, row := range matrix {
		for _, val := range row {
			sb.WriteString(fmt.Sprintf("%.6f ", val))
		}
		sb.WriteString("\n")
	}
	return sb.String()
}

// Charlotte Barretto
func formatARGTree(tree ARGTree) string {
	var sb strings.Builder
	for i, node := range tree {
		sb.WriteString(fmt.Sprintf("Node %d: %s", i, node.GenomeStrain))
		if node.ARGpointer != nil {
			sb.WriteString(fmt.Sprintf(" %s %v", node.ARGpointer.Nucleotides, node.ARGpointer.IndexInGenome))
		}
		if node.Descendents1 != nil && node.Descendents2 != nil {
			sb.WriteString(fmt.Sprintf(" -> [%s,%s]", node.Descendents1.GenomeStrain, node.Descendents2.GenomeStrain))
		}
		sb.WriteString("\n")
	}
	return sb.String()
}

func TestCreateProteinTree(t *testing.T) {
	inputDir := "Tests/CreateProteinTree/Input"
	outputDir := "Tests/CreateProteinTree/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			matrix, genomeStrains, _, _, err := readARGInputFile(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Read AA sequences
			aaSequences, err := readAASequences(inputFile)
			if err != nil {
				t.Fatalf("Failed to read AA sequences: %v", err)
			}

			// Read ARG tree from input file
			argTree, err := readARGTreeFromInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read ARG tree: %v", err)
			}

			percentIdentity := 0.95

			// Create Protein Tree
			resultTree := CreateProteinTree(matrix, aaSequences, genomeStrains, genomeStrains, argTree, percentIdentity)

			// Validate resultTree
			if resultTree == nil {
				t.Fatal("CreateProteinTree returned nil")
			}

			// Read expected output
			expectedTree, err := readProteinTreeOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Compare results
			if !compareProteinTrees(resultTree, expectedTree) {
				t.Errorf("Protein tree result does not match expected output for file %s", inputFile)
				t.Errorf("Input Matrix:\n%v", formatMatrix(matrix))
				t.Errorf("Result Tree:\n%v", formatProteinTree(resultTree))
				t.Errorf("Expected Tree:\n%v", formatProteinTree(expectedTree))
			}
		})
	}
}

func readAASequences(filename string) ([]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var sequences []string
	scanner := bufio.NewScanner(file)
	section := ""

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "//") {
			continue
		}

		if strings.HasPrefix(line, "#") {
			section = strings.TrimSpace(strings.TrimPrefix(line, "#"))
			continue
		}

		if section == "AA Sequences" {
			sequences = append(sequences, line)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	if len(sequences) == 0 {
		return nil, fmt.Errorf("no AA sequences found in file")
	}

	return sequences, nil
}

func readARGTreeFromInput(filename string) (ARGTree, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var tree ARGTree
	scanner := bufio.NewScanner(file)
	section := ""
	inARGTree := false

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "//") {
			continue
		}

		if strings.HasPrefix(line, "#") {
			section = strings.TrimSpace(strings.TrimPrefix(line, "#"))
			inARGTree = section == "ARG Tree"
			continue
		}

		if inARGTree {
			node := &ARGNode{
				ARGpointer: &ARG{
					IndexInGenome: [2]int{-1, -1},
				},
			}

			parts := strings.Split(line, ":")
			if len(parts) != 2 {
				continue
			}

			strainInfo := strings.TrimSpace(parts[1])
			if strings.Contains(strainInfo, "->") {
				ancestorParts := strings.Split(strainInfo, "->")
				node.GenomeStrain = strings.TrimSpace(ancestorParts[0])
				node.ARGpointer = nil

				descendantStr := strings.Trim(ancestorParts[1], "[] ")
				descendants := strings.Split(descendantStr, ",")
				if len(descendants) == 2 {
					desc1, desc2 := strings.TrimSpace(descendants[0]), strings.TrimSpace(descendants[1])
					node.Descendents1 = &ARGNode{GenomeStrain: desc1}
					node.Descendents2 = &ARGNode{GenomeStrain: desc2}
				}
			} else {
				parts = strings.Fields(strainInfo)
				node.GenomeStrain = parts[0]
				if len(parts) > 1 {
					node.ARGpointer.GenomeStrain = parts[0]
					node.ARGpointer.Nucleotides = parts[1]
					if len(parts) > 2 {
						indicesStr := strings.Trim(parts[2], "[]")
						indices := strings.Split(indicesStr, ",")
						for i := 0; i < len(indices) && i < 2; i++ {
							num, err := strconv.Atoi(indices[i])
							if err != nil {
								continue
							}
							node.ARGpointer.IndexInGenome[i] = num
						}
					}
				}
			}
			tree = append(tree, node)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	if len(tree) == 0 {
		return nil, fmt.Errorf("no ARG tree found in file")
	}

	return tree, nil
}

func readProteinTreeOutput(filename string) (ProteinTree, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var tree ProteinTree
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		// Split on first colon only
		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			return nil, fmt.Errorf("invalid protein tree format: %s", line)
		}

		//nodeStr := strings.TrimSpace(parts[0])     // "Node X"
		nodeDetails := strings.TrimSpace(parts[1]) // Rest of the line

		node := &ProteinNode{
			Protein: &Protein{},
		}

		if strings.Contains(nodeDetails, "Ancestor Protein") {
			// Handle ancestor node
			node.Label = nodeDetails
			if strings.Contains(nodeDetails, "->") {
				// Don't need to process descendants as they're handled in comparison
				node.Descendants1 = &ProteinNode{}
				node.Descendants2 = &ProteinNode{}
			}
		} else {
			// Handle leaf node
			parts = strings.Fields(nodeDetails)
			if len(parts) > 0 {
				node.Protein.AminoAcids = parts[0]
				node.Protein.ARG = &ARG{
					Nucleotides:   parts[0],
					IndexInGenome: [2]int{-1, -1},
				}

				if len(parts) > 1 {
					indicesStr := strings.Trim(parts[1], "[]")
					indices := strings.Split(indicesStr, ",")
					if len(indices) == 2 {
						start, _ := strconv.Atoi(indices[0])
						end, _ := strconv.Atoi(indices[1])
						node.Protein.ARG.IndexInGenome = [2]int{start, end}
					}
				}
			}
		}

		tree = append(tree, node)
	}

	return tree, nil
}

func compareProteinTrees(tree1, tree2 ProteinTree) bool {
	if len(tree1) != len(tree2) {
		return false
	}

	for i := range tree1 {
		if !compareProteinNodes(tree1[i], tree2[i]) {
			return false
		}
	}
	return true
}

func compareProteinNodes(node1, node2 *ProteinNode) bool {
	if (node1 == nil) != (node2 == nil) {
		return false
	}
	if node1 == nil {
		return true
	}

	// Compare labels for ancestor nodes
	if strings.Contains(node1.Label, "Ancestor") || strings.Contains(node2.Label, "Ancestor") {
		return node1.Label == node2.Label
	}

	// Compare proteins for leaf nodes
	return compareProteins(node1.Protein, node2.Protein)
}

func compareProteins(p1, p2 *Protein) bool {
	if p1 == nil || p2 == nil {
		return p1 == p2
	}

	return p1.AminoAcids == p2.AminoAcids &&
		compareARGsProteinTree(p1.ARG, p2.ARG)
}

func compareARGsProteinTree(arg1, arg2 *ARG) bool {
	if arg1 == nil || arg2 == nil {
		return arg1 == arg2
	}

	return arg1.Nucleotides == arg2.Nucleotides &&
		arg1.IndexInGenome == arg2.IndexInGenome
}

func formatProteinTree(tree ProteinTree) string {
	var sb strings.Builder
	for i, node := range tree {
		if strings.Contains(node.Label, "Ancestor") {
			fmt.Fprintf(&sb, "Node %d: %s\n", i, node.Label)
		} else {
			indices := ""
			if node.Protein != nil && node.Protein.ARG != nil {
				indices = fmt.Sprintf(" [%d,%d]",
					node.Protein.ARG.IndexInGenome[0],
					node.Protein.ARG.IndexInGenome[1])
			}
			fmt.Fprintf(&sb, "Node %d: %s%s\n",
				i, node.Protein.AminoAcids, indices)
		}
	}
	return sb.String()
}

func TestSaveArgRegions(t *testing.T) {
	inputDir := "Tests/SaveArgRegions/Input"
	outputDir := "Tests/SaveArgRegions/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			genome, args, err := readSaveArgRegionsInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Explicitly print debug information
			t.Logf("Input Genome: %s", genome)
			t.Logf("Input Args (Before): %v", args)

			// Call SaveArgRegions function
			resultRegions := SaveArgRegions(genome, args)

			// Read expected output
			expectedRegions, err := readSaveArgRegionsOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Explicitly print more debug information
			t.Logf("Input Args (After): %v", args)
			t.Logf("Result Regions: %v", resultRegions)
			t.Logf("Expected Regions: %v", expectedRegions)

			// Compare results
			if !reflect.DeepEqual(resultRegions, expectedRegions) {
				t.Errorf("SaveArgRegions result does not match expected output for file %s", inputFile)
			}
		})
	}
}

func readSaveArgRegionsInput(filename string) (string, [][2]int, error) {
	file, err := os.Open(filename)
	if err != nil {
		return "", nil, err
	}
	defer file.Close()

	var genome string
	var args [][2]int

	scanner := bufio.NewScanner(file)
	section := ""

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "//") {
			continue
		}

		if strings.HasPrefix(line, "#") {
			section = strings.TrimSpace(strings.TrimPrefix(line, "#"))
			continue
		}

		switch section {
		case "Genome":
			genome = line
		case "Args":
			// Ensure only one set of args is processed
			args = [][2]int{}

			coord := strings.Split(line, ",")
			if len(coord) != 2 {
				return "", nil, fmt.Errorf("invalid args format: %s", line)
			}

			start, err := strconv.Atoi(strings.TrimSpace(coord[0]))
			if err != nil {
				return "", nil, fmt.Errorf("invalid start coordinate: %s", coord[0])
			}

			end, err := strconv.Atoi(strings.TrimSpace(coord[1]))
			if err != nil {
				return "", nil, fmt.Errorf("invalid end coordinate: %s", coord[1])
			}

			args = append(args, [2]int{start, end})
		}
	}

	if err := scanner.Err(); err != nil {
		return "", nil, err
	}

	// Ensure genome and args are not empty
	if genome == "" || len(args) == 0 {
		return "", nil, fmt.Errorf("missing genome or args")
	}

	return genome, args, nil
}

func readSaveArgRegionsOutput(filename string) ([]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var regions []string
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line != "" {
			regions = append(regions, line)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return regions, nil
}

func TestFindARGInTree(t *testing.T) {
	// Create test cases with ARGTree structures
	tests := []struct {
		name         string
		searchStrain string
		argTree      ARGTree
		expected     *ARG
	}{
		{
			name:         "Find existing strain",
			searchStrain: "Strain1",
			argTree: ARGTree{
				&ARGNode{
					ARGpointer: &ARG{
						Nucleotides:   "ATCG",
						GenomeStrain:  "Strain1",
						IndexInGenome: [2]int{0, 24},
					},
				},
				&ARGNode{
					ARGpointer: &ARG{
						Nucleotides:   "GCTA",
						GenomeStrain:  "Strain2",
						IndexInGenome: [2]int{0, -1},
					},
				},
			},
			expected: &ARG{
				Nucleotides:   "ATCG",
				GenomeStrain:  "Strain1",
				IndexInGenome: [2]int{0, 24},
			},
		},
		{
			name:         "Strain not found",
			searchStrain: "NonexistentStrain",
			argTree: ARGTree{
				&ARGNode{
					ARGpointer: &ARG{
						Nucleotides:   "ATCG",
						GenomeStrain:  "Strain1",
						IndexInGenome: [2]int{0, 24},
					},
				},
			},
			expected: &ARG{
				ArgLabel: "Ancestor Strain",
			},
		},
		{
			name:         "Empty tree",
			searchStrain: "Strain1",
			argTree:      ARGTree{},
			expected: &ARG{
				ArgLabel: "Ancestor Strain",
			},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := FindARGInTree(tt.searchStrain, tt.argTree)

			// Compare result with expected
			if !compareARGsFindARG(result, tt.expected) {
				t.Errorf("FindARGInTree() = %v, want %v", formatARG(result), formatARG(tt.expected))
			}
		})
	}
}

func compareARGsFindARG(arg1, arg2 *ARG) bool {
	if arg1 == nil || arg2 == nil {
		return arg1 == arg2
	}

	// For ancestor strain case
	if arg1.ArgLabel == "Ancestor Strain" || arg2.ArgLabel == "Ancestor Strain" {
		return arg1.ArgLabel == arg2.ArgLabel
	}

	return arg1.GenomeStrain == arg2.GenomeStrain &&
		arg1.Nucleotides == arg2.Nucleotides &&
		arg1.IndexInGenome == arg2.IndexInGenome
}

func formatARG(arg *ARG) string {
	if arg == nil {
		return "nil"
	}

	if arg.ArgLabel == "Ancestor Strain" {
		return fmt.Sprintf("ARG{ArgLabel: %s}", arg.ArgLabel)
	}

	return fmt.Sprintf("ARG{GenomeStrain: %s, Nucleotides: %s, IndexInGenome: [%d,%d]}",
		arg.GenomeStrain, arg.Nucleotides, arg.IndexInGenome[0], arg.IndexInGenome[1])
}

func TestInitializeProteinTree(t *testing.T) {
	inputDir := "Tests/InitializeProteinTree/Input"
	outputDir := "Tests/InitializeProteinTree/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			aaSequences, genomeStrains, genomeNames, argTree, err := readInitializeProteinTreeInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Initialize Protein Tree
			resultTree := InitializeProteinTree(aaSequences, genomeStrains, genomeNames, argTree)

			// Read expected output
			expectedTree, err := readInitializeProteinTreeOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Compare results
			if compareProteinTreesInitialize(resultTree, expectedTree) {
				t.Errorf("InitializeProteinTree result does not match expected output for file %s", inputFile)
				t.Errorf("Result Tree:\n%s", formatProteinTreeInitialize(resultTree))
				t.Errorf("Expected Tree:\n%s", formatProteinTreeInitialize(expectedTree))
			}

		})
	}
}

func readInitializeProteinTreeInput(filename string) ([]string, []string, []string, ARGTree, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, nil, nil, err
	}
	defer file.Close()

	var aaSequences, genomeStrains, genomeNames []string
	var argTree ARGTree
	scanner := bufio.NewScanner(file)
	section := ""

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "//") {
			continue
		}

		if strings.HasPrefix(line, "#") {
			section = strings.TrimSpace(strings.TrimPrefix(line, "#"))
			continue
		}

		switch section {
		case "AA Sequences":
			aaSequences = append(aaSequences, line)
		case "Genome Strains":
			genomeStrains = append(genomeStrains, line)
		case "Genome Names":
			genomeNames = append(genomeNames, line)
		case "ARG Tree":
			node := &ARGNode{
				ARGpointer: &ARG{
					IndexInGenome: [2]int{-1, -1},
				},
			}

			parts := strings.Split(line, ":")
			if len(parts) != 2 {
				continue
			}

			strainInfo := strings.TrimSpace(parts[1])
			if strings.Contains(strainInfo, "->") {
				// Ancestor node
				ancestorParts := strings.Split(strainInfo, "->")
				node.GenomeStrain = strings.TrimSpace(ancestorParts[0])
				node.ARGpointer = nil

				descendantStr := strings.Trim(ancestorParts[1], "[] ")
				descendants := strings.Split(descendantStr, ",")
				if len(descendants) == 2 {
					desc1, desc2 := strings.TrimSpace(descendants[0]), strings.TrimSpace(descendants[1])
					node.Descendents1 = &ARGNode{GenomeStrain: desc1}
					node.Descendents2 = &ARGNode{GenomeStrain: desc2}
				}
			} else {
				// Leaf node
				parts = strings.Fields(strainInfo)
				if len(parts) >= 2 {
					node.GenomeStrain = parts[0]
					node.ARGpointer.GenomeStrain = parts[0]
					node.ARGpointer.Nucleotides = parts[1]
					if len(parts) > 2 {
						indicesStr := strings.Trim(parts[2], "[]")
						indices := strings.Split(indicesStr, ",")
						for i := 0; i < len(indices) && i < 2; i++ {
							num, err := strconv.Atoi(indices[i])
							if err != nil {
								continue
							}
							node.ARGpointer.IndexInGenome[i] = num
						}
					}
				}
			}
			argTree = append(argTree, node)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, nil, nil, nil, err
	}

	return aaSequences, genomeStrains, genomeNames, argTree, nil
}

func readInitializeProteinTreeOutput(filename string) (ProteinTree, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var tree ProteinTree
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			return nil, fmt.Errorf("invalid protein tree format: %s", line)
		}

		nodeDetails := strings.TrimSpace(parts[1])
		node := &ProteinNode{
			Protein: &Protein{},
		}

		if strings.Contains(nodeDetails, "Ancestor Protein") {
			node.Label = nodeDetails
		} else {
			parts = strings.Fields(nodeDetails)
			if len(parts) >= 1 {
				node.Protein.AminoAcids = parts[0]
				if len(parts) > 1 {
					indicesStr := strings.Trim(parts[1], "[]")
					indices := strings.Split(indicesStr, ",")
					if len(indices) == 2 {
						start, _ := strconv.Atoi(indices[0])
						end, _ := strconv.Atoi(indices[1])
						node.Protein.ARG = &ARG{
							IndexInGenome: [2]int{start, end},
						}
					}
				}
			}
		}

		tree = append(tree, node)
	}

	return tree, nil
}

func compareProteinTreesInitialize(tree1, tree2 ProteinTree) bool {
	if len(tree1) != len(tree2) {
		return false
	}

	for i := range tree1 {
		if !compareProteinNodesInitialize(tree1[i], tree2[i]) {
			return false
		}
	}
	return true
}

func compareProteinNodesInitialize(node1, node2 *ProteinNode) bool {
	if node1 == nil || node2 == nil {
		return node1 == node2
	}

	if node1.Label != node2.Label {
		return false
	}

	return compareProteinsInitialize(node1.Protein, node2.Protein)
}

func compareProteinsInitialize(p1, p2 *Protein) bool {
	if p1 == nil || p2 == nil {
		return p1 == p2
	}

	if p1.AminoAcids != p2.AminoAcids {
		return false
	}

	return compareARGsInitialize(p1.ARG, p2.ARG)
}

func compareARGsInitialize(arg1, arg2 *ARG) bool {
	if arg1 == nil || arg2 == nil {
		return arg1 == arg2
	}

	return arg1.IndexInGenome == arg2.IndexInGenome
}

func formatProteinTreeInitialize(tree ProteinTree) string {
	var sb strings.Builder
	for i, node := range tree {
		if node.Label != "" {
			fmt.Fprintf(&sb, "Node %d: %s\n", i, node.Label)
		} else {
			indices := ""
			if node.Protein != nil && node.Protein.ARG != nil {
				indices = fmt.Sprintf(" [%d,%d]",
					node.Protein.ARG.IndexInGenome[0],
					node.Protein.ARG.IndexInGenome[1])
			}
			fmt.Fprintf(&sb, "Node %d: %s%s\n",
				i, node.Protein.AminoAcids, indices)
		}
	}
	return sb.String()
}

func TestCountLeavesProtein(t *testing.T) {
	inputDir := "Tests/CountLeavesProtein/Input"
	outputDir := "Tests/CountLeavesProtein/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			proteinTree, err := readCountLeavesInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Count leaves
			result := CountLeavesProtein(proteinTree)

			// Read expected output
			expected, err := readCountLeavesOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Compare results
			if result != expected {
				t.Errorf("CountLeavesProtein = %v, want %v", result, expected)
			}
		})
	}
}

func readCountLeavesInput(filename string) (*ProteinNode, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var nodes []*ProteinNode
	nodeMap := make(map[int]*ProteinNode)
	scanner := bufio.NewScanner(file)

	// First pass: Create all nodes
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}

		nodeStr := strings.TrimSpace(parts[0])
		nodeNum, _ := strconv.Atoi(strings.TrimPrefix(nodeStr, "Node "))
		nodeDetails := strings.TrimSpace(parts[1])

		var node *ProteinNode
		if strings.Contains(nodeDetails, "Ancestor") {
			node = &ProteinNode{
				Label: nodeDetails,
			}
		} else {
			node = &ProteinNode{
				Protein: &Protein{
					AminoAcids: nodeDetails,
				},
			}
		}
		nodeMap[nodeNum] = node
		nodes = append(nodes, node)
	}

	// Reset file for second pass
	file.Seek(0, 0)
	scanner = bufio.NewScanner(file)

	// Second pass: Connect descendants
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if !strings.Contains(line, "->") {
			continue
		}

		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}

		nodeStr := strings.TrimSpace(parts[0])
		nodeNum, _ := strconv.Atoi(strings.TrimPrefix(nodeStr, "Node "))

		// Get descendants
		nodeDetails := strings.TrimSpace(parts[1])
		descendantsPart := strings.Split(nodeDetails, "->")[1]
		descendantsPart = strings.Trim(descendantsPart, "[] ")
		descendants := strings.Split(descendantsPart, ",")

		if len(descendants) == 2 {
			desc1, _ := strconv.Atoi(strings.TrimSpace(descendants[0]))
			desc2, _ := strconv.Atoi(strings.TrimSpace(descendants[1]))

			if node, exists := nodeMap[nodeNum]; exists {
				node.Descendants1 = nodeMap[desc1]
				node.Descendants2 = nodeMap[desc2]
			}
		}
	}

	// Return the root (last node)
	if len(nodes) > 0 {
		return nodes[len(nodes)-1], nil
	}
	return nil, fmt.Errorf("no nodes found in input file")
}

func readCountLeavesOutput(filename string) (int, error) {
	file, err := os.Open(filename)
	if err != nil {
		return 0, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	if scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		count, err := strconv.Atoi(line)
		if err != nil {
			return 0, fmt.Errorf("invalid output format: %v", err)
		}
		return count, nil
	}

	return 0, fmt.Errorf("empty output file")
}

func TestInitializeClustersForProtein(t *testing.T) {
	inputDir := "Tests/InitializeClustersForProtein/Input"
	outputDir := "Tests/InitializeClustersForProtein/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			proteinTree, err := readInitializeClustersInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Initialize clusters
			resultClusters := proteinTree.InitializeClustersForProtein()

			// Read expected output
			expectedClusters, err := readInitializeClustersOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Compare results
			if !compareProteinClusters(resultClusters, expectedClusters) {
				t.Errorf("InitializeClustersForProtein result does not match expected output for file %s", inputFile)
				t.Errorf("Result Clusters:\n%s", formatProteinClusters(resultClusters))
				t.Errorf("Expected Clusters:\n%s", formatProteinClusters(expectedClusters))
			}
		})
	}
}

func readInitializeClustersInput(filename string) (ProteinTree, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var tree ProteinTree
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}

		nodeDetails := strings.TrimSpace(parts[1])
		node := &ProteinNode{}

		if strings.Contains(nodeDetails, "Ancestor") {
			node.Label = nodeDetails
		} else {
			parts = strings.Fields(nodeDetails)
			if len(parts) >= 1 {
				node.Protein = &Protein{
					AminoAcids: parts[0],
				}
				if len(parts) > 1 {
					indicesStr := strings.Trim(parts[1], "[]")
					indices := strings.Split(indicesStr, ",")
					if len(indices) == 2 {
						start, _ := strconv.Atoi(indices[0])
						end, _ := strconv.Atoi(indices[1])
						node.Protein.ARG = &ARG{
							IndexInGenome: [2]int{start, end},
						}
					}
				}
			}
		}
		tree = append(tree, node)
	}

	return tree, nil
}

func readInitializeClustersOutput(filename string) ([]*ProteinNode, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var clusters []*ProteinNode
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}

		nodeDetails := strings.TrimSpace(parts[1])
		node := &ProteinNode{}

		parts = strings.Fields(nodeDetails)
		if len(parts) >= 1 {
			node.Protein = &Protein{
				AminoAcids: parts[0],
			}
			if len(parts) > 1 {
				indicesStr := strings.Trim(parts[1], "[]")
				indices := strings.Split(indicesStr, ",")
				if len(indices) == 2 {
					start, _ := strconv.Atoi(indices[0])
					end, _ := strconv.Atoi(indices[1])
					node.Protein.ARG = &ARG{
						IndexInGenome: [2]int{start, end},
					}
				}
			}
		}
		clusters = append(clusters, node)
	}

	return clusters, nil
}

func compareProteinClusters(clusters1, clusters2 []*ProteinNode) bool {
	if len(clusters1) != len(clusters2) {
		return false
	}

	for i := range clusters1 {
		if !compareProteinNodes(clusters1[i], clusters2[i]) {
			return false
		}
	}
	return true
}

func formatProteinClusters(clusters []*ProteinNode) string {
	var sb strings.Builder
	for i, node := range clusters {
		if node.Protein != nil {
			indices := ""
			if node.Protein.ARG != nil {
				indices = fmt.Sprintf(" [%d,%d]",
					node.Protein.ARG.IndexInGenome[0],
					node.Protein.ARG.IndexInGenome[1])
			}
			fmt.Fprintf(&sb, "Cluster %d: %s%s\n",
				i, node.Protein.AminoAcids, indices)
		}
	}
	return sb.String()
}

func TestCountLeavesARG(t *testing.T) {
	inputDir := "Tests/CountLeavesARG/Input"
	outputDir := "Tests/CountLeavesARG/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			argNode, err := readCountLeavesARGInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Count leaves
			result := CountLeavesARG(argNode)

			// Read expected output
			expected, err := readCountLeavesARGOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Compare results
			if result != expected {
				t.Errorf("CountLeavesARG = %v, want %v", result, expected)
			}
		})
	}
}

func readCountLeavesARGInput(filename string) (*ARGNode, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var nodes []*ARGNode
	nodeMap := make(map[int]*ARGNode)
	scanner := bufio.NewScanner(file)

	// First pass: Create all nodes
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}

		nodeStr := strings.TrimSpace(parts[0])
		nodeNum, _ := strconv.Atoi(strings.TrimPrefix(nodeStr, "Node "))
		nodeDetails := strings.TrimSpace(parts[1])

		var node *ARGNode
		if strings.Contains(nodeDetails, "->") {
			// Ancestor node
			node = &ARGNode{}
		} else {
			// Leaf node
			parts = strings.Fields(nodeDetails)
			node = &ARGNode{
				ARGpointer: &ARG{
					GenomeStrain: parts[0],
					Nucleotides:  parts[1],
				},
			}
			if len(parts) > 2 {
				indicesStr := strings.Trim(parts[2], "[]")
				indices := strings.Split(indicesStr, ",")
				if len(indices) == 2 {
					start, _ := strconv.Atoi(indices[0])
					end, _ := strconv.Atoi(indices[1])
					node.ARGpointer.IndexInGenome = [2]int{start, end}
				}
			}
		}
		nodeMap[nodeNum] = node
		nodes = append(nodes, node)
	}

	// Reset file for second pass
	file.Seek(0, 0)
	scanner = bufio.NewScanner(file)

	// Second pass: Connect descendants
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if !strings.Contains(line, "->") {
			continue
		}

		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}

		nodeStr := strings.TrimSpace(parts[0])
		nodeNum, _ := strconv.Atoi(strings.TrimPrefix(nodeStr, "Node "))

		// Get descendants
		nodeDetails := strings.TrimSpace(parts[1])
		descendantsPart := strings.Split(nodeDetails, "->")[1]
		descendantsPart = strings.Trim(descendantsPart, "[] ")
		descendants := strings.Split(descendantsPart, ",")

		if len(descendants) == 2 {
			desc1, _ := strconv.Atoi(strings.TrimSpace(descendants[0]))
			desc2, _ := strconv.Atoi(strings.TrimSpace(descendants[1]))

			if node, exists := nodeMap[nodeNum]; exists {
				node.Descendents1 = nodeMap[desc1]
				node.Descendents2 = nodeMap[desc2]
			}
		}
	}

	// Return the root (last node)
	if len(nodes) > 0 {
		return nodes[len(nodes)-1], nil
	}
	return nil, fmt.Errorf("no nodes found in input file")
}

func readCountLeavesARGOutput(filename string) (int, error) {
	file, err := os.Open(filename)
	if err != nil {
		return 0, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	if scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		count, err := strconv.Atoi(line)
		if err != nil {
			return 0, fmt.Errorf("invalid output format: %v", err)
		}
		return count, nil
	}

	return 0, fmt.Errorf("empty output file")
}

func TestInitializeClustersForARG(t *testing.T) {
	inputDir := "Tests/InitializeClustersForARG/Input"
	outputDir := "Tests/InitializeClustersForARG/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			argTree, err := readInitializeClustersARGInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Initialize clusters
			resultClusters := argTree.InitializeClustersForARG()

			// Read expected output
			expectedClusters, err := readInitializeClustersARGOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Compare results
			if !compareARGClusters(resultClusters, expectedClusters) {
				t.Errorf("InitializeClustersForARG result does not match expected output for file %s", inputFile)
				t.Errorf("Result Clusters:\n%s", formatARGClusters(resultClusters))
				t.Errorf("Expected Clusters:\n%s", formatARGClusters(expectedClusters))
			}
		})
	}
}

func readInitializeClustersARGInput(filename string) (ARGTree, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var tree ARGTree
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}

		nodeDetails := strings.TrimSpace(parts[1])
		node := &ARGNode{}

		if strings.Contains(nodeDetails, "->") {
			// Ancestor node
			node.ARGpointer = nil
		} else {
			// Leaf node
			parts = strings.Fields(nodeDetails)
			if len(parts) >= 2 {
				node.ARGpointer = &ARG{
					GenomeStrain: parts[0],
					Nucleotides:  parts[1],
				}
				if len(parts) > 2 {
					indicesStr := strings.Trim(parts[2], "[]")
					indices := strings.Split(indicesStr, ",")
					if len(indices) == 2 {
						start, _ := strconv.Atoi(indices[0])
						end, _ := strconv.Atoi(indices[1])
						node.ARGpointer.IndexInGenome = [2]int{start, end}
					}
				}
			}
		}
		tree = append(tree, node)
	}

	return tree, nil
}

func readInitializeClustersARGOutput(filename string) ([]*ARGNode, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var clusters []*ARGNode
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}

		nodeDetails := strings.TrimSpace(parts[1])
		parts = strings.Fields(nodeDetails)

		node := &ARGNode{
			ARGpointer: &ARG{
				GenomeStrain: parts[0],
				Nucleotides:  parts[1],
			},
		}

		if len(parts) > 2 {
			indicesStr := strings.Trim(parts[2], "[]")
			indices := strings.Split(indicesStr, ",")
			if len(indices) == 2 {
				start, _ := strconv.Atoi(indices[0])
				end, _ := strconv.Atoi(indices[1])
				node.ARGpointer.IndexInGenome = [2]int{start, end}
			}
		}
		clusters = append(clusters, node)
	}

	return clusters, nil
}

func compareARGClusters(clusters1, clusters2 []*ARGNode) bool {
	if len(clusters1) != len(clusters2) {
		return false
	}

	for i := range clusters1 {
		if !compareARGNodesCluster(clusters1[i], clusters2[i]) {
			return false
		}
	}
	return true
}

func compareARGNodesCluster(node1, node2 *ARGNode) bool {
	if node1 == nil || node2 == nil {
		return node1 == node2
	}

	if node1.ARGpointer == nil || node2.ARGpointer == nil {
		return node1.ARGpointer == node2.ARGpointer
	}

	return node1.ARGpointer.GenomeStrain == node2.ARGpointer.GenomeStrain &&
		node1.ARGpointer.Nucleotides == node2.ARGpointer.Nucleotides &&
		node1.ARGpointer.IndexInGenome == node2.ARGpointer.IndexInGenome
}

func formatARGClusters(clusters []*ARGNode) string {
	var sb strings.Builder
	for i, node := range clusters {
		if node.ARGpointer != nil {
			indices := fmt.Sprintf("[%d,%d]",
				node.ARGpointer.IndexInGenome[0],
				node.ARGpointer.IndexInGenome[1])
			fmt.Fprintf(&sb, "Cluster %d: %s %s %s\n",
				i, node.ARGpointer.GenomeStrain,
				node.ARGpointer.Nucleotides, indices)
		}
	}
	return sb.String()
}

func TestMinIdentityDifference(t *testing.T) {
	inputDir := "Tests/MinIdentityDifference/Input"
	outputDir := "Tests/MinIdentityDifference/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			matrix, err := readMinIdentityDifferenceInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Get minimum identity difference
			resultRow, resultCol, resultDiff := MinIdentityDifference(matrix)

			// Read expected output
			expectedRow, expectedCol, expectedDiff, err := readMinIdentityDifferenceOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Compare results
			if resultRow != expectedRow || resultCol != expectedCol || math.Abs(resultDiff-expectedDiff) > 1e-6 {
				t.Errorf("MinIdentityDifference = (%d, %d, %f), want (%d, %d, %f)",
					resultRow, resultCol, resultDiff,
					expectedRow, expectedCol, expectedDiff)
			}
		})
	}
}

func readMinIdentityDifferenceInput(filename string) (DistanceMatrix, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var matrix DistanceMatrix
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		var row []float64
		values := strings.Fields(line)
		for _, val := range values {
			num, err := strconv.ParseFloat(val, 64)
			if err != nil {
				return nil, fmt.Errorf("invalid number format: %v", err)
			}
			row = append(row, num)
		}
		matrix = append(matrix, row)
	}

	return matrix, nil
}

func readMinIdentityDifferenceOutput(filename string) (int, int, float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return 0, 0, 0, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	if !scanner.Scan() {
		return 0, 0, 0, fmt.Errorf("empty output file")
	}

	line := strings.TrimSpace(scanner.Text())
	parts := strings.Fields(line)
	if len(parts) != 3 {
		return 0, 0, 0, fmt.Errorf("invalid output format")
	}

	row, err := strconv.Atoi(parts[0])
	if err != nil {
		return 0, 0, 0, fmt.Errorf("invalid row number: %v", err)
	}

	col, err := strconv.Atoi(parts[1])
	if err != nil {
		return 0, 0, 0, fmt.Errorf("invalid column number: %v", err)
	}

	diff, err := strconv.ParseFloat(parts[2], 64)
	if err != nil {
		return 0, 0, 0, fmt.Errorf("invalid difference value: %v", err)
	}

	return row, col, diff, nil
}

func TestCheckSquareMatrixAndSameNumSpecies(t *testing.T) {
	inputDir := "Tests/CheckSquareMatrix/Input"
	outputDir := "Tests/CheckSquareMatrix/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			matrix, strains, err := readCheckSquareMatrixInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Read expected output (whether it should panic)
			shouldPanic, err := readCheckSquareMatrixOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Test the function
			testCheckSquareMatrix(t, matrix, strains, shouldPanic)
		})
	}
}

// should panic is the bool thats seen in the output file
func testCheckSquareMatrix(t *testing.T, matrix DistanceMatrix, strains []string, shouldPanic bool) {
	defer func() {
		r := recover()
		if shouldPanic && r == nil {
			t.Error("Expected panic but got none")
		}
		if !shouldPanic && r != nil {
			t.Errorf("Unexpected panic: %v", r)
		}
	}()

	CheckSquareMatrixAndSameNumSpecies(matrix, strains)
}

func readCheckSquareMatrixInput(filename string) (DistanceMatrix, []string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, err
	}
	defer file.Close()

	var matrix DistanceMatrix
	var strains []string
	scanner := bufio.NewScanner(file)
	section := ""

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		if strings.HasPrefix(line, "#") {
			section = strings.TrimSpace(strings.TrimPrefix(line, "#"))
			continue
		}

		switch section {
		case "Distance Matrix":
			var row []float64
			values := strings.Fields(line)
			for _, val := range values {
				num, err := strconv.ParseFloat(val, 64)
				if err != nil {
					return nil, nil, fmt.Errorf("invalid number format: %v", err)
				}
				row = append(row, num)
			}
			matrix = append(matrix, row)
		case "Genome Strains":
			strains = append(strains, line)
		}
	}

	return matrix, strains, nil
}

func readCheckSquareMatrixOutput(filename string) (bool, error) {
	file, err := os.Open(filename)
	if err != nil {
		return false, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	if !scanner.Scan() {
		return false, fmt.Errorf("empty output file")
	}

	line := strings.TrimSpace(scanner.Text())
	shouldPanic, err := strconv.ParseBool(line)
	if err != nil {
		return false, fmt.Errorf("invalid output format: %v", err)
	}

	return shouldPanic, nil
}

func TestInitializeTree(t *testing.T) {
	inputDir := "Tests/InitializeTree/Input"
	outputDir := "Tests/InitializeTree/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			genomeStrains, genomeNames, searchString, argName, percentIdentity, err := readInitializeTreeInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Initialize tree
			resultTree := InitializeTree(genomeStrains, genomeNames, searchString, argName, percentIdentity)

			// Read expected output
			expectedTree, err := readInitializeTreeOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Compare results
			if compareARGTreesInitialize(resultTree, expectedTree) {
				t.Errorf("InitializeTree result does not match expected output for file %s", inputFile)
				t.Errorf("Result Tree:\n%s", formatARGTree(resultTree))
				t.Errorf("Expected Tree:\n%s", formatARGTree(expectedTree))
			}
		})
	}
}

func readInitializeTreeInput(filename string) ([]string, []string, string, string, float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, nil, "", "", 0, err
	}
	defer file.Close()

	var genomeStrains, genomeNames []string
	var searchString, argName string
	var percentIdentity float64
	scanner := bufio.NewScanner(file)
	section := ""

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		if strings.HasPrefix(line, "#") {
			section = strings.TrimSpace(strings.TrimPrefix(line, "#"))
			continue
		}

		switch section {
		case "Genome Strains":
			genomeStrains = append(genomeStrains, line)
		case "Genome Names":
			genomeNames = append(genomeNames, line)
		case "Search String":
			searchString = line
		case "ARG Name":
			argName = line
		case "Percent Identity":
			percentIdentity, err = strconv.ParseFloat(line, 64)
			if err != nil {
				return nil, nil, "", "", 0, fmt.Errorf("invalid percent identity: %v", err)
			}
		}
	}

	return genomeStrains, genomeNames, searchString, argName, percentIdentity, nil
}

func readInitializeTreeOutput(filename string) (ARGTree, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var tree ARGTree
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		parts := strings.SplitN(line, ":", 2)
		if len(parts) != 2 {
			continue
		}

		nodeDetails := strings.TrimSpace(parts[1])
		node := &ARGNode{
			ARGpointer: &ARG{
				IndexInGenome: [2]int{-1, -1},
			},
		}

		// Split the details into sequence and indices
		fields := strings.Split(nodeDetails, "[")
		if len(fields) == 2 {
			// Get sequence
			sequence := strings.TrimSpace(fields[0])
			node.ARGpointer.Nucleotides = sequence

			// Get indices
			indicesStr := strings.Trim(fields[1], "[]")
			indices := strings.Split(indicesStr, " ")
			if len(indices) == 2 {
				start, _ := strconv.Atoi(indices[0])
				end, _ := strconv.Atoi(indices[1])
				node.ARGpointer.IndexInGenome = [2]int{start, end}
			}
		}

		tree = append(tree, node)
	}

	return tree, nil
}

func compareARGTreesInitialize(tree1, tree2 ARGTree) bool {
	if len(tree1) != len(tree2) {
		return false
	}

	for i := range tree1 {
		if !compareARGNodes(tree1[i], tree2[i]) {
			return false
		}
	}
	return true
}

func formatARGTreeInitialize(tree ARGTree) string {
	var sb strings.Builder
	for i, node := range tree {
		if node.ARGpointer != nil {
			if strings.Contains(node.ARGpointer.GenomeStrain, "Ancestor") {
				fmt.Fprintf(&sb, "Node %d: %s\n", i, node.ARGpointer.GenomeStrain)
			} else {
				fmt.Fprintf(&sb, "Node %d: %s %s [%d,%d] %s\n",
					i, node.ARGpointer.GenomeStrain,
					node.ARGpointer.Nucleotides,
					node.ARGpointer.IndexInGenome[0],
					node.ARGpointer.IndexInGenome[1],
					node.ARGpointer.ArgLabel)
			}
		}
	}
	return sb.String()
}

func TestAddRowCol(t *testing.T) {
	inputDir := "Tests/AddRowCol/Input"
	outputDir := "Tests/AddRowCol/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			matrix, row, col, clusterSize1, clusterSize2, err := readAddRowColInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Add row and column
			resultMatrix := AddRowCol(row, col, clusterSize1, clusterSize2, matrix)

			// Read expected output
			expectedMatrix, err := readAddRowColOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Compare results
			if !compareMatrices(resultMatrix, expectedMatrix) {
				t.Errorf("AddRowCol result does not match expected output for file %s", inputFile)
				t.Errorf("Result Matrix:\n%s", formatMatrix(resultMatrix))
				t.Errorf("Expected Matrix:\n%s", formatMatrix(expectedMatrix))
			}
		})
	}
}

func readAddRowColInput(filename string) (DistanceMatrix, int, int, int, int, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, 0, 0, 0, 0, err
	}
	defer file.Close()

	var matrix DistanceMatrix
	var row, col, clusterSize1, clusterSize2 int
	scanner := bufio.NewScanner(file)
	section := ""

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		if strings.HasPrefix(line, "#") {
			section = strings.TrimSpace(strings.TrimPrefix(line, "#"))
			continue
		}

		switch section {
		case "Matrix":
			var rowValues []float64
			values := strings.Fields(line)
			for _, val := range values {
				num, err := strconv.ParseFloat(val, 64)
				if err != nil {
					return nil, 0, 0, 0, 0, fmt.Errorf("invalid number format: %v", err)
				}
				rowValues = append(rowValues, num)
			}
			matrix = append(matrix, rowValues)
		case "Parameters":
			params := strings.Fields(line)
			if len(params) != 4 {
				return nil, 0, 0, 0, 0, fmt.Errorf("invalid parameters format")
			}
			row, _ = strconv.Atoi(params[0])
			col, _ = strconv.Atoi(params[1])
			clusterSize1, _ = strconv.Atoi(params[2])
			clusterSize2, _ = strconv.Atoi(params[3])
		}
	}

	return matrix, row, col, clusterSize1, clusterSize2, nil
}

func readAddRowColOutput(filename string) (DistanceMatrix, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var matrix DistanceMatrix
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		var rowValues []float64
		values := strings.Fields(line)
		for _, val := range values {
			num, err := strconv.ParseFloat(val, 64)
			if err != nil {
				return nil, fmt.Errorf("invalid number format: %v", err)
			}
			rowValues = append(rowValues, num)
		}
		matrix = append(matrix, rowValues)
	}

	return matrix, nil
}

func compareMatrices(matrix1, matrix2 DistanceMatrix) bool {
	if len(matrix1) != len(matrix2) {
		return false
	}

	for i := range matrix1 {
		if len(matrix1[i]) != len(matrix2[i]) {
			return false
		}
		for j := range matrix1[i] {
			if math.Abs(matrix1[i][j]-matrix2[i][j]) > 1e-6 {
				return false
			}
		}
	}
	return true
}

func formatMatrixAdd(matrix DistanceMatrix) string {
	var sb strings.Builder
	for _, row := range matrix {
		for _, val := range row {
			fmt.Fprintf(&sb, "%.6f ", val)
		}
		sb.WriteString("\n")
	}
	return sb.String()
}

func TestDeleteClustersARG(t *testing.T) {
	inputDir := "Tests/DeleteClustersARG/Input"
	outputDir := "Tests/DeleteClustersARG/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			clusters, row, col, err := readDeleteClustersInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Delete clusters
			resultClusters := DeleteClustersARG(clusters, row, col)

			// Read expected output
			expectedClusters, err := readDeleteClustersOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Compare results
			if !compareARGClusters(resultClusters, expectedClusters) {
				t.Errorf("DeleteClustersARG result does not match expected output for file %s", inputFile)
				t.Errorf("Result Clusters:\n%s", formatARGClusters(resultClusters))
				t.Errorf("Expected Clusters:\n%s", formatARGClusters(expectedClusters))
			}
		})
	}
}

func readDeleteClustersInput(filename string) ([]*ARGNode, int, int, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, 0, 0, err
	}
	defer file.Close()

	var clusters []*ARGNode
	var row, col int
	scanner := bufio.NewScanner(file)
	section := ""

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "//") {
			continue
		}

		if strings.HasPrefix(line, "#") {
			section = strings.TrimSpace(strings.TrimPrefix(line, "#"))
			continue
		}

		switch section {
		case "Clusters":
			parts := strings.Fields(line)
			if len(parts) >= 3 {
				node := &ARGNode{
					ARGpointer: &ARG{
						GenomeStrain: parts[0],
						Nucleotides:  parts[1],
					},
				}
				indicesStr := strings.Trim(parts[2], "[]")
				indices := strings.Split(indicesStr, ",")
				if len(indices) == 2 {
					start, _ := strconv.Atoi(indices[0])
					end, _ := strconv.Atoi(indices[1])
					node.ARGpointer.IndexInGenome = [2]int{start, end}
				}
				clusters = append(clusters, node)
			}
		case "Parameters":
			params := strings.Fields(line)
			if len(params) == 2 {
				row, _ = strconv.Atoi(params[0])
				col, _ = strconv.Atoi(params[1])
			}
		}
	}

	return clusters, row, col, nil
}

func readDeleteClustersOutput(filename string) ([]*ARGNode, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var clusters []*ARGNode
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		parts := strings.Fields(line)
		if len(parts) >= 3 {
			node := &ARGNode{
				ARGpointer: &ARG{
					GenomeStrain: parts[0],
					Nucleotides:  parts[1],
				},
			}
			indicesStr := strings.Trim(parts[2], "[]")
			indices := strings.Split(indicesStr, ",")
			if len(indices) == 2 {
				start, _ := strconv.Atoi(indices[0])
				end, _ := strconv.Atoi(indices[1])
				node.ARGpointer.IndexInGenome = [2]int{start, end}
			}
			clusters = append(clusters, node)
		}
	}

	return clusters, nil
}

func TestDeleteClustersProtein(t *testing.T) {
	inputDir := "Tests/DeleteClustersProtein/Input"
	outputDir := "Tests/DeleteClustersProtein/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			clusters, row, col, err := readDeleteClustersProteinInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Delete clusters
			resultClusters := DeleteClustersProtein(clusters, row, col)

			// Read expected output
			expectedClusters, err := readDeleteClustersProteinOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Compare results
			if !compareProteinClusters(resultClusters, expectedClusters) {
				t.Errorf("DeleteClustersProtein result does not match expected output for file %s", inputFile)
				t.Errorf("Result Clusters:\n%s", formatProteinClusters(resultClusters))
				t.Errorf("Expected Clusters:\n%s", formatProteinClusters(expectedClusters))
			}
		})
	}
}

func readDeleteClustersProteinInput(filename string) ([]*ProteinNode, int, int, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, 0, 0, err
	}
	defer file.Close()

	var clusters []*ProteinNode
	var row, col int
	scanner := bufio.NewScanner(file)
	section := ""

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "//") {
			continue
		}

		if strings.HasPrefix(line, "#") {
			section = strings.TrimSpace(strings.TrimPrefix(line, "#"))
			continue
		}

		switch section {
		case "Clusters":
			parts := strings.Fields(line)
			if len(parts) >= 3 {
				node := &ProteinNode{
					Protein: &Protein{
						AminoAcids: parts[0],
						ARG: &ARG{
							IndexInGenome: [2]int{-1, -1},
						},
					},
				}
				indicesStr := strings.Trim(parts[2], "[]")
				indices := strings.Split(indicesStr, ",")
				if len(indices) == 2 {
					start, _ := strconv.Atoi(indices[0])
					end, _ := strconv.Atoi(indices[1])
					node.Protein.ARG.IndexInGenome = [2]int{start, end}
				}
				clusters = append(clusters, node)
			}
		case "Parameters":
			params := strings.Fields(line)
			if len(params) == 2 {
				row, _ = strconv.Atoi(params[0])
				col, _ = strconv.Atoi(params[1])
			}
		}
	}

	return clusters, row, col, nil
}

func readDeleteClustersProteinOutput(filename string) ([]*ProteinNode, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var clusters []*ProteinNode
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		parts := strings.Fields(line)
		if len(parts) >= 3 {
			node := &ProteinNode{
				Protein: &Protein{
					AminoAcids: parts[0],
					ARG: &ARG{
						IndexInGenome: [2]int{-1, -1},
					},
				},
			}
			indicesStr := strings.Trim(parts[2], "[]")
			indices := strings.Split(indicesStr, ",")
			if len(indices) == 2 {
				start, _ := strconv.Atoi(indices[0])
				end, _ := strconv.Atoi(indices[1])
				node.Protein.ARG.IndexInGenome = [2]int{start, end}
			}
			clusters = append(clusters, node)
		}
	}

	return clusters, nil
}

func TestDeleteRowCol(t *testing.T) {
	inputDir := "Tests/DeleteRowCol/Input"
	outputDir := "Tests/DeleteRowCol/Output"

	inputFiles, err := filepath.Glob(filepath.Join(inputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find input files: %v", err)
	}

	outputFiles, err := filepath.Glob(filepath.Join(outputDir, "*.txt"))
	if err != nil {
		t.Fatalf("Failed to find output files: %v", err)
	}

	for i, inputFile := range inputFiles {
		testName := filepath.Base(inputFile)
		t.Run(testName, func(t *testing.T) {
			// Read input file
			matrix, row, col, err := readDeleteRowColInput(inputFile)
			if err != nil {
				t.Fatalf("Failed to read input file %s: %v", inputFile, err)
			}

			// Delete row and column
			resultMatrix := DeleteRowCol(matrix, row, col)

			// Read expected output
			expectedMatrix, err := readDeleteRowColOutput(outputFiles[i])
			if err != nil {
				t.Fatalf("Failed to read output file %s: %v", outputFiles[i], err)
			}

			// Compare results
			if !compareMatricesDelete(resultMatrix, expectedMatrix) {
				t.Errorf("DeleteRowCol result does not match expected output for file %s", inputFile)
				t.Errorf("Result Matrix:\n%s", formatMatrix(resultMatrix))
				t.Errorf("Expected Matrix:\n%s", formatMatrix(expectedMatrix))
			}
		})
	}
}

func readDeleteRowColInput(filename string) (DistanceMatrix, int, int, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, 0, 0, err
	}
	defer file.Close()

	var matrix DistanceMatrix
	var row, col int
	scanner := bufio.NewScanner(file)
	section := ""

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "//") {
			continue
		}

		if strings.HasPrefix(line, "#") {
			section = strings.TrimSpace(strings.TrimPrefix(line, "#"))
			continue
		}

		switch section {
		case "Matrix":
			var rowValues []float64
			values := strings.Fields(line)
			for _, val := range values {
				num, err := strconv.ParseFloat(val, 64)
				if err != nil {
					return nil, 0, 0, fmt.Errorf("invalid number format: %v", err)
				}
				rowValues = append(rowValues, num)
			}
			matrix = append(matrix, rowValues)
		case "Parameters":
			params := strings.Fields(line)
			if len(params) != 2 {
				return nil, 0, 0, fmt.Errorf("invalid parameters format")
			}
			row, _ = strconv.Atoi(params[0])
			col, _ = strconv.Atoi(params[1])
		}
	}

	return matrix, row, col, nil
}

func readDeleteRowColOutput(filename string) (DistanceMatrix, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	var matrix DistanceMatrix
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" {
			continue
		}

		var rowValues []float64
		values := strings.Fields(line)
		for _, val := range values {
			num, err := strconv.ParseFloat(val, 64)
			if err != nil {
				return nil, fmt.Errorf("invalid number format: %v", err)
			}
			rowValues = append(rowValues, num)
		}
		matrix = append(matrix, rowValues)
	}

	return matrix, nil
}

func compareMatricesDelete(matrix1, matrix2 DistanceMatrix) bool {
	if len(matrix1) != len(matrix2) {
		return false
	}

	for i := range matrix1 {
		if len(matrix1[i]) != len(matrix2[i]) {
			return false
		}
		for j := range matrix1[i] {
			if math.Abs(matrix1[i][j]-matrix2[i][j]) > 1e-6 {
				return false
			}
		}
	}
	return true
}

// TestTranslateProteinSequences tests the TranslateProteinSequences function (Shruthi)
func TestTranslateProteinSequences(t *testing.T) {
	tests, err := ReadProteinTranslationTests("Tests/TranslateProteinSequences/")
	if err != nil {
		t.Fatalf("Failed to read test cases: %v", err)
	}

	for _, test := range tests {
		result := TranslateProteinSequences(test.dnaFragments)
		// Check result
		if !SlicesAreEqual(result, test.result) {
			t.Errorf("TranslateProteinSequences(%v) = %v, want %v", test.dnaFragments, result, test.result)
		}
	}
}

// TestRNAtoAminoAcid tests the RNAtoAminoAcid function (Shruthi)
func TestRNAtoAminoAcid(t *testing.T) {
	tests, err := ReadRNAtoAminoAcid("Tests/RNAtoAminoAcid/")
	if err != nil {
		t.Fatalf("Failed to read test cases: %v", err)
	}

	for _, test := range tests {
		result := RNAtoAminoAcid(test.codon)
		// Check result
		if result != test.result {
			t.Errorf("RNAtoAminoAcid(%v) = %v, want %v", test.codon, result, test.result)
		}
	}
}

// SlicesAreEqual checks if two slices of strings are equal (Shruthi)
func SlicesAreEqual(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// TestCalculatePreservationOfNucleotides tests the CalculatePreservationOfNucleotides method (Shruthi)
func TestCalculatePreservationOfNucleotides(t *testing.T) {
	tests, err := ReadNucleotidePreservationTests("Tests/CalculatePreservationOfNucleotides/")
	if err != nil {
		t.Fatalf("Failed to read test cases: %v", err)
	}

	for _, test := range tests {
		test.nucleotideTree.CalculatePreservationOfNucleotides(test.wildtypeSeq)

		// Check preservation values
		for i, node := range test.nucleotideTree {
			if node.Preservation != test.expectedPreservation[i] {
				t.Errorf("Node %d Preservation = %v, want %v", i, node.Preservation, test.expectedPreservation[i])
			}
		}
	}
}

// TestCalculatePreservationOfProteins tests the CalculatePreservationOfProteins method (Shruthi)
func TestCalculatePreservationOfProteins(t *testing.T) {
	tests, err := ReadProteinPreservationTests("Tests/CalculatePreservationOfProteins/")
	if err != nil {
		t.Fatalf("Failed to read test cases: %v", err)
	}

	for _, test := range tests {
		test.proteinTree.CalculatePreservationOfProteins(test.wildtypeSeq)

		// Check preservation values
		for i, node := range test.proteinTree {
			if node.Preservation != test.expectedPreservation[i] {
				t.Errorf("Node %d Preservation = %v, want %v", i, node.Preservation, test.expectedPreservation[i])
			}
		}
	}
}

// TestCalculateLikelihood tests the CalculateLikelihood method (Shruthi)
func TestCalculateLikelihood(t *testing.T) {
	tests, err := ReadLikelihoodTests("Tests/CalculateLikelihood/")
	if err != nil {
		t.Fatalf("Failed to read test cases: %v", err)
	}

	for _, test := range tests {
		test.nucleotideTree.CalculateLikelihood(test.threshold)

		// Check likelihood values
		for i, node := range test.nucleotideTree {
			if node.Likelihood != test.expectedLikelihood[i] {
				t.Errorf("Node %d Likelihood = %v, want %v", i, node.Likelihood, test.expectedLikelihood[i])
			}
		}
	}
}

// ReadProteinTranslationTests takes as input a directory and returns a slice of ProteinTranslationTest objects (Shruthi)
func ReadProteinTranslationTests(directory string) ([]TranslateProteinSequencesTest, error) {
	inputFiles := ReadDirectory(directory + "/Input")
	var tests []TranslateProteinSequencesTest

	// Reading all input files
	for _, inputFile := range inputFiles {
		f, err := os.Open(directory + "Input/" + inputFile.Name())
		if err != nil {
			return nil, err
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		for scanner.Scan() {
			line := scanner.Text()
			parts := strings.Split(line, ",") // Assume input DNA fragments are comma-separated
			tests = append(tests, TranslateProteinSequencesTest{
				dnaFragments: parts,
			})
		}
	}

	// Reading all output files (expected results)
	outputFiles := ReadDirectory(directory + "/Output")
	if len(outputFiles) != len(tests) {
		return nil, fmt.Errorf("number of input and output files do not match")
	}

	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "Output/" + outputFile.Name())
		if err != nil {
			return nil, err
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		for scanner.Scan() {
			line := scanner.Text()
			parts := strings.Split(line, ",") // Assume output amino acids are comma-separated
			tests[i].result = parts
		}
	}

	return tests, nil
}

// ReadRNAtoAminoAcid takes as input a directory and returns a slice of RNAtoAminoAcidTest objects (Shruthi)
func ReadRNAtoAminoAcid(directory string) ([]RNAtoAminoAcidTest, error) {
	inputFiles := ReadDirectory(directory + "/Input")
	var tests []RNAtoAminoAcidTest

	// Reading all input files
	for _, inputFile := range inputFiles {
		f, err := os.Open(directory + "Input/" + inputFile.Name())
		if err != nil {
			return nil, err
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		for scanner.Scan() {
			line := scanner.Text()
			tests = append(tests, RNAtoAminoAcidTest{
				codon: line,
			})
		}
	}

	// Reading all output files (expected results)
	outputFiles := ReadDirectory(directory + "/Output")
	if len(outputFiles) != len(tests) {
		return nil, fmt.Errorf("number of input and output files do not match")
	}

	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "Output/" + outputFile.Name())
		if err != nil {
			return nil, err
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		for scanner.Scan() {
			line := scanner.Text()
			tests[i].result = line
		}
	}

	return tests, nil
}

// ReadNucleotidePreservationTests takes as input a directory and returns a slice of PreservationOfNucleotidesTest objects (Shruthi)
func ReadNucleotidePreservationTests(directory string) ([]PreservationOfNucleotidesTest, error) {
	inputFiles := ReadDirectory(directory + "/Input")
	var tests []PreservationOfNucleotidesTest

	// Reading all input files
	for _, inputFile := range inputFiles {
		f, err := os.Open(directory + "Input/" + inputFile.Name())
		if err != nil {
			return nil, err
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		var tree ARGTree
		var wildtypeSeq string
		lines := []string{}
		for scanner.Scan() {
			lines = append(lines, scanner.Text())
		}

		if len(lines) < 1 {
			return nil, fmt.Errorf("invalid input format: not enough lines")
		}

		// Wildtype sequence is the last line
		wildtypeSeq = lines[len(lines)-1]

		// Process all other lines as nucleotide sequences for nodes
		for _, line := range lines[:len(lines)-1] {
			node := &ARGNode{
				ARGpointer: &ARG{
					Nucleotides: line,
				},
				Preservation: 0.0, // Initial preservation value
			}
			tree = append(tree, node)
		}

		// Add to tests without expectedPreservation initially
		tests = append(tests, PreservationOfNucleotidesTest{
			nucleotideTree: tree,
			wildtypeSeq:    wildtypeSeq,
		})
	}

	// Reading all output files (expected preservation values)
	outputFiles := ReadDirectory(directory + "/Output")
	if len(outputFiles) != len(tests) {
		return nil, fmt.Errorf("number of input and output files do not match")
	}

	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "Output/" + outputFile.Name())
		if err != nil {
			return nil, err
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		var expected []float64
		for scanner.Scan() {
			line := scanner.Text()
			value, err := strconv.ParseFloat(line, 64)
			if err != nil {
				return nil, fmt.Errorf("failed to parse preservation value: %v", err)
			}
			expected = append(expected, value)
		}
		tests[i].expectedPreservation = expected
	}

	return tests, nil
}

// ReadProteinPreservationTests takes as input a directory and returns a slice of PreservationOfProteinsTest objects (Shruthi)
func ReadProteinPreservationTests(directory string) ([]PreservationOfProteinsTest, error) {
	inputFiles := ReadDirectory(directory + "/Input")
	var tests []PreservationOfProteinsTest

	// Reading all input files
	for _, inputFile := range inputFiles {
		f, err := os.Open(directory + "Input/" + inputFile.Name())
		if err != nil {
			return nil, err
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		var tree ProteinTree
		var wildtypeSeq string
		// var wildtypeLine string
		lines := []string{}
		for scanner.Scan() {
			lines = append(lines, scanner.Text())
		}

		if len(lines) < 1 {
			return nil, fmt.Errorf("invalid input format: not enough lines")
		}

		// Wildtype sequence is the last line
		wildtypeSeq = lines[len(lines)-1]

		// Process all other lines as nucleotide sequences for nodes
		for _, line := range lines[:len(lines)-1] {
			node := &ProteinNode{
				Protein: &Protein{
					AminoAcids: line,
				},
				Preservation: 0.0, // Initial preservation value
			}
			tree = append(tree, node)
		}

		// Add to tests without expectedPreservation initially
		tests = append(tests, PreservationOfProteinsTest{
			proteinTree: tree,
			wildtypeSeq: wildtypeSeq,
		})
	}

	// Reading all output files (expected preservation values)
	outputFiles := ReadDirectory(directory + "/Output")
	if len(outputFiles) != len(tests) {
		return nil, fmt.Errorf("number of input and output files do not match")
	}

	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "Output/" + outputFile.Name())
		if err != nil {
			return nil, err
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		var expected []float64
		for scanner.Scan() {
			line := scanner.Text()
			value, err := strconv.ParseFloat(line, 64)
			if err != nil {
				return nil, fmt.Errorf("failed to parse preservation value: %v", err)
			}
			expected = append(expected, value)
		}
		tests[i].expectedPreservation = expected
	}

	return tests, nil
}

// ReadLikelihoodTests takes as input a directory and returns a slice of LikelihoodTest objects (Shruthi)
func ReadLikelihoodTests(directory string) ([]CalculateLikelihoodTest, error) {
	inputFiles := ReadDirectory(directory + "/Input")
	var tests []CalculateLikelihoodTest

	// Reading all input files
	for _, inputFile := range inputFiles {
		f, err := os.Open(directory + "Input/" + inputFile.Name())
		if err != nil {
			return nil, err
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		var tree ARGTree
		var threshold float64
		for scanner.Scan() {
			line := scanner.Text()
			// Parse the ARGTree structure and threshold (assuming tab-separated format: Nucleotides\tPreservation\tThreshold)
			parts := strings.Split(line, " ")
			if len(parts) < 3 {
				return nil, fmt.Errorf("invalid input format: %s", line)
			}

			// Parse Nucleotides and Preservation
			preservation, err := strconv.ParseFloat(parts[1], 64)
			if err != nil {
				return nil, fmt.Errorf("failed to parse preservation value: %v", err)
			}
			node := &ARGNode{
				Preservation: preservation,
			}
			tree = append(tree, node)

			// Parse threshold if available
			threshold, err = strconv.ParseFloat(parts[2], 64)
			if err != nil {
				return nil, fmt.Errorf("failed to parse threshold value: %v", err)
			}
		}

		// Add to tests without expectedLikelihood initially
		tests = append(tests, CalculateLikelihoodTest{
			nucleotideTree: tree,
			threshold:      threshold,
		})
	}

	// Reading all output files (expected likelihood values)
	outputFiles := ReadDirectory(directory + "/Output")
	if len(outputFiles) != len(tests) {
		return nil, fmt.Errorf("number of input and output files do not match")
	}

	for i, outputFile := range outputFiles {
		f, err := os.Open(directory + "Output/" + outputFile.Name())
		if err != nil {
			return nil, err
		}
		defer f.Close()

		scanner := bufio.NewScanner(f)
		var expected []float64
		for scanner.Scan() {
			line := scanner.Text()
			value, err := strconv.ParseFloat(line, 64)
			if err != nil {
				return nil, fmt.Errorf("failed to parse likelihood value: %v", err)
			}
			expected = append(expected, value)
		}
		tests[i].expectedLikelihood = expected
	}

	return tests, nil
}
