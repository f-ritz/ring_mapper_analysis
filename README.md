# Ring Mapper Analysis

This is a ring mapper algorithm.

Compiling instructions:
- Download the zip
- Find and enter the directory "ring_mapper_analysis"
- "cmake -G Ninja"
- "ninja"
- To run the program follow the instructions below.

Input arguments in the command line are as follows:
"./ring_mapper_analysis directory_of_your_test_file src/unittest_resurces/structure_input.txt number_of_reads i_to_retrieve j_to_retrieve"

The directory and the number of reads are required for the program to function.
The structure input directory should remain as is, only the text file itself should be edited.
In structure_input.txt, the first line should always be the ATCG chain, while the second should be the dot-bracket notation.
You have to manually input ATCG chain and dot-bracket notation.
If you don't want to retrieve anything just put random numbers for i and j, it won't affect the rest of the program.

Example command:
./ring_mapper_analysis src/unittest_resources/minittr-6-2HP-ref_bitvectors.txt src/unittest_resurces/structure_input.txt 1000 8 16

ring_mapper_results.csv delivers the raw results
ring_mapper_histogram.csv prints a histogram
ring_mapper_retrieval.csv retrieves the data in the specified position
ring_mapper_counter.csv displays the number of mutations in each position
ring_mapper_pairmap.csv shows the types of mutated basepairs

ring_mapper_analysis.csv shows everything except the histogram and retrieval

everything else is debugging rubbish