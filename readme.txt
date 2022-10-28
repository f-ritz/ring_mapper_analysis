This is a ring mapper algorithm.

Input arguments in the command line are as follows:

"./ring_mapper_analysis <directory of your test file> <number of reads> <i to retrieve> <j to retrieve>"

The directory and the number of reads are required for the program to function.
The last two arguments (i and j) can be omitted if you don't want to retrieve anything.

ring_mapper_results.csv delivers the raw results
ring_mapper_histogram.csv prints a histogram
ring_mapper_retrieval.csv retrieves the data in the specified position (not currently functional)
everything else is debugging rubbish

*** note: just put 0 for i and j each for now as it's currently segfaulting otherwise, the program will function just fine
*** note_2: the retrieval algorithm is not currently working properly