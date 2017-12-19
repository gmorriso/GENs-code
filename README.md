# GENs-code
Code to compute the Generalized Erdos Numbers.  This code was written by Levi Dudte, implementing the algorithms described in 
Morrison G, Mahadevan L (2012) Discovering Communities through Friendship. PLOS ONE 7(7): e38704. https://doi.org/10.1371/journal.pone.0038704



Executive summary:

This source code will run the GENs community detection algorithm on a variety of file types.  The algorithm outputs the GENs as well as the community structure at varying levels of hierarchy.  

To run the algorithm on the provided example, run on the command line:
g++ ./*.cpp -O3 -o network
./network example.txt u

In general, to run the algorithm:

(1)  Store the network in a tab-delimited .txt file.
(2)  Cd to the code directory and type "g++ ./*.cpp -O3 -o network" (to generate an executable called "network")
(3)  Type "./network [network file] [network file type]" to analyze the network, where the network file type is: 'u' for an unweighted edges file, 'w' for a weighted edges file, 'a' for NxN weight matrix

Three file types will be output in the same directory as the network file: 
[network file]_GENs.txt (the closeness between nodes)
[network file]_OrderedEdges.txt (re-ordering the nodes into the detected community structure)
[network file]_CommunityTier_[i].txt (a set of files each corresponding to a tier of the hierarchical community structure: fine-grained communities @ i = 0, entire network @ i = # tiers - 1)  

Detailed summary:

This program is intended to detect community structure in networks using the closest unpopular friend algorithm with the Generalized Erdos Numbers as a measure of closeness.  The program requires a network to analyzed, which is stored in a plain text file.  The network is assumed to have N nodes and M edges, with the edges possibly weighted.  The network file must be formatted in one of two ways:

(a)  a sparse edge list of length M.  The correct format is node 1, node 2, and weight, tab delimited.  For an unweighted network, the weight value can be neglected.

(b) a NxN weight matrix.  The connection between each node is explicitly listed.  The weight between node 1 and every node in the network is listed in line 1, the connection between node 2 and every node in the network is listed in line 2, and so on. The first line must contain the number of nodes (N) in the NxN matrix.

EACH LINE OF THE FILES MUST BE TAB-DELIMITED.  Comma-separated files will not be read properly.

For large networks, the code may be relatively slow.  On a single core, we found a network of ~4000 nodes to take 1.5 hours to partition.

EXAMPLE (valid on Mac and Linux)

1)  Place each of the source files and the 'example.txt' (main.cpp, Network.h, NetworkIO.h) in the same directory.

2)  Cd to this directory and type "g++ ./*.cpp -O3 -o network"

	'g++' is the compiler, './*cpp' tells it to run on the files in your current directory, '-O3' optimizes the compilation so you have a faster executable, and '-o network' names the executable 'network'

3)  Type "./network example.txt u"

	'./network' calls the executable you just compiled, 'example.txt' contains a small network, and 'u' informs the program that 'example.txt' contains unweighted edges

When the code finishes you should see this output in the terminal:

Constructing theNetwork
Reading Unweighted Edges File [example.txt]
removeRedundantEdges()
organizeNodes()
initNodeToEdgesMap()
initConnectivity()
	# Nodes: 512
	# Edges: 7701
Computing GENs
	GENs -> 511/511
GENs Converged
Find Friends
Build FPs
Build Fine-Grained Communities: 17 Found
Repair Fractured Fine-Grained Communities: 16 Final
Build Coarse-Grained Communities:
	# communities at tier 0: 16
	# communities at tier 1: 4
	# communities at tier 2: 1
Writing Output Files

And you should see 5 new files in the directory:

example_GENs.txt (each row corresponds to a node and contains N entries, each of which is a GEN describing how close that column node feels to the row node)
example_OrderedEdges.txt (in the same format as example.txt, but the nodes have been renamed based on their order in the community hierarchy)
example_CommunityTier_0.txt (each row contains the nodes in a community)
example_CommunityTier_1.txt
example_CommunityTier_2.txt

Please report bugs / send comments or questions to ldudte@seas.harvard.edu and gcmorrison@uh.edu. Suggestions for improvement very welcome!


