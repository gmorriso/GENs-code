// includes
#include <iostream>
#include <fstream>
#include "Network.h"
#include "NetworkIO.h"

// globals
Network theNetwork;
string connectionsFile;
char fileType;

// network initialization
void initNetwork( const string& fileName, char fileType )
{	
	// declare structures
	vector < Network::Edge > tempEdges;
	vector < double > tempWeights;
	
	// read data and build theNetwork
	cout << "Constructing theNetwork" << endl;
	switch (fileType) {
		case 'u':
		{
			cout << "Reading Unweighted Edges File [" << fileName << "]" << endl;
			readEdgesUnweighted( fileName, tempEdges );
			theNetwork = Network( tempEdges, true );
			break;
		}
		case 'w':
		{
			cout << "Reading Weighted Edges File [" << fileName << "]" << endl;
			readEdgesWeighted( fileName, tempEdges, tempWeights );
			theNetwork = Network( tempEdges, tempWeights, true );
			break;
		}
		case 'a':
		{
			cout << "Reading Adjacency Matrix File [" << fileName << "]" << endl;
			readAdjacencyMatrix( fileName, tempEdges, tempWeights );
			theNetwork = Network( tempEdges, tempWeights, true );
			break;
		}
		default:
		{
			cout << "Bad File Type: " << fileType << endl;
			break;
		}
	}
}

// main
int main(int argc, char** argv)
{   
	// determine input file + type from command line parameters
	if(argc != 3)
	{
		cout << "Command line usage: [executable] [file name] [file type]" << endl;
		return 0;
	}
	else
	{
		stringstream ss;
		ss << argv[1];
		ss >> connectionsFile;
		fileType = argv[2][0];
	}

	// initialize network
	initNetwork(connectionsFile, fileType);
	
	// display network stats
	cout << "	# Nodes: " << theNetwork.getNumNodes() << endl;
	cout << "	# Edges: " << theNetwork.getNumEdges() << endl;
	
	// communities
	cout << "Computing GENs" << endl;
	theNetwork.computeGENs(theNetwork.getGENs());
	
	cout << "Find Friends" << endl;
	theNetwork.findCUFs(theNetwork.getGENs(), theNetwork.getFriends());
	
	cout << "Build FPs" << endl;
	theNetwork.buildFPs(theNetwork.getFriends(), theNetwork.getFriendPaths());
	
	cout << "Build Fine-Grained Communities: ";
	theNetwork.buildFineCommunities();
	cout << theNetwork.getNumCommunities(0) << " Found" << endl;
	
	cout << "Repair Fractured Fine-Grained Communities: ";
	theNetwork.repairFracturedCommunities(0);
	cout << theNetwork.getNumCommunities(0) << " Final" << endl;
	
	cout << "Build Coarse-Grained Communities:" << endl;
	theNetwork.buildCoarseCommunities(theNetwork.getGENs(), 0);
	int numTiers = theNetwork.getNumCommunityTiers();
	for(int i = 0; i < numTiers; i++)
		cout << "	# communities at tier " << i << ": " << theNetwork.getNumCommunities(i) << endl;

	cout << "Writing Output Files" << endl;
	string outputRoot = connectionsFile;
	outputRoot.erase(outputRoot.begin() + outputRoot.rfind("."), outputRoot.end());
	writeReorderedEdges(theNetwork, outputRoot + "_OrderedEdges.txt");
	writeCommunities(theNetwork, outputRoot + "_CommunityTier");
	writeGENs(outputRoot + "_GENs.txt", theNetwork.getGENs());
}
