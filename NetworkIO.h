////////////////////////////////////////////////////////////////////////////////
// 
//	Written by Levi Dudte, Mahadevan Group @ Harvard, Spring 2012
//
////////////////////////////////////////////////////////////////////////////////

// includes
#include "Network.h"
#include <sstream>

// return a string of an int (for output files)
string convertInt(int number)
{
	ostringstream ss; //create a stringstream
	ss << number; //add number to the stream
	return ss.str(); //return a string with the contents of the stream
}

// read unweighted edges file
void readEdgesUnweighted( const string& fileName, vector< Network::Edge >& edges )
{
	// open the file stream
	ifstream connectionsStream;
	connectionsStream.open( fileName.c_str() );
	if (!connectionsStream) 
	{
		cout << "Unable to open connections file." << endl;
		exit(1);
	}
		
	// read in connections
	while( true )
	{
		// read in value
		Network::Edge tempEdge;
		connectionsStream >> tempEdge.first;
		
		// check for eof
		if( connectionsStream.eof() )
			break;
		// add a new edge
		else
		{
			connectionsStream >> tempEdge.second;
			edges.push_back( tempEdge );
		}
	}
		
	// close connections file stream
	connectionsStream.close();
}

// read unweighted edges file
void readEdgesWeighted( const string& fileName, vector< Network::Edge >& edges, vector < double >& weights )
{
	// open the file stream
	ifstream connectionsStream;
	connectionsStream.open( fileName.c_str() );
	if (!connectionsStream) 
	{
		cout << "Unable to open connections file." << endl;
		exit(1);
	}
	
	// read in connections
	while( true )
	{
		// read in value
		Network::Edge tempEdge;
		double tempWeight;
		connectionsStream >> tempEdge.first;
		
		// check for eof
		if( connectionsStream.eof() )
			break;
		// add a new edge
		else
		{
			connectionsStream >> tempEdge.second;
			connectionsStream >> tempWeight;
			edges.push_back( tempEdge );
			weights.push_back( tempWeight );
		}
	}
	
	// close connections file stream
	connectionsStream.close();
}

// read a weighted adjacency matrix
void readAdjacencyMatrix( const string& fileName, vector< Network::Edge >& edges, vector < double >& weights )
{
	// open the file stream
	ifstream connectionsStream;
	connectionsStream.open( fileName.c_str() );
	if (!connectionsStream) 
	{
		cout << "Unable to open connections file." << endl;
		exit(1);
	}
	
	// get number of nodes
	int numNodes;
	connectionsStream >> numNodes;
	int numEntries = (int)pow(numNodes, 2.0);
	
	// read in connections
	vector < double > adjMatrix;
	for(int i = 0; i < numEntries; i++)
	{
		// read in value
		double tempValue;
		connectionsStream >> tempValue;
		
		// add it to the matrix
		adjMatrix.push_back(tempValue);
	}
	
	// get the edges
	double epsilon = .0000001;
	for(int i = 0; i < numNodes - 1; i++)
		for(int j = i + 1; j < numNodes; j++)
		{
			// get the indices into adjMatrix
			int thisIndex = i*numNodes + j, symmetricIndex = j*numNodes + i;
			
			// make sure the values are within an epsilon 
			// and that both are non-zero
			if(abs(adjMatrix[thisIndex] - adjMatrix[symmetricIndex]) > epsilon
			   || abs(adjMatrix[thisIndex]) < epsilon
			   || abs(adjMatrix[symmetricIndex]) < epsilon)
				continue;
			else
			{
				// add this edge
				Network::Edge newEdge(i, j);
				double newWeight = adjMatrix[thisIndex];
				edges.push_back(newEdge);
				weights.push_back(newWeight);
			}
		}
	
	// close connections file stream
	connectionsStream.close();
}

// write the GENs
void writeGENs(const string& fileName, const vector < vector < double > >& GENs)
{
	// don't write anything if GENs is not square
	int numRows = GENs.size();
	for(int i = 0; i < numRows; i++)
		if(GENs[i].size() != numRows)
		{
			cout << "Invalid GENs dimensions in writeGENs." << endl;
			return;
		}

	// write the values
	ofstream outStream;
	outStream.open( fileName.c_str() );
	for(int i = 0; i < numRows; i++)
	{
		for(int j = 0; j < numRows; j++)
			outStream << GENs[i][j] << " ";
		outStream << endl;
	}
	outStream.close();
}

// write the edges to display as an adj matrix
void writeReorderedEdges(Network& theNetwork, const string& fileName)
{
	// get stats
	int numNodes = theNetwork.getNumNodes();
	int numEdges = theNetwork.getNumEdges();
	int numTiers = theNetwork.getNumCommunityTiers();
	
	// get the nodes in the order they appear from the top of the community hierarchy
	vector < int > orderedNodes = theNetwork.getNodesInCommunity(0, numTiers - 1);
	vector < int > newIndices(numNodes, -1);
	for(int i = 0; i < numNodes; i++)
		newIndices[orderedNodes[i]] = i;
	
	// compute scores 
	int origScore = 0, newScore = 0;
	for(int i = 0; i < numEdges; i++)
	{
		origScore += abs(theNetwork.getEdge(i).first - theNetwork.getEdge(i).second);
		newScore += abs(newIndices[theNetwork.getEdge(i).first] - newIndices[theNetwork.getEdge(i).second]);
	}	
	
	// now write the edges, replacing each node index with its effective index (entry in effectiveIndices)
	ofstream outStream;
	outStream.open(fileName.c_str());
	for(int i = 0; i < numEdges; i++)
	{
		int nodeA = newIndices[theNetwork.getEdge(i).first], nodeB = newIndices[theNetwork.getEdge(i).second];
		Network::Edge newEdge(nodeA < nodeB ? nodeA : nodeB, nodeA < nodeB ? nodeB : nodeA);
		outStream << newEdge.first << " " << newEdge.second << endl;
	}	
	outStream.close();
}

// write the communities
void writeCommunities(Network& theNetwork, const string& fileNameRoot)
{
	// get stats
	int numTiers = theNetwork.getNumCommunityTiers();
	
	// write each community
	int numCommunities;
	for(int tier = 0; tier < numTiers; tier++)
	{
		numCommunities = theNetwork.getNumCommunities(tier);
		ofstream outStream;
		outStream.open((fileNameRoot + "_" + convertInt(tier) + ".txt").c_str());
		for(int i = 0; i < numCommunities; i++)
		{
			// get the nodes in this community
			vector < int > thisCommunity = theNetwork.getNodesInCommunity(i, tier);
			for(int j = 0; j < thisCommunity.size(); j++)
				outStream << thisCommunity[j] << " ";
			outStream << endl;
		}
		outStream.close();
	}
}
