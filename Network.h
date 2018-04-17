////////////////////////////////////////////////////////////////////////////////
// 
//	Written by Levi Dudte, Mahadevan Group @ Harvard, Spring 2012
//
////////////////////////////////////////////////////////////////////////////////

// prevent class redefinition
#ifndef NETWORK_H
#define NETWORK_H

// includes
#include <vector>
#include <iostream>
#include <cmath>
#include <stdlib.h>

// namespace
using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Network
////////////////////////////////////////////////////////////////////////////////

class Network
{
	public:
		// define type of pair: Edge
		typedef pair< int, int > Edge;

		// constructors
		Network( bool createEdges = false );
		Network( const vector< Network::Edge >& edges, bool createEdges );
		Network( const vector< Network::Edge >& edges, const vector< double >& weights, bool createEdges );

		// getters (network data)
		int getNumNodes() const;
		int getNumEdges() const;
		int getDegree(int nodeIndex) const;
		bool sameConnectedComponent(int nodeA, int nodeB) const;
		double getGEN(int i, int j) const;
		int getNumCommunities(int tier) const;
		int getNumCommunityTiers() const;
	
		// getters (references / network structures)
		vector < vector < double > >& getGENs();
		vector < int >& getFriends();
		vector < vector < int > >& getFriendPaths();
		vector < vector < int > >& getCommunities(int tier);
		vector < int >& getCommunityIDs(int tier);
		Edge getEdge(int edgeIndex);
		
		// initialize additional topological information
		void initNodeToEdgesMap();
		void initConnectivity();
		
		// organize network
		void removeRedundantEdges();
		void organizeNodes();
		void renameNodesInEdges();
		
		// weights
		void initEdgeWeights();
		
		// compute GENs (conventional or modified)
		void computeGENs(vector < vector < double > >& GENs);
		void computeModifiedGENs(vector < vector < double > >& GENs);
	
		// random walking
		void seedRandomWalk();
		int randomInt(int OutOf);
		int getNextRandomNode(int nodeIndex);
		int takeRandomWalk(int startNode, int destNode);
		double MFPT(int startNode, int destNode, int numWalks);
	
		// closest friends
		void findCFs(const vector < vector < double > >& GENs, vector < int >& CFs);
		void findCUFs(const vector < vector < double > >& GENs, vector < int >& CUFs);
		void buildFPs(const vector < int >& friends, vector < vector < int > >& friendPaths);
	
		// community detection
		void buildCommunities(const vector < vector < int > >& friendPaths, 
							  vector < vector < int > >& communities,
							  vector < int >& communityIDs);
		void buildFineCommunities();
		void repairFracturedCommunities(int tier);
		void buildCoarseCommunities(const vector < vector < double > >& GENs, int thisTier);
		vector < int > getNodesInCommunity(int communityID, int tier);
		int getCommunityOfNode(int nodeIndex, int tier);
	
		// output
		void writeGENs(const string& fileName, const vector < vector < double > >& GENs);
		void writeReorderedEdges(const string& fileName);
		void writeCommunities(const string& fileNameRoot);
		
	private:
		// nodes
		vector < int > m_originalNodes;
	
		// topology
		bool m_createEdges;
		vector < Edge > m_edges;
		vector < vector < int > > m_nodeToEdgesMap;
		vector < int > m_connectivity;
		
		// weights
		vector < double > m_edgeWeights;
		
		// Generalized Erdos Numbers (conventional or modified)
		vector < vector< double > > m_GENs;

		// closest friends (fine grained)
		vector < int > m_friends;
		vector < vector< int > > m_friendPaths;
		
		// communities
		vector < vector < int > > m_communityIDs;
		vector < vector < vector < int > > > m_communities;
	
		// hierarchical coarse-grained communities
		vector < vector < int > > m_coarseCommunityIDs;
		vector < vector < vector < int > > > m_coarseCommunities;
	
		// booleans (tracks which definitions were used to compute various network structures)
		bool m_usedConventionalGENs, m_usedModifiedGENs, m_usedCFs, m_usedCUFs;
		
};

////////////////////////////////////////////////////////////////////////////////
// Constructors
////////////////////////////////////////////////////////////////////////////////

Network::Network( bool createEdges )
:
m_createEdges( createEdges )
{
}

Network::Network( const vector< Network::Edge >& edges, bool createEdges )
:
m_createEdges( createEdges ),m_edges( edges )
{
	removeRedundantEdges();
	organizeNodes();
	renameNodesInEdges();
	initNodeToEdgesMap();
	initConnectivity();
	initEdgeWeights();
}
Network::Network( const vector< Network::Edge >& edges, const vector< double >& weights, bool createEdges )
:
m_createEdges( createEdges ),m_edges( edges ),m_edgeWeights( weights )
{
	removeRedundantEdges();
	organizeNodes();
	renameNodesInEdges();
	initNodeToEdgesMap();
	initConnectivity();
}


////////////////////////////////////////////////////////////////////////////////
// Getters
////////////////////////////////////////////////////////////////////////////////

int Network::getNumNodes() const{return m_originalNodes.size();}
int Network::getNumEdges() const{return m_edges.size();}
int Network::getDegree(int nodeIndex) const{return m_nodeToEdgesMap[nodeIndex].size();}
double Network::getGEN(int i, int j) const{return m_GENs[i][j];}
int Network::getNumCommunities(int tier) const{return m_communities[tier].size();}
int Network::getNumCommunityTiers() const{return m_communities.size();}

vector < vector < double > >& Network::getGENs(){return m_GENs;}
vector < int >& Network::getFriends(){return m_friends;}
vector < vector < int > >& Network::getFriendPaths(){return m_friendPaths;}
vector < vector < int > >& Network::getCommunities(int tier){return m_communities[tier];}
vector < int >& Network::getCommunityIDs(int tier){return m_communityIDs[tier];}
Network::Edge Network::getEdge(int edgeIndex){return m_edges[edgeIndex];}

////////////////////////////////////////////////////////////////////////////////
// Additional Topology
////////////////////////////////////////////////////////////////////////////////

// initialize m_vertToEdgesMap
void Network::initNodeToEdgesMap()
{
	// inform user
	cout << "initNodeToEdgesMap()" << endl;
	
	// nums
	int numNodes = getNumNodes();
	int numEdges = getNumEdges();
	vector < int > empty;
	
	// put numVerts empty vectors into these maps
	for(int i = 0; i < numNodes; i++)
		m_nodeToEdgesMap.push_back(empty);
	
	// fill m_vertToEdgesMap
	for(int i = 0; i < numEdges; i++)
	{
		m_nodeToEdgesMap[m_edges[i].first].push_back(i);
		m_nodeToEdgesMap[m_edges[i].second].push_back(i);
	}
}

// determine connectivity information
void Network::initConnectivity()
{
	// inform user
	cout << "initConnectivity()" << endl;
	
	// number of nodes
	int numNodes = getNumNodes();
	
	// initialize m_connectivity
	for(int i = 0; i < numNodes; i++)
		m_connectivity.push_back(-1);
	
	// iterate through connectivity groups
	int remainingNodes = numNodes;
	for(int i = 0; remainingNodes > 0; i++)
	{
		// structures needed to build this disjoint subset of the graph
		vector< int > thisWave, nextWave;
		
		// set thisWave to the first unfound node index
		for(int j = 0; j < numNodes; j++)
			if(m_connectivity[j] == -1)
			{
				thisWave.push_back(j);
				break;
			}
		
		// add new waves until a nextWave contributes no new nodes to m_connectivity
		int numNewNodes=0;//, totalFound = 0;
		do
		{
			// store all nodes connected to thisWave in nextWave
			for(int j = 0; j < thisWave.size(); j++)
			{
				// iterate through m_nodeToEdgesMap[thisWave[j]]
				for(int k = 0; k <  m_nodeToEdgesMap[thisWave[j]].size(); k++)
				{
					// add the other node index to next wave
					if(m_edges[m_nodeToEdgesMap[thisWave[j]][k]].first == thisWave[j] && 
					   m_connectivity[m_edges[m_nodeToEdgesMap[thisWave[j]][k]].second] < 0)
						nextWave.push_back(m_edges[m_nodeToEdgesMap[thisWave[j]][k]].second);
					else if(m_edges[m_nodeToEdgesMap[thisWave[j]][k]].second == thisWave[j] && 
							m_connectivity[m_edges[m_nodeToEdgesMap[thisWave[j]][k]].first] < 0)
						nextWave.push_back(m_edges[m_nodeToEdgesMap[thisWave[j]][k]].first);
				}
			}
			
			// set the nodes in nextWave to i (disjoint subset label) in m_connectivity
			numNewNodes = 0;
			for(int j = 0; j < nextWave.size(); j++)
				if(m_connectivity[nextWave[j]] < 0) // an unfound node
				{
					m_connectivity[nextWave[j]] = i;
					numNewNodes++;
				}
			
			// set thisWave to nextWave
			thisWave = nextWave;
			nextWave.clear();
		}
		while(numNewNodes > 0);
		
		// count unfound nodes
		remainingNodes = 0;
		for(int j = 0; j < numNodes; j++)
			if(m_connectivity[j] < 0)
				remainingNodes++;
	}
}

// return true if nodes are in the same disjoint subset
bool Network::sameConnectedComponent(int nodeA, int nodeB) const
{
	if(m_connectivity[nodeA] == m_connectivity[nodeB])
		return true;
	else
		return false;
}

////////////////////////////////////////////////////////////////////////////////
// Clean Up
////////////////////////////////////////////////////////////////////////////////

void Network::removeRedundantEdges()
{	
	// inform user
	cout << "removeRedundantEdges()" << endl;
	
	// remove any duplicate edges in m_edges
	//int numRemoved = 0;
	for(int i = 0; i < m_edges.size(); i++)
	{
		// get this edge
		Network::Edge thisEdge = m_edges[i];
		if(thisEdge.first == thisEdge.second) 
		{
			m_edges.erase(m_edges.begin() + i); // erase edges with the same node for each index
			i--; // shift i down one for next iteration
		}
		
		else
		{
			// compare it with every other edge in m_edges
			for(int j = 0; j < m_edges.size(); j++)
			{
				if(i != j) // don't check this edge against itself
				{
					if( (thisEdge.first == m_edges[j].first && thisEdge.second == m_edges[j].second) ||
					   (thisEdge.first == m_edges[j].second && thisEdge.second == m_edges[j].first) )
					{
						m_edges.erase(m_edges.begin() + j); // erase this edge
						j--; // shift j down for next iteration
						if(j < i) 
							i--; // shift i down for next iteration if j was beneath it in m_edges
					}
				}
			}
		}	
	}
}

// fixes holes in the indices count
void Network::organizeNodes()
{
	// inform user
	cout << "organizeNodes()" << endl;
	
	// num edges
	int numEdges = getNumEdges();
	
	// each unique node
	vector < int > nodes;
	
	// fill nodes with each unique index drawn from raw m_edges
	for(int i = 0; i < numEdges; i++)
	{
		// get the nodes of this edge
		int firstNode = m_edges[i].first;
		int secondNode = m_edges[i].second;
		
		// check all preceding nodes against the two current nodes
		bool haveFirstNode = false;
		bool haveSecondNode = false;
		for(int j = 0; j < nodes.size(); j++)
		{
			if(nodes[j] == firstNode)
				haveFirstNode = true;
			if(nodes[j] == secondNode)
				haveSecondNode = true;
		}
		
		// add these two nodes, if not found
		if(!haveFirstNode)
			nodes.push_back(firstNode);
		if(!haveSecondNode)
			nodes.push_back(secondNode);
	}
	
	// sort these nodes (selection)
	int numNodes = nodes.size();
	for(int i = 0; i < numNodes - 1; i++)
	{	
		// find the smallest index in the remaining subset
		int smallest = i;
		for(int j = i+1; j < numNodes; j++)
			if(nodes[j] < nodes[smallest])
				smallest = j;
		int swap = nodes[smallest];
		nodes[smallest] = nodes[i];
		nodes[i] = swap;
	}
	
	// set m_originalNodes to nodes (sorted)
	// this allows us to rename all of the nodes internally in case of any gaps in the node indicies from the data
	m_originalNodes = nodes;
}

// rename the node indices in the edges for each node index that doesn't match its position in m_originalNodes
void Network::renameNodesInEdges()
{
	// nums
	int numNodes = getNumNodes();
	int numEdges = getNumEdges();
	
	// check each index in m_originalNodes
	for(int i = 0; i < numNodes; i++)
		if(m_originalNodes[i] != i)
			// go through m_edges, renaming each instance of the problem node
			for(int j = 0; j < numEdges; j++)
			{
				if(m_edges[j].first == m_originalNodes[i])
					m_edges[j].first = i;
				if(m_edges[j].second == m_originalNodes[i])
					m_edges[j].second = i;
			}
}

////////////////////////////////////////////////////////////////////////////////
// Edge Weights
////////////////////////////////////////////////////////////////////////////////
void Network::initEdgeWeights()
{	
	// set edge weights to 1
	int numEdges = getNumEdges();
	for(int i = 0; i < numEdges; i++)
		m_edgeWeights.push_back(1.0);
}

////////////////////////////////////////////////////////////////////////////////
// GENs
////////////////////////////////////////////////////////////////////////////////

// iteratively compute GENs
void Network::computeGENs(vector < vector < double > >& GENs)
{	
	// num nodes
	int numNodes = getNumNodes();
	
	// compute W for each node (its strength: the sum of its incident edge weights)
	vector< double > nodeWs;
	for(int i = 0; i < numNodes; i++)
	{
		double thisNodeW = 0.0;
		for(int j = 0; j < (int)m_nodeToEdgesMap[i].size(); j++)
			thisNodeW += m_edgeWeights[m_nodeToEdgesMap[i][j]];
		nodeWs.push_back(thisNodeW);
	}
	
	// define a GEN convergence threshold
	double convergenceThreshold = .005;
	
	// iteratively compute GENs
	int allHaveConverged, maxGENcount = 6000; // track iterations
	
	// GENs(t)_ij are given by GENs(t-1)_ij
	for(int i = 0; i < numNodes; i++) // i: node index
	{	
		// allocate temporary storage for GENs iterations
		vector< double > currentGENs;
		vector< double > nextGENs;
		for(int j = 0; j < numNodes; j++)
		{
			if(i != j)
			{
				currentGENs.push_back(1.0);
				nextGENs.push_back(1.0);
			}
			else
			{
				currentGENs.push_back(0.0);
				nextGENs.push_back(0.0);
			}
		}
		
		// compute this row of GENs
		for(int GENcounter = 0; GENcounter < maxGENcount; GENcounter++)
		{
			// check for convergence in this iteration
			allHaveConverged = 1;
			
			// update GEN row from every other node in the network to node i
			for(int j = 0; j < numNodes; j++) // j: node index
			{
				// assign -1 GEN to nodes in disjoint subsets of the network
				if(m_connectivity[i] != m_connectivity[j])
					nextGENs[j] = -1;
				
				// update all GENs
				else if(i != j)
				{
					// number of incident edges at this node
					int localDegree = m_nodeToEdgesMap[j].size();
					
					// compute the summation term in the RHS of GEN definition
					double flowerSum = 0.0;
					for(int k = 0; k < localDegree; k++) // k: m_nodeToEdgesMap[j] index
					{
						// edge index
						int thisEdgeIndex = m_nodeToEdgesMap[j][k];
						
						// determine the index of the node connected to node i
						int otherNodeIndex=0;
						if(m_edges[thisEdgeIndex].first == j)
							otherNodeIndex = m_edges[thisEdgeIndex].second;
						else
							otherNodeIndex = m_edges[thisEdgeIndex].first;
						
						// add this node's contribution to flowerSum
						flowerSum += m_edgeWeights[thisEdgeIndex]/(currentGENs[otherNodeIndex] + 1.0/m_edgeWeights[thisEdgeIndex]);	
					}
					
					// scale the flowerSum by 1/nodeWs[i]
					nextGENs[j] = nodeWs[j]/flowerSum;
					
					// check for convergence
					if(abs(currentGENs[j] - nextGENs[j]) < convergenceThreshold)
						allHaveConverged *= 1; // will remain "true" if all GENs have converged below the threshold
					else
						allHaveConverged *= 0; // will switch to "false" if any GEN has not converged beneath the threshold*/
				}
			}
			
			// set currentGENs to nextGENs
			currentGENs = nextGENs;
			
			// break if all GENs have converged to convergenceThreshold
			if(allHaveConverged)
				break;
		}
		
		// store the temp GENs
		GENs.push_back(currentGENs);
		
		// inform user
		cout << "	GENs -> " << i << "/" << numNodes - 1 << "\r" << flush;
	}
	
	// inform user of convered GENs
	cout << endl << "GENs Converged" << endl;
	
	// track GEN definition
	m_usedConventionalGENs = true;
	m_usedModifiedGENs = false;
}

// compute the modified definition of the GENs (for comparision to the original definition)
void Network::computeModifiedGENs(vector < vector < double > >& GENs)
{
	// num nodes
	int numNodes = getNumNodes();
	
	// compute W for each node (its strength: the sum of its incident edge weights)
	vector< double > nodeWs;
	for(int i = 0; i < numNodes; i++)
	{
		double thisNodeW = 0.0;
		for(int j = 0; j < (int)m_nodeToEdgesMap[i].size(); j++)
			thisNodeW += m_edgeWeights[m_nodeToEdgesMap[i][j]];
		nodeWs.push_back(thisNodeW);
	}
	
	// allocate temporary storage for GENs iterations
	vector < vector < double > > currentGENs;
	vector < vector < double > > nextGENs;
	for(int i = 0; i < numNodes; i++)
	{
		// make new row
		vector < double > emptyRow;
		currentGENs.push_back(emptyRow);
		nextGENs.push_back(emptyRow);
		
		// fill row
		for(int j = 0; j < numNodes; j++)
		{
			if(i != j)
			{
				currentGENs[i].push_back(1.0);
				nextGENs[i].push_back(1.0);
			}
			else
			{
				currentGENs[i].push_back(0.0);
				nextGENs[i].push_back(0.0);
			}
		}
	}
	
	// define a GEN convergence threshold
	double convergenceThreshold = .005;
	
	// iteratively compute GENs
	int allHaveConverged, maxGENcount = 2000; // track iterations
	
	// compute new iteration of GENs
	for(int GENcounter = 0; GENcounter < maxGENcount; GENcounter++)
	{
		// check for convergence in this iteration
		allHaveConverged = 1;
		
		// iterate the entire network of GENs
		for(int i = 0; i < numNodes; i++) // i: node index
		{
			// update GEN row from every other node in the network to node i
			for(int j = 0; j < numNodes; j++) // j: node index
			{
				// assign -1 GEN to nodes in disjoint subsets of the network
				if(!sameConnectedComponent(i, j))
					nextGENs[i][j] = -1;
				
				// update all GENs
				else if(i != j)
				{
					// number of incident edges at this node
					int localDegree = m_nodeToEdgesMap[j].size();
					
					// compute the summation term in the RHS of GEN definition
					double flowerSum = 0.0, extraWeight;
					bool addExtraWeight = false;
					for(int k = 0; k < localDegree; k++) // k: m_nodeToEdgesMap[j] index
					{
						// edge index
						int thisEdgeIndex = m_nodeToEdgesMap[j][k];
						
						// determine the index of the node connected to node i
						int otherNodeIndex;
						if(m_edges[thisEdgeIndex].first == j)
							otherNodeIndex = m_edges[thisEdgeIndex].second;
						else
							otherNodeIndex = m_edges[thisEdgeIndex].first;
						
						// add this node's contribution to flowerSum
						if(i != otherNodeIndex)
							flowerSum += m_edgeWeights[thisEdgeIndex]/(currentGENs[i][otherNodeIndex] + currentGENs[otherNodeIndex][j]);	
						else
						{
							addExtraWeight = true;
							extraWeight = pow(m_edgeWeights[thisEdgeIndex], 2.0);
						}
					}
					
					// add the extra weight if i was in the flower of j
					if(addExtraWeight)
						nextGENs[i][j] = nodeWs[j]/(flowerSum + extraWeight);
					
					// scale the flowerSum by 1/nodeWs[i]
					else
						nextGENs[i][j] = nodeWs[j]/flowerSum;
					
					// check for convergence
					if(abs(currentGENs[i][j] - nextGENs[i][j]) < convergenceThreshold)
						allHaveConverged *= 1; // will remain "true" if all GENs have converged below the threshold
					else
						allHaveConverged *= 0; // will switch to "false" if any GEN has not converged beneath the threshold*/
				}
			}
		}
		
		// set currentGENs to nextGENs
		currentGENs = nextGENs;
		
		// inform user
		cout << "	# Iterations -> " << GENcounter + 1 << "\r" << flush;
		
		// break if all GENs have converged to convergenceThreshold
		if(allHaveConverged)
		{
			cout << "Modified GENs Converged (" << GENcounter + 1 << " Iterations)" << endl;
			break;
		}
	}
	
	// store the temp GENs
	GENs = currentGENs;
	
	// track GEN definition
	m_usedConventionalGENs = false;
	m_usedModifiedGENs = true;
}

////////////////////////////////////////////////////////////////////////////////
// Random Walking
////////////////////////////////////////////////////////////////////////////////

void Network::seedRandomWalk()
{
	// seed for random walking
	srand(time(NULL));
}
 
// return a perfectly uniform random int in [0, OutOf - 1]
inline
int Network::randomInt(int OutOf)
{
	// compute out of bounds
	int leftOver = RAND_MAX % OutOf;
	
	// compute random int
	int random;
	do{random = rand();}while(random >= RAND_MAX - leftOver);
	
	// return the random int
	return random%OutOf;
}

// returns a random node (equally weighted) connected to nodeIndex
inline
int Network::getNextRandomNode(int nodeIndex)
{
	// get the list of incident edges
	vector < int > theseEdges = m_nodeToEdgesMap[nodeIndex];
	
	// choose a pseudorandom edge index from theseEdges
	Edge randomEdge = m_edges[theseEdges[randomInt(theseEdges.size())]];
	
	// get the next node
	int otherNodeIndex=0;
	if(randomEdge.first == nodeIndex)
		otherNodeIndex = randomEdge.second;
	else if(randomEdge.second == nodeIndex)
		otherNodeIndex = randomEdge.first;
	
	// return the next node
	return otherNodeIndex;
}

// returns the number of steps a random walk takes to get from startNode to endNode
// returns the -1 if nodes are disconnected
int Network::takeRandomWalk(int startNode, int destNode)
{
	// check connectivity
	if(!sameConnectedComponent(startNode, destNode))
		return -1;
	
	// get the next random node until reaching the destNode
	int currentNode = startNode, numSteps = 0;
	do
	{
		currentNode = getNextRandomNode(currentNode);
		numSteps++;
	}while(currentNode != destNode);
	
	// you've reached your destination!
	return numSteps;
}

// returns the mean (over numWalks) first passage time from startNode to destNode
// returns -1 if the nodes are in disjoint subsets of the network
double Network::MFPT(int startNode, int destNode, int numWalks)
{
	// check connectivity
	if(!sameConnectedComponent(startNode, destNode))
		return -1;
	
	// take numWalks random walks and sum the total number of steps required
	int totalNumSteps = 0;
	for(int i = 0; i < numWalks; i++)
		totalNumSteps += takeRandomWalk(startNode, destNode);
	
	// return totalNumSteps / numWalks ( = mean FPT)
	return totalNumSteps/(double)numWalks;
}

////////////////////////////////////////////////////////////////////////////////
// Community Detection
////////////////////////////////////////////////////////////////////////////////

// closest friends - works with any square matrix of closeness measure (i.e. with coarse detection as well)
void Network::findCFs(const vector < vector < double > >& GENs, vector < int >& CFs)
{
	// num nodes
	int numNodes = GENs.size();
	
	// find closest friend of this node (smallest GEN from this node to friend node)
	for(int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
	{
		double smallestGEN;
		int closestFriend;
		if(nodeIndex == 0)
		{
			smallestGEN = 1000000000.0;
			closestFriend = 1;
		}
		else
		{
			smallestGEN = 1000000000.0;
			closestFriend = 0;
		}
		
		// find closest friend of this node (smallest GEN from this node to friend node)
		for(int i = 0; i < numNodes; i++)
			// don't check this node to itself
			if(i != nodeIndex)
				if(GENs[i][nodeIndex] < smallestGEN && GENs[i][nodeIndex] > 0)
				{
					smallestGEN = GENs[i][nodeIndex];
					closestFriend = i;
				}
		
		// return index of closest friend
		CFs.push_back(closestFriend);
	}
	
	// track friend definitions
	m_usedCFs = true;
	m_usedCUFs = false;
}

// closest unpopular friends: the closest friend with degree than or equal to the next closest friend
void Network::findCUFs(const vector < vector < double > >& GENs, vector < int >& CUFs)
{
	// num nodes
	int numNodes = GENs.size();
	
	// find closest friend of this node (smallest GEN from this node to friend node)
	for(int nodeIndex = 0; nodeIndex < numNodes; nodeIndex++)
	{
		// get the GENs from this node to the rest of the network
		vector < int > originalSort;
		vector < double > theseGENs;
		for(int i = 0; i < numNodes; i++)
			if(i != nodeIndex) // don't consider the GEN of this node to itself
			{
				double thisGEN = GENs[i][nodeIndex];
				if(thisGEN > 0) // don't consider GENs to disconnected nodes
				{
					originalSort.push_back(i);
					theseGENs.push_back(thisGEN);
				}
			}

		// sort theseGENs (selection)
		bool foundCUF = false;
		int numConnected = theseGENs.size();
		for(int i = 0; i < numConnected - 1; i++)
		{
			// find the smallest remaining GEN
			int smallestIndex = i;
			for(int j = i + 1; j < numConnected; j++)
				if(theseGENs[j] < theseGENs[smallestIndex] && theseGENs[j] >= 0.0)
					smallestIndex = j;
			
			// swap it to position i
			int swapOriginalIndex = originalSort[i];
			double swapGEN = theseGENs[i];
			originalSort[i] = originalSort[smallestIndex];
			theseGENs[i] = theseGENs[smallestIndex];
			originalSort[smallestIndex] = swapOriginalIndex;
			theseGENs[smallestIndex] = swapGEN;
			
			// checking the degree of each newly sorted node relative to the previously sorted node
			if(i > 0 && getDegree(originalSort[i]) <= getDegree(originalSort[i - 1]))
			{
				// you've found the CUF
				CUFs.push_back(originalSort[i]);
				foundCUF = true;
				break;
			}
			
			// if you're on the last i, check the i + 1 value against the i value
			if(i == numConnected - 2 && getDegree(originalSort[i + 1]) <= getDegree(originalSort[i]))
			{
				// you've found the CUF
				CUFs.push_back(originalSort[i + 1]);
				foundCUF = true;
				break;
			}
		}
		
		// if you didn't find a CUF because this node has only one connection, use the CF as the CUF
		if(!foundCUF)
			CUFs.push_back(originalSort[0]);
	}
	
	// track friend definitions
	m_usedCFs = false;
	m_usedCUFs = true;
}

// find the friend path of every node, given a single 'friend' for each node
void Network::buildFPs(const vector < int >& friends, vector < vector < int > >& friendPaths)
{	
	// get num 'nodes'
	int numNodes = friends.size();
	
	// build paths
	for(int i = 0; i < numNodes; i++)
	{
		// this path
		vector < int > thisPath;
		thisPath.push_back(i);
		int nextFriend = friends[i];
		
		// determine the path
		bool addToPath = true;
		while(addToPath)
		{
			// is nextFriend already in this path?
			for(int j = 0; j < thisPath.size(); j++)
				if(thisPath[j] == nextFriend)
				{
					addToPath = false;
					break;
				}
			
			// extend the path
			if(addToPath)
			{
				// add nextFriend to thisPath
				thisPath.push_back(nextFriend);
				
				// get the closest friend of the last node added to the path
				nextFriend = friends[nextFriend];
			}
		}
		
		// add thisPath to m_FPs
		friendPaths.push_back(thisPath);
	}
}

// reduce the set of CFPs (each associated with a node) into fine-grained communities
void Network::buildCommunities(const vector < vector < int > >& friendPaths, 
							   vector < vector < int > >& communities,
							   vector < int >& communityIDs)
{
	// number of nodes
	int numNodes = friendPaths.size();
	
	// set all community IDs to -1 in m_communities initially
	for(int i = 0; i < numNodes; i++)
		communityIDs.push_back(-1);
	
	// assume each friendPath[i] is a new community
	int numCommunities = 0;
	for(int i = 0; i < numNodes; i++)
	{
		// get this path size
		int thisPathLength = friendPaths[i].size();
		
		// assign new community IDs to nodes in this path
		// unless a pre-established community is found, if so: start over with the found ID
		int thisCommunityID = numCommunities;
		bool foundID = false;
		int foundStartIndex;
		for(int j = 0; j < thisPathLength; j++)
		{
			// check for exisiting community
			if(communityIDs[friendPaths[i][j]] >= 0)
			{
				thisCommunityID = communityIDs[friendPaths[i][j]];
				foundID = true;
				foundStartIndex = j;
				break;
			}
			// otherwise assume a new community
			else
				communityIDs[friendPaths[i][j]] = thisCommunityID;
		}
		
		// if an existing community in this path was found, assign exisint ID from the beginning of the path
		if(foundID)
			for(int j = 0; j < foundStartIndex; j++)
				communityIDs[friendPaths[i][j]] = thisCommunityID;
		// otherwise the new IDs are set and a new community kernel was found
		else
			numCommunities++;
	}
	
	// set up the vector of community vectors m_communities
	vector < int > emptyCommunity;
	int largestCommunityID = -1;
	for(int i = 0; i < numNodes; i++)
		if(communityIDs[i] > largestCommunityID)
		{
			communities.push_back(emptyCommunity);
			largestCommunityID++;
		}
	
	// add new node indices to each community
	for(int i = 0; i < numNodes; i++)
		communities[communityIDs[i]].push_back(i);
}

// initialize m_communities and m_communityIDs
// call buildCommunities on the first tier
void Network::buildFineCommunities()
{
	// initialize community structures
	vector < vector < int > > tempCommunities;
	vector < int > tempCommunityIDs;
	m_communities.push_back(tempCommunities);
	m_communityIDs.push_back(tempCommunityIDs);
	
	// build the first tier of communities
	buildCommunities(getFriendPaths(), getCommunities(0), getCommunityIDs(0));
}

// unite communities that may have been fractured during detection
void Network::repairFracturedCommunities(int tier)
{
	// get # edges
	int numEdges = getNumEdges();
	
	// do this until there are no more fractures found
	while(true)
	{
		// get initial # communities
		int numCommunities = getNumCommunities(tier);
		
		// find # edges between all pairs of communities
		vector < vector < int > > numEdgesBetweenCommunities;
		for(int g = 0; g < numCommunities; g++)
		{
			vector < int > theseNumEdgesBetweenCommunities;
			numEdgesBetweenCommunities.push_back(theseNumEdgesBetweenCommunities);
			for(int h = 0; h <= g; h++)
				numEdgesBetweenCommunities[g].push_back(0);
		}  //initialized the array of edges 
		
		// double count the number of edges between all pairs of communities
		for(int i = 0; i < numEdges; i++)  //looking at all edges
		{
			// get the communities across this edge
			Edge thisEdge(m_edges[i]);
			Edge theseCommunities(getCommunityOfNode(thisEdge.first, tier), getCommunityOfNode(thisEdge.second, tier));  //figure out the communities of the two nodes
			
			// swap the order of theseCommunities if first is smaller than second
			if(theseCommunities.first < theseCommunities.second)
			{
				int swap = theseCommunities.first;
				theseCommunities.first = theseCommunities.second;
				theseCommunities.second = swap;
			}
			//order the edges such that the communities are increasing 
			//for bookkeeping, doesn't affect method
			
			
			// count this weight in numEdgesBetweenCommunities
			numEdgesBetweenCommunities[theseCommunities.first][theseCommunities.second]+=m_edgeWeights[i];
			
		}
		
		// divide the double counted edges
		int totalEdges = 0;
		for(int i = 0; i < numCommunities; i++)
			for(int j = 0; j < numEdgesBetweenCommunities[i].size(); j++)
			{
				totalEdges += numEdgesBetweenCommunities[i][j];
			}
		//this is the total weight in the network
		
		// find the largest fractured score
		double largestFracturedScore = 0.0;
		int fracturedGroupA, fracturedGroupB; // A's index should always be < B's index
		bool foundFracture = false;
		for(int i = 0; i < numCommunities; i++)
			for(int j = 0; j < i; j++)
			{
				// compute this unique community pair's "fractured score" k_(g->h) / max(n_g, n_h)
				double thisFracturedScore = numEdgesBetweenCommunities[i][j] / 
					pow((double)max(getNodesInCommunity(i, tier).size(), getNodesInCommunity(j, tier).size()), 2);
				
				// determine whether this meets the fractured condition
				bool thisFracturedCondition = numEdgesBetweenCommunities[i][j] >
					min(numEdgesBetweenCommunities[i][i], numEdgesBetweenCommunities[j][j]);
				
				// if thisFracturedScore is the largest score thus far, record the fractured community pair
				if(thisFracturedScore > largestFracturedScore && thisFracturedCondition)
				{
					largestFracturedScore = thisFracturedScore;
					fracturedGroupA = j;
					fracturedGroupB = i;
					if(!foundFracture)
						foundFracture = true;
				}
			}
		
		// determine whether the location of the largest fractured score meets the fractured condition
		if(foundFracture)
		{
			// unite fractured community
			for(int i = fracturedGroupB + 1; i < numCommunities; i++)
			{
				// subtract 1 from the indices of all m_communityIDs larger than fracturedGroupB
				int thisCommSize = m_communities[tier][i].size();
				for(int j = 0; j < thisCommSize; j++)
					m_communityIDs[tier][m_communities[tier][i][j]]--;
			}
			int fracturedGroupB_Size = m_communities[tier][fracturedGroupB].size();
			for(int i = 0; i < fracturedGroupB_Size; i++)
				// change each community ID of the nodes in fracturedGroupB
				m_communityIDs[tier][m_communities[tier][fracturedGroupB][i]] = fracturedGroupA;
			
			// move all nodes in m_communities[fracturedGroupB] into m_communities[fracturedGroupA]
			m_communities[tier][fracturedGroupA].insert(m_communities[tier][fracturedGroupA].end(),
														m_communities[tier][fracturedGroupB].begin(),
														m_communities[tier][fracturedGroupB].end());
			
			// delete m_communities[fracturedGroupB]
			m_communities[tier].erase(m_communities[tier].begin() + fracturedGroupB);
		}
		
		// break if no fracture found 
		// i.e. the pair of communities with the largest fractured score doesn't meet the fractured condition
		else
			break;
	}		
}

// recursively build the coarse-grained community structure
// call this once on the GENs and the fine community structure at tier 0
void Network::buildCoarseCommunities(const vector < vector < double > >& GENs, int thisTier)
{
	// get # communities at the current tier
	int numCommunities = getNumCommunities(thisTier);
	
	// break recursion if the # of communities at thisTier is 1 (i.e. the entire network)
	if(numCommunities == 1)
		return;
	
	// otherwise compute the next tier of communities
	else
	{
		// initialize next tier
		vector < vector < int > > tempCommunities;
		vector < int > tempCommunityIDs;
		m_communities.push_back(tempCommunities);
		m_communityIDs.push_back(tempCommunityIDs);
		
		// community closeness at thisTier
		vector < vector < double > > communityCloseness;
		for(int i = 0; i < numCommunities; i++)
		{
			vector < double > tempCommunityCloseness;
			communityCloseness.push_back(tempCommunityCloseness);
			for(int j = 0; j < numCommunities; j++)
				communityCloseness[i].push_back(0.0);
		}
		
		// sum the closenesses felt from each community to each other community
		for(int i = 0; i < numCommunities; i++)
			for(int j = 0; j < numCommunities; j++)
				if(i != j)
				{
					// get the nodes and size of each community
					vector < int > nodesInI = getNodesInCommunity(i, thisTier);
					vector < int > nodesInJ = getNodesInCommunity(j, thisTier);
					int sizeOfI = nodesInI.size();
					int sizeOfJ = nodesInJ.size();
					
					// sum the closeness felt by each node in J -> each node in I
					for(int a = 0; a < sizeOfI; a++)
						for(int b = 0; b < sizeOfJ; b++)
							communityCloseness[i][j] += 1/GENs[nodesInI[a]][nodesInJ[b]];
					
					// invert these closeness sums and multiply each by sizeOfI * sizeOfJ
					communityCloseness[i][j] = sizeOfI * sizeOfJ / communityCloseness[i][j];
				}
		
		// find the closest friends
		vector < int > friends;
		findCFs(communityCloseness, friends);
		
		// build the friend paths
		vector < vector < int > > friendPaths;
		buildFPs(friends, friendPaths);
		
		// reduce the friend paths into the next tier of communities
		buildCommunities(friendPaths, m_communities[thisTier + 1], m_communityIDs[thisTier + 1]);
		
		// repair fractures in the next tier of communities
		repairFracturedCommunities(thisTier + 1);
		
		// recursively compute the next tier
		buildCoarseCommunities(GENs, thisTier + 1);
	}
}

// recursively aggregate the actual network nodes that belong to a community at any level of the hierarchy
vector < int > Network::getNodesInCommunity(int communityID, int tier)
{
	// recursive breaking condition -> if we're actually at the level above the nodes, return a vector of nodes
	if(tier == 0)
		return m_communities[0][communityID];
	
	// otherwise, call this function on all of the 'nodes' in this community
	else
	{
		vector < int > theseNodes;
		int numNodes = m_communities[tier][communityID].size();
		for(int i = 0; i < numNodes; i++)
		{
			vector < int > nextNodes = getNodesInCommunity(m_communities[tier][communityID][i] , tier - 1);
			theseNodes.insert(theseNodes.end(), nextNodes.begin(), nextNodes.end());
		}					  
		return theseNodes;
	}
}

// move up the community hierarchy from a node to the given tier, returning the community ID at that tier
int Network::getCommunityOfNode(int nodeIndex, int tier)
{
	// check that this tier exists
	if(m_communities.size() < tier + 1)
		return -1;
	
	// otherwise move up the community IDs at each tier
	else
	{
		int theCommunityID = m_communityIDs[0][nodeIndex];
		for(int i = 1; i <= tier; i++)
			theCommunityID = m_communityIDs[i][theCommunityID];
		return theCommunityID;
	}
}

#endif
