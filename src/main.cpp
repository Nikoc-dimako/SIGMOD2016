#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <unistd.h>
#include <algorithm>
#include <unordered_map>
#include <string.h>
#include <queue>
#include <pthread.h>
#include <set>
#include <limits.h>

#include <omp.h>

#include "threadpool11/threadpool11.hpp"

// #include <fstream>

using namespace std;
using namespace __gnu_cxx;

// One of the threads will be the master thread
#define NUMOFTHREADS 25

// Size of the node arrays
#define MAXN ((size_t)1<<26)

// Max value of the node number
#define NAMECHANGESIZE ((size_t)1<<28)		// INT_MAX

// Values for the proprocess to find quick queries
#define TOTALNODES 5
#define MAXDEPTH 5

/*
 * This is the factor to decide if a  graph is a star or not
 * In order to be a star the graph must have:
 * averageNeighbors*STARFACTOR < maxNeighbors
 * where maxNeighbors is the most neighbors a node has.
 */
#define STARFACTOR 5000

/*
 * If nodes/MERGEFACTOR > totalDeletes then merge
 */
#define MERGEFACTOR 10

// Type of a node
typedef unsigned int Node;

/* Visitor is used to save a state of a Node. Visited or not.
 * This value is not cleaned at every query so the value must increase.
 * The max number of queries before reset is the value of Visitor/2.
 */
typedef unsigned short int Visitor;

// This type of node is used for preprocessing
typedef struct GraphNode {
   vector<Node> nodes;
   int children;

   GraphNode(){
	   children = 0;
   }
} GraphNode;

// This type of node is used for star type Graphs where we need the children for predictability
typedef struct StarGNode {
	unsigned int start;
	unsigned int children	: 30,
				 addition	: 1,
				 deletion	: 1;

   StarGNode(){
	   start = 0;
	   addition = 0;
	   deletion = 0;
	   children = 0;
   }
} StarGNode;

// This type of node is at graphs where there is not a big variation in the number of neighbors
// We have better memory locality using this type of node
typedef struct FastGNode {
	unsigned int start	: 30,
				 addition	: 1,
				 deletion	: 1;
   FastGNode(){
	   start = 0;
	   addition=false;
	   deletion=false;
   }
} FastGNode;


// Struct used to save the query parameters, in order to enqueue the queries
typedef struct SP_Params{
	unsigned int a;
	unsigned int b;
	unsigned int versionCounter;
	unsigned int resultsCounter;
} SP_Params;


struct MyPair{
	Node b;
	int version;
};

struct Operation_info{
	char op;
	Node a,b;
};

// Variables for the Queues
StarGNode *StarForwardG, *StarBackwardG;
vector<Node> Nodes;

FastGNode *FastForwardG, *FastBackwardG;

GraphNode *ForwardGraph, *BackwardGraph;
Visitor **visited_global;
Node ***fQueue_global;
Node ***bQueue_global;
Visitor *visitedCounter_global;

threadpool11::Pool query_pool;


//Variables for the node's nameChange
Node *nameChange;
Node *nameReorder;
Node *nameChangeNewPositions;	// This array will hold the positions of the node's new position
Node *nodePositions;			// This array holds the value that is in this position
unsigned int nameChangeCounter = 1;
unsigned int initialNodes;


// Variables for the multiversion
vector<MyPair> *Forward_add,*Forward_del,*Backward_add,*Backward_del;
int results[100000];

int threadCounter = 0;
unordered_map <pthread_t, int> threadIds;
mutex mtx;

class Mycomparison
{
  bool reverse;
public:
  bool operator() (const SP_Params& q1, const SP_Params& q2) const
  {
	  int q1smaller = ForwardGraph[nameChange[q1.a]].children > BackwardGraph[nameChange[q1.b]].children ? nameChange[q1.a] : nameChange[q1.b];
	  int q2smaller = ForwardGraph[nameChange[q2.a]].children > BackwardGraph[nameChange[q2.b]].children ? nameChange[q2.a] : nameChange[q2.b];
	  return q1smaller > q2smaller;
  }
};


/*
 * This function is run at the beginning from every thread to link it with a number from 0 to NUMOFTHREADS
 * The thread will use that number to run it's jobs to the proper structs
 */
void initThread(){
	mtx.lock();
	unordered_map<pthread_t, int>::const_iterator um_it = threadIds.find(pthread_self());
	threadIds.emplace(pthread_self(), threadCounter);

	threadCounter++;
	mtx.unlock();

	sleep(1);
}


void shortest_path_star(Node tempa, Node tempb, int localVersion, int resultsCounter) {
	// Maybe the query is for a node that doesn't exist
	if(nameChange[tempa] == 0 || nameChange[tempb] == 0){
		if(tempa != tempb)
			results[resultsCounter] = -1;
		else
			results[resultsCounter] = 0;
		return;
	}

	Node a = nameChange[tempa];
	Node b = nameChange[tempb];

	if(a == b){
		results[resultsCounter] = 0;
		return;
	}

	// Get the proper structs to run the query on
	int thread_id = threadIds.find(pthread_self())->second;
	Visitor *visited = visited_global[thread_id];
	Node *fFrontQueue = fQueue_global[thread_id][0];
	Node *fNextQueue = fQueue_global[thread_id][1];
	Node *fNextQueue_iter;

	Node *bFrontQueue = bQueue_global[thread_id][0];
	Node *bNextQueue = bQueue_global[thread_id][1];
	Node *bNextQueue_iter;

	// If a thread completed a lot of queries maybe a reset is needed at the visited array
	if(visitedCounter_global[thread_id] == 65534){
		visitedCounter_global[thread_id] = 0;
		memset(visited, 0, MAXN);
	}

	/*
	 * In order not to reset the visited array at every query we have a number to represent 0
	 * We increase that number at every query, so we don't need to reset the array
	 */
	visitedCounter_global[thread_id] += 2;
	const unsigned int visitedCounter = visitedCounter_global[thread_id];

	const StarGNode *const FGraph = StarForwardG;
	const StarGNode *const BGraph = StarBackwardG;

	// Initializing Queues and Variables
	unsigned int fChildrenCount = FGraph[a+1].start - FGraph[a].start;
	unsigned int fCurrentNodes = 1;
	unsigned int bChildrenCount = BGraph[b+1].start - BGraph[b].start;
	unsigned int bCurrentNodes = 1;

	fFrontQueue[0] = a;
	bFrontQueue[0] = b;

	// Initializing Distances and visited Nodes
	unsigned int fGraphDistance = 0;
	unsigned int bGraphDistance = 0;
	visited[a] = visitedCounter;
	visited[b] = visitedCounter+1;

	// The queue saves the currently explored nodes
	// The children counters contain the number of children that the already explored nodes have.

	// Main Loop
	while(1){

		if(fChildrenCount <= bChildrenCount){	// Move forward, there are less children there
			fChildrenCount = 0;
			fNextQueue_iter = fNextQueue;
			fGraphDistance++;	// Going to the next distance

			// Reading all the children from the nodes in the queue
			for(int n=fCurrentNodes-1; n >= 0; n--){
				const unsigned int currentFather = fFrontQueue[n];
				fChildrenCount += FGraph[currentFather].children;

				// Check if this currentFather node has any deletions done at him
				if(FGraph[currentFather].deletion == 0){	// Just check the children nodes

					// Reading the children of the current node in the queue
					const unsigned int end = FGraph[currentFather+1].start;
					for(unsigned int j=FGraph[currentFather].start; j < end; j++){
						const Node currentChild = Nodes[j];
						if(visited[currentChild] < visitedCounter){	// Explored by the other side
							visited[currentChild] = visitedCounter;
							*fNextQueue_iter = currentChild;
							fNextQueue_iter++;
						}
						else{
							if(visited[currentChild] == visitedCounter+1){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
						}
					}
				}else{		// Deletion == 1
					// We have to check at every child if it was deleted
					const unsigned int end = FGraph[currentFather+1].start;
					for(unsigned int j=FGraph[currentFather].start; j < end; j++){
						const Node currentChild = Nodes[j];

						// Check if the node was deleted
						// If it was continue to the next child;

						std::vector<MyPair>:: iterator it;
						for(it=Forward_del[currentFather].begin(); it!=Forward_del[currentFather].end(); it++){
							if(it->b==currentChild){
								break;
							}
						}
						if(it!=Forward_del[currentFather].end() && it->version<localVersion){
							continue;
						}

						if(visited[currentChild] < visitedCounter){	// Explored by the other side
							visited[currentChild] = visitedCounter;
							*fNextQueue_iter = currentChild;
							fNextQueue_iter++;
						}
						else{
							if(visited[currentChild] == visitedCounter+1){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
						}
					}
				}


				// Check if there were any additions in the current father node
				if(FGraph[currentFather].addition == 1){
					if(thread_id==0){
						query_pool.postWork<void>([tempa, tempb, localVersion, resultsCounter] {  shortest_path_star(tempa, tempb, localVersion, resultsCounter);  });
						return;
					}
					// For every node added
					if(FGraph[currentFather].deletion == 0){
						for(std::vector<MyPair> ::iterator it=Forward_add[currentFather].begin(); it!=Forward_add[currentFather].end(); it++){

							// Just add the nodes to the queue
							const Node child=it->b;

							// The addition must have appeared before the query
							if(it->version > localVersion)
								continue;

							if(visited[child] < visitedCounter){	// Explored by the other side
								visited[child]=visitedCounter;
								*fNextQueue_iter = child;
								fNextQueue_iter++;
							}
							else{
								if(visited[child] == visitedCounter+1){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
								}
							}
						}
					}else{		// Deletion == 1
						// You have to check if it was deleted first
						for(std::vector<MyPair> ::iterator it=Forward_add[currentFather].begin(); it!=Forward_add[currentFather].end(); it++){

							const Node child=it->b;

							// The addition must have appeared before the query
							if(it->version > localVersion)
								continue;

							std::vector<MyPair>:: iterator itD;
							for(itD=Forward_del[currentFather].begin(); itD!=Forward_del[currentFather].end(); itD++){
								if(itD->b == child && itD->version > it->version){		// Delete must have happened after the addition
									break;
								}
							}
							if(itD!=Forward_del[currentFather].end() && itD->version<localVersion){
								continue;
							}

							if(visited[child] < visitedCounter){	// Explored by the other side
								visited[child]=visitedCounter;
								*fNextQueue_iter = child;
								fNextQueue_iter++;
							}
							else{
								if(visited[child] == visitedCounter+1){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
								}
							}
						}
					}
				}
			}

			if(fNextQueue == fNextQueue_iter){
				results[resultsCounter] = -1;
				return;
			}
			fCurrentNodes = (unsigned int)(fNextQueue_iter - fNextQueue);
			Node * temp = fNextQueue;
			fNextQueue = fFrontQueue;
			fFrontQueue = temp;

		}else{			// bChildrenCounter > fChildrenCounter
			bChildrenCount = 0;
			bNextQueue_iter = bNextQueue;
			bGraphDistance++;	// Going to the next distance

			// Reading all the children from the nodes in the queue
			for(int n=bCurrentNodes-1; n >= 0; n--){
				const unsigned int currentFather = bFrontQueue[n];
				bChildrenCount += BGraph[currentFather].children;

				// Check if this currentFather node has any deletions done at him
				if(BGraph[currentFather].deletion == 0){	// Just check the children nodes

					// Reading the children of the current node in the queue
					const unsigned int end = BGraph[currentFather+1].start;
					for(unsigned int j=BGraph[currentFather].start; j < end; j++){
						const Node currentChild = Nodes[j];
						if(visited[currentChild] < visitedCounter){	// Explored by the other side
							visited[currentChild] = visitedCounter + 1;
							*bNextQueue_iter = currentChild;
							bNextQueue_iter++;
							//iterator++;
						}
						else{
							if(visited[currentChild] == visitedCounter){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
						}
					}
				}else{		// Deletion == 1
					// We have to check at every child if it was deleted
					const unsigned int end = BGraph[currentFather+1].start;
					for(unsigned int j=BGraph[currentFather].start; j < end; j++){
						const Node currentChild = Nodes[j];

						// Check if the node was deleted
						// If it was continue to the next child;

						std::vector<MyPair>:: iterator it;
						for(it=Backward_del[currentFather].begin(); it!=Backward_del[currentFather].end(); it++){
							if(it->b==currentChild){
								break;
							}
						}
						if(it!=Backward_del[currentFather].end() && it->version<localVersion){
							continue;
						}

						if(visited[currentChild] < visitedCounter){	// Explored by the other side
							visited[currentChild] = visitedCounter + 1;
							*bNextQueue_iter = currentChild;
							bNextQueue_iter++;
						}
						else{
							if(visited[currentChild] == visitedCounter){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
						}
					}
				}


				// Check if there were any additions in the current father node
				if(BGraph[currentFather].addition == 1){
					if(thread_id==0){
						query_pool.postWork<void>([tempa, tempb, localVersion, resultsCounter] {  shortest_path_star(tempa, tempb, localVersion, resultsCounter);  });
						return;
					}
					// For every node added
					if(BGraph[currentFather].deletion == 0){
						for(std::vector<MyPair> ::iterator it=Backward_add[currentFather].begin(); it!=Backward_add[currentFather].end(); it++){

							// Just add the nodes to the queue
							const Node child=it->b;

							// The addition must have appeared before the query
							if(it->version > localVersion)
								continue;

							if(visited[child] < visitedCounter){	// Explored by the other side
								visited[child]=visitedCounter + 1;
								*bNextQueue_iter = child;
								bNextQueue_iter++;
							}
							else{
								if(visited[child] == visitedCounter){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
								}
							}
						}
					}else{		// Deletion == 1
						// You have to check if it was deleted first
						for(std::vector<MyPair> ::iterator it=Backward_add[currentFather].begin(); it!=Backward_add[currentFather].end(); it++){

							const Node child=it->b;

							// The addition must have appeared before the query
							if(it->version > localVersion)
								continue;

							std::vector<MyPair>:: iterator itD;
							for(itD=Backward_del[currentFather].begin(); itD!=Backward_del[currentFather].end(); itD++){
								if(itD->b == child && itD->version > it->version){		// Delete must have happened after the addition
									break;
								}
							}
							if(itD!=Backward_del[currentFather].end() && itD->version<localVersion){
								continue;
							}

							if(visited[child] < visitedCounter){	// Explored by the other side
								visited[child]=visitedCounter + 1;
								*bNextQueue_iter = child;
								bNextQueue_iter++;
							}
							else{
								if(visited[child] == visitedCounter){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
								}
							}
						}
					}
				}
			}

			if(bNextQueue_iter == bNextQueue){
				results[resultsCounter] = -1;
				return;
			}
			bCurrentNodes = (unsigned int)(bNextQueue_iter - bNextQueue);
			Node * temp = bNextQueue;
			bNextQueue = bFrontQueue;
			bFrontQueue = temp;
		}
	}
}


void shortest_path_fast(Node tempa, Node tempb, int localVersion, int resultsCounter) {
	// Maybe the query is for a node that doesn't exist
	if(nameChange[tempa] == 0 || nameChange[tempb] == 0){
		if(tempa != tempb)
			results[resultsCounter] = -1;
		else
			results[resultsCounter] = 0;
		return;
	}

	Node a = nameChange[tempa];
	Node b = nameChange[tempb];

	if(a == b){
		results[resultsCounter] = 0;
		return;
	}

	// Get the proper structs to run the query on
	int thread_id = threadIds.find(pthread_self())->second;
	Visitor *visited = visited_global[thread_id];
	Node *fFrontQueue = fQueue_global[thread_id][0];
	Node *fNextQueue = fQueue_global[thread_id][1];

	Node *bFrontQueue = bQueue_global[thread_id][0];
	Node *bNextQueue = bQueue_global[thread_id][1];

	// If a thread completed a lot of queries maybe a reset is needed at the visited array
	if(visitedCounter_global[thread_id] == 65534){
		visitedCounter_global[thread_id] = 0;
		memset(visited, 0, MAXN);
	}

	/*
	 * In order not to reset the visited array at every query we have a number to represent 0
	 * We increase that number at every query, so we don't need to reset the array
	 */
	visitedCounter_global[thread_id] += 2;
	const unsigned int visitedCounter = visitedCounter_global[thread_id];

	const FastGNode *const FGraph = FastForwardG;
	const FastGNode *const BGraph = FastBackwardG;

	// Initialize Queues and Variables
	unsigned int fCurrentNodes = 1;
	unsigned int bCurrentNodes = 1;

	fFrontQueue[0] = a;
	bFrontQueue[0] = b;

	// Initializing Distances and visited Nodes
	unsigned int fGraphDistance = 0;
	unsigned int bGraphDistance = 0;
	visited[a] = visitedCounter;
	visited[b] = visitedCounter+1;

	// The queue saves the currently explore1d nodes
	// The children counters contain the number of children that the already explored nodes have.

	// Main Loop
	while(1){
		unsigned int iterator = 0;

		if(fCurrentNodes <= bCurrentNodes){	// Move forward, there are less children there
			fGraphDistance++;	// Going to the next distance

			// Reading all the children from the nodes in the queue
			for(int n=fCurrentNodes-1; n >= 0; n--){
				const unsigned int currentFather = fFrontQueue[n];

				// Check if this currentFather node has any deletions done at him
				if(FGraph[currentFather].deletion == 0){	// Just check the children nodes

					// Reading the children of the current node in the queue
					const Node start = FGraph[currentFather].start;
					const unsigned int size = FGraph[currentFather+1].start - start;
					for(unsigned int j=0; j < size; j++){
						const Node currentChild = Nodes[start+j];
						if(visited[currentChild] < visitedCounter){	// Explored by the other side
							visited[currentChild] = visitedCounter;
							fNextQueue[iterator] = currentChild;
							iterator++;
						}
						else{
							if(visited[currentChild] == visitedCounter+1){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
						}
					}
				}else{		// Deletion == 1
					// We have to check at every child if it was deleted
					const Node start = FGraph[currentFather].start;
					const unsigned int size = FGraph[currentFather+1].start - start;
					for(unsigned int j=0; j < size; j++){
						const Node currentChild = Nodes[start+j];

						// Check if the node was deleted
						// If it was continue to the next child;

						std::vector<MyPair>:: iterator it;
						for(it=Forward_del[currentFather].begin(); it!=Forward_del[currentFather].end(); it++){
							if(it->b==currentChild){
								break;
							}
						}
						if(it!=Forward_del[currentFather].end() && it->version<localVersion){
							continue;
						}

						if(visited[currentChild] < visitedCounter){	// Explored by the other side
							visited[currentChild] = visitedCounter;
							fNextQueue[iterator] = currentChild;
							iterator++;
						}
						else{
							if(visited[currentChild] == visitedCounter+1){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
						}
					}
				}


				// Check if there were any additions in the current father node
				if(FGraph[currentFather].addition == 1){
					// For every node added
					if(FGraph[currentFather].deletion == 0){
						for(std::vector<MyPair> ::iterator it=Forward_add[currentFather].begin(); it!=Forward_add[currentFather].end(); it++){

							// Just add the nodes to the queue
							const Node child=it->b;

							// The addition must have appeared before the query
							if(it->version > localVersion)
								continue;

							if(visited[child] < visitedCounter){	// Explored by the other side
								visited[child]=visitedCounter;
								fNextQueue[iterator] = child;
								iterator++;
							}
							else{
								if(visited[child] == visitedCounter+1){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
								}
							}
						}
					}else{		// Deletion == 1
						// You have to check if it was deleted first
						for(std::vector<MyPair> ::iterator it=Forward_add[currentFather].begin(); it!=Forward_add[currentFather].end(); it++){

							const Node child=it->b;

							// The addition must have appeared before the query
							if(it->version > localVersion)
								continue;

							std::vector<MyPair>:: iterator itD;
							for(itD=Forward_del[currentFather].begin(); itD!=Forward_del[currentFather].end(); itD++){
								if(itD->b == child && itD->version > it->version){		// Delete must have happened after the addition
									break;
								}
							}
							if(itD!=Forward_del[currentFather].end() && itD->version<localVersion){
								continue;
							}

							if(visited[child] < visitedCounter){	// Explored by the other side
								visited[child]=visitedCounter;
								fNextQueue[iterator] = child;
								iterator++;
							}
							else{
								if(visited[child] == visitedCounter+1){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
								}
							}
						}
					}
				}
			}

			if(iterator == 0){
				results[resultsCounter] = -1;
				return;
			}
			fCurrentNodes = iterator;
			Node * temp = fNextQueue;
			fNextQueue = fFrontQueue;
			fFrontQueue = temp;

		}else{			// bChildrenCounter > fChildrenCounter
			bGraphDistance++;	// Going to the next distance

			// Reading all the children from the nodes in the queue
			for(int n=bCurrentNodes-1; n >= 0; n--){
				const unsigned int currentFather = bFrontQueue[n];

				// Check if this currentFather node has any deletions done at him
				if(BGraph[currentFather].deletion == 0){	// Just check the children nodes

					// Reading the children of the current node in the queue
					const Node start = BGraph[currentFather].start;
					const unsigned int size = BGraph[currentFather+1].start - start;
					for(unsigned int j=0; j < size; j++){
						const Node currentChild = Nodes[start+j];
						if(visited[currentChild] < visitedCounter){	// Explored by the other side
							visited[currentChild] = visitedCounter + 1;
							bNextQueue[iterator] = currentChild;
							iterator++;
						}
						else{
							if(visited[currentChild] == visitedCounter){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
						}
					}
				}else{		// Deletion == 1
					// We have to check at every child if it was deleted
					const Node start = BGraph[currentFather].start;
					const unsigned int size = BGraph[currentFather+1].start - start;
					for(unsigned int j=0; j < size; j++){
						const Node currentChild = Nodes[start+j];

						// Check if the node was deleted
						// If it was continue to the next child;

						std::vector<MyPair>:: iterator it;
						for(it=Backward_del[currentFather].begin(); it!=Backward_del[currentFather].end(); it++){
							if(it->b==currentChild){
								break;
							}
						}
						if(it!=Backward_del[currentFather].end() && it->version<localVersion){
							continue;
						}

						if(visited[currentChild] < visitedCounter){	// Explored by the other side
							visited[currentChild] = visitedCounter + 1;
							bNextQueue[iterator] = currentChild;
							iterator++;
						}
						else{
							if(visited[currentChild] == visitedCounter){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
						}
					}
				}


				// Check if there were any additions in the current father node
				if(BGraph[currentFather].addition == 1){
					// For every node added
					if(BGraph[currentFather].deletion == 0){
						for(std::vector<MyPair> ::iterator it=Backward_add[currentFather].begin(); it!=Backward_add[currentFather].end(); it++){

							// Just add the nodes to the queue
							const Node child=it->b;

							// The addition must have appeared before the query
							if(it->version > localVersion)
								continue;

							if(visited[child] < visitedCounter){	// Explored by the other side
								visited[child]=visitedCounter + 1;
								bNextQueue[iterator] = child;
								iterator++;
							}
							else{
								if(visited[child] == visitedCounter){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
								}
							}
						}
					}else{		// Deletion == 1
						// You have to check if it was deleted first
						for(std::vector<MyPair> ::iterator it=Backward_add[currentFather].begin(); it!=Backward_add[currentFather].end(); it++){

							const Node child=it->b;

							// The addition must have appeared before the query
							if(it->version > localVersion)
								continue;

							std::vector<MyPair>:: iterator itD;
							for(itD=Backward_del[currentFather].begin(); itD!=Backward_del[currentFather].end(); itD++){
								if(itD->b == child && itD->version > it->version){		// Delete must have happened after the addition
									break;
								}
							}
							if(itD!=Backward_del[currentFather].end() && itD->version<localVersion){
								continue;
							}

							if(visited[child] < visitedCounter){	// Explored by the other side
								visited[child]=visitedCounter + 1;
								bNextQueue[iterator] = child;
								iterator++;
							}
							else{
								if(visited[child] == visitedCounter){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
								}
							}
						}
					}
				}
			}

			if(iterator == 0){
				results[resultsCounter] = -1;
				return;
			}
			bCurrentNodes = iterator;
			Node * temp = bNextQueue;
			bNextQueue = bFrontQueue;
			bFrontQueue = temp;
		}
	}
}


void add_edge_final(Node a, Node b) {
  ForwardGraph[a].nodes.push_back(b);
  BackwardGraph[b].nodes.push_back(a);
}

int _delete_edge(Node a, vector<Node> &E) {
  int size = E.size(), newsize = 0;
  for(int i=0; i < size; i++){
    E[newsize] = E[i];
    if(E[i] != a)
		newsize++;
  }
  if(newsize == size)
	  return 0;
  E.resize(newsize);
  return 1;
}

void delete_edge_final(Node a,Node b) {
	if(!_delete_edge(b,ForwardGraph[a].nodes) || !_delete_edge(a,BackwardGraph[b].nodes) )
		return;
}

void add_edge_star(Node a, Node b, int versionCounter){
	MyPair fPair, bPair;
	fPair.version=versionCounter;
	fPair.b=b;

	bPair.version=versionCounter;
	bPair.b=a;

	Forward_add[a].push_back(fPair);
	StarForwardG[a].addition=true;

	Backward_add[b].push_back(bPair);
	StarBackwardG[b].addition=true;
}

void delete_edge_star(Node a, Node b, int versionCounter){
	MyPair fPair, bPair;
	fPair.version=versionCounter;
	fPair.b=b;

	bPair.version=versionCounter;
	bPair.b=a;

	Forward_del[a].push_back(fPair);
	StarForwardG[a].deletion=true;

	Backward_del[b].push_back(bPair);
	StarBackwardG[b].deletion=true;
}

void add_edge_fast(Node a, Node b, int versionCounter){
	MyPair fPair, bPair;
	fPair.version=versionCounter;
	fPair.b=b;

	bPair.version=versionCounter;
	bPair.b=a;

	Forward_add[a].push_back(fPair);
	FastForwardG[a].addition=true;

	Backward_add[b].push_back(bPair);
	FastBackwardG[b].addition=true;
}

void delete_edge_fast(Node a, Node b, int versionCounter){
	MyPair fPair, bPair;
	fPair.version=versionCounter;
	fPair.b=b;

	bPair.version=versionCounter;
	bPair.b=a;

	Forward_del[a].push_back(fPair);
	FastForwardG[a].deletion=true;

	Backward_del[b].push_back(bPair);
	FastBackwardG[b].deletion=true;
}

// Count the children of a node
void preprocess(){
	for(unsigned int i=1; i < nameChangeCounter; i++)
		for(unsigned int j=0; j < ForwardGraph[i].nodes.size(); j++)
			ForwardGraph[i].children += ForwardGraph[ForwardGraph[i].nodes[j]].nodes.size();

	for(unsigned int i=1; i < nameChangeCounter; i++)
		for(unsigned int j=0; j < BackwardGraph[i].nodes.size(); j++)
			BackwardGraph[i].children += BackwardGraph[BackwardGraph[i].nodes[j]].nodes.size();
}

// Order the list of vectors, the one with most children first
bool sortF(Node a,Node b) { return (ForwardGraph[a].nodes.size() > ForwardGraph[b].nodes.size()); }
bool sortB(Node a,Node b) { return (BackwardGraph[a].nodes.size() > BackwardGraph[b].nodes.size()); }

void preprocess2() {
  for(unsigned int a=1; a < nameChangeCounter; a++) {
    sort(ForwardGraph[a].nodes.begin(),ForwardGraph[a].nodes.end(),sortF);
    sort(BackwardGraph[a].nodes.begin(),BackwardGraph[a].nodes.end(),sortB);
  }
}

bool *evaluationF,*evaluationB;

bool neighbours5;
bool totalNodes5;
set<int> visited;
int sumNodes;

void dfs4(int x, int depth) {

    visited.insert(x);
    sumNodes++;

    if(depth > MAXDEPTH) {
        neighbours5 = true;
        return;
    }

    if(sumNodes > TOTALNODES) {
        totalNodes5 = true;
        return;
    }

    unsigned int size = ForwardGraph[x].nodes.size();
    for(unsigned int i=0;i<size;i++) {
        if(visited.find(ForwardGraph[x].nodes[i]) != visited.end()) continue;

        dfs4(ForwardGraph[x].nodes[i], depth+1);
        if(neighbours5 || totalNodes5) return;
    }
}

void preprocess4(){

	evaluationF[0] = false;

	for(unsigned int i=1; i<nameChangeCounter; i++){
        totalNodes5 = false;
		neighbours5 = false;
        sumNodes = 0;
        visited.clear();

        dfs4(i,0);

		if(neighbours5 || totalNodes5)
			evaluationF[i] = false;
		else
			evaluationF[i] = true;
	}
}
void dfs5(int x, int depth) {

    visited.insert(x);
    sumNodes++;

    if(depth > MAXDEPTH) {
        neighbours5 = true;
        return;
    }

    if(sumNodes > TOTALNODES) {
        totalNodes5 = true;
        return;
    }

    unsigned int size = BackwardGraph[x].nodes.size();
    for(unsigned int i=0;i<size;i++) {
        if(visited.find(BackwardGraph[x].nodes[i]) != visited.end()) continue;

        dfs5(BackwardGraph[x].nodes[i], depth+1);
        if(neighbours5 || totalNodes5) return;
    }
}

void preprocess5(){

	evaluationB[0] = false;

	for(unsigned int i=1; i<nameChangeCounter; i++){
        totalNodes5 = false;
		neighbours5 = false;
        sumNodes = 0;
        visited.clear();

        dfs5(i,0);

		if(neighbours5 || totalNodes5)
			evaluationB[i] = false;
		else
			evaluationB[i] = true;
	}
}


unsigned int getMaxNeighbor(){
	unsigned int max = 0;
	for(unsigned int i=1; i<nameChangeCounter; i++)
		max = (max>ForwardGraph[i].nodes.size())? max : ForwardGraph[i].nodes.size();

	for(unsigned int i=1; i<nameChangeCounter; i++)
		max = (max>BackwardGraph[i].nodes.size())? max : BackwardGraph[i].nodes.size();

	return max;
}

void merge_stargraph(){
	// Forward Graph
#pragma omp for schedule(static)
	for(unsigned int node=0; node<initialNodes; node++){
		// Delete the proper nodes
		const unsigned int end = StarForwardG[node+1].start;
		for(unsigned int child=StarForwardG[node].start; child < end; child++){
			const Node currentChild = Nodes[child];

			for(unsigned int del_node=0; del_node < Forward_del[node].size(); del_node++)
				if(Forward_del[node][del_node].b == currentChild){
					StarForwardG[node].children -= StarForwardG[child+1].start - StarForwardG[child].start;
					Nodes[child] = 0;
					break;
				}
		}

		// Remove additions that were deleted
		for(unsigned int add_node=0; add_node < Forward_add[node].size(); add_node++){
			Node current_add = Forward_add[node][add_node].b;
			int add_version = Forward_add[node][add_node].version;

			for(unsigned int del_node=0; del_node < Forward_del[node].size(); del_node++)
				if(Forward_del[node][del_node].b == current_add && Forward_del[node][del_node].version > add_version){
					Forward_add[node].erase(  Forward_add[node].begin()+add_node );
					add_node--;
					break;
				}
		}

		// Add additions to the graph if there is empty space
		bool allAdded = false;
		for(unsigned int add_node=0; add_node < Forward_add[node].size(); add_node++){
			Node current_add = Forward_add[node][add_node].b;
			allAdded = false;

			for(unsigned int child=StarForwardG[node].start; child < end; child++){
				if(Nodes[child] == 0){
					allAdded = true;
					StarForwardG[child].children += StarForwardG[child+1].start - StarForwardG[child].start;
					Nodes[child] = current_add;
					break;
				}
			}

			// If we get here, there are no more empty nodes to insert an element
			if(allAdded == false)
				break;
		}

		StarForwardG[node].deletion = false;
		if(allAdded)
			StarForwardG[node].addition = false;
	}


	// Backward Graph
#pragma omp for schedule(static)
	for(unsigned int node=0; node<initialNodes; node++){
		// Delete the proper nodes
		const unsigned int end = StarBackwardG[node+1].start;
		for(unsigned int child=StarBackwardG[node].start; child < end; child++){
			const Node currentChild = Nodes[child];

			for(unsigned int del_node=0; del_node < Backward_del[node].size(); del_node++)
				if(Backward_del[node][del_node].b == currentChild){
					StarBackwardG[node].children -= StarBackwardG[child+1].start - StarBackwardG[child].start;
					Nodes[child] = 0;
					break;
				}
		}

		// Remove additions that were deleted
		for(unsigned int add_node=0; add_node < Backward_add[node].size(); add_node++){
			Node current_add = Backward_add[node][add_node].b;
			int add_version = Backward_add[node][add_node].version;

			for(unsigned int del_node=0; del_node < Backward_del[node].size(); del_node++)
				if(Backward_del[node][del_node].b == current_add && Backward_del[node][del_node].version > add_version){
					Backward_add[node].erase(  Backward_add[node].begin()+add_node );
					add_node--;
					break;
				}
		}

		// Add additions to the graph if there is empty space
		bool allAdded = false;
		for(unsigned int add_node=0; add_node < Backward_add[node].size(); add_node++){
			Node current_add = Backward_add[node][add_node].b;
			allAdded = false;

			for(unsigned int child=StarBackwardG[node].start; child < end; child++){
				if(Nodes[child] == 0){
					allAdded = true;
					StarBackwardG[child].children += StarBackwardG[child+1].start - StarBackwardG[child].start;
					Nodes[child] = current_add;
					break;
				}
			}

			// If we get here, there are no more empty nodes to insert an element
			if(allAdded == false)
				break;
		}

		StarBackwardG[node].deletion = false;
		if(allAdded)
			StarBackwardG[node].addition = false;
	}
}

void merge_fastgraph(){
	// Forward Graph
#pragma omp for schedule(static)
	for(unsigned int node=0; node<initialNodes; node++){
		// Delete the proper nodes
		const unsigned int end = FastForwardG[node+1].start;
		for(unsigned int child=FastForwardG[node].start; child < end; child++){
			const Node currentChild = Nodes[child];

			for(unsigned int del_node=0; del_node < Forward_del[node].size(); del_node++)
				if(Forward_del[node][del_node].b == currentChild){
					Nodes[child] = 0;
					break;
				}
		}

		// Remove additions that were deleted
		for(unsigned int add_node=0; add_node < Forward_add[node].size(); add_node++){
			Node current_add = Forward_add[node][add_node].b;
			int add_version = Forward_add[node][add_node].version;

			for(unsigned int del_node=0; del_node < Forward_del[node].size(); del_node++)
				if(Forward_del[node][del_node].b == current_add && Forward_del[node][del_node].version > add_version){
					Forward_add[node].erase(  Forward_add[node].begin()+add_node );
					add_node--;
					break;
				}
		}

		// Add additions to the graph if there is empty space
		bool allAdded = false;
		for(unsigned int add_node=0; add_node < Forward_add[node].size(); add_node++){
			Node current_add = Forward_add[node][add_node].b;
			allAdded = false;

			for(unsigned int child=FastForwardG[node].start; child < end; child++){
				if(Nodes[child] == 0){
					allAdded = true;
					Nodes[child] = current_add;
					break;
				}
			}

			// If we get here, there are no more empty nodes to insert an element
			if(allAdded == false)
				break;
		}

		FastForwardG[node].deletion = false;
		if(allAdded)
			FastForwardG[node].addition = false;
	}


	// Backward Graph
#pragma omp for schedule(static)
	for(unsigned int node=0; node<initialNodes; node++){
		// Delete the proper nodes
		const unsigned int end = FastBackwardG[node+1].start;
		for(unsigned int child=FastBackwardG[node].start; child < end; child++){
			const Node currentChild = Nodes[child];

			for(unsigned int del_node=0; del_node < Backward_del[node].size(); del_node++)
				if(Backward_del[node][del_node].b == currentChild){
					Nodes[child] = 0;
					break;
				}
		}

		// Remove additions that were deleted
		for(unsigned int add_node=0; add_node < Backward_add[node].size(); add_node++){
			Node current_add = Backward_add[node][add_node].b;
			int add_version = Backward_add[node][add_node].version;

			for(unsigned int del_node=0; del_node < Backward_del[node].size(); del_node++)
				if(Backward_del[node][del_node].b == current_add && Backward_del[node][del_node].version > add_version){
					Backward_add[node].erase(  Backward_add[node].begin()+add_node );
					add_node--;
					break;
				}
		}

		// Add additions to the graph if there is empty space
		bool allAdded = false;
		for(unsigned int add_node=0; add_node < Backward_add[node].size(); add_node++){
			Node current_add = Backward_add[node][add_node].b;
			allAdded = false;

			for(unsigned int child=FastBackwardG[node].start; child < end; child++){
				if(Nodes[child] == 0){
					allAdded = true;
					Nodes[child] = current_add;
					break;
				}
			}

			// If we get here, there are no more empty nodes to insert an element
			if(allAdded == false)
				break;
		}

		FastBackwardG[node].deletion = false;
		if(allAdded)
			FastBackwardG[node].addition = false;
	}
}

int main() {
	ios_base::sync_with_stdio(false);
	cin.tie(nullptr);

	int versionCounter=0;
	Node a,b;
	char c;
	int resultsCounter = 0;

	// For profiling
	//ifstream input;
	//input.open("/home/thanasis/Desktop/test-harness/inputDenser.txt");
	//ofstream output;
	//output.open("/home/thanasis/Desktop/test-harness/output.txt");

	ForwardGraph = (GraphNode*)calloc(MAXN, sizeof(GraphNode));
	BackwardGraph = (GraphNode*)calloc(MAXN, sizeof(GraphNode));
	nameChange = (Node*)calloc(NAMECHANGESIZE, sizeof(Node));
	Forward_del = (vector<MyPair>*)calloc(MAXN, sizeof(vector<MyPair>));
	Forward_add = (vector<MyPair>*)calloc(MAXN, sizeof(vector<MyPair>));
	Backward_del = (vector<MyPair>*)calloc(MAXN, sizeof(vector<MyPair>));
	Backward_add = (vector<MyPair>*)calloc(MAXN, sizeof(vector<MyPair>));

	// Creating structures for each thread
	fQueue_global = (Node***)calloc(NUMOFTHREADS, sizeof(Node**));
	bQueue_global = (Node***)calloc(NUMOFTHREADS, sizeof(Node**));
	visited_global = (Visitor**)calloc(NUMOFTHREADS, sizeof(Visitor*));
	visitedCounter_global = (Visitor*)calloc(NUMOFTHREADS, sizeof(Visitor));

	for(int i=0; i<NUMOFTHREADS; i++){
		fQueue_global[i] = (Node**)calloc(2, sizeof(Node*));
		bQueue_global[i] = (Node**)calloc(2, sizeof(Node*));
		visited_global[i] = (Visitor*)calloc(MAXN, sizeof(Visitor));
	}

	for(int i=0; i<NUMOFTHREADS; i++)
		for(int j=0; j<2; j++){
			fQueue_global[i][j] = (Node*)calloc(MAXN, sizeof(Node));
			bQueue_global[i][j] = (Node*)calloc(MAXN, sizeof(Node));
		}

	unsigned int counter = 0;
	while(cin >> a >> b){
		if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
		if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
		add_edge_final(nameChange[a],nameChange[b]);
		//if(counter == 30622564){		// Input normal: 3232855, Input Denser: 30622564,	for profiling
		//	break;
		//}
		counter++;
	}
	initialNodes = nameChangeCounter;

	for(unsigned int i=0; i<initialNodes+100; i++){
		Forward_del[i].reserve(15);
		Forward_add[i].reserve(15);
		Backward_del[i].reserve(15);
		Backward_add[i].reserve(15);
	}

	Nodes.reserve(counter*2);

	evaluationF = (bool *)calloc(nameChangeCounter,sizeof(bool));
	evaluationB = (bool *)calloc(nameChangeCounter,sizeof(bool));

	preprocess();			// Sorting the biggest Nodes to the front of the lists
	preprocess2();			// Count children


	unsigned int maxNeighbors = getMaxNeighbor();
	double averageNeighbors = (counter)/nameChangeCounter;

	bool star = false;
	if(averageNeighbors*STARFACTOR<maxNeighbors)
		star=true;


	if(star){
		// Create the new Graph
		StarForwardG = (StarGNode*)calloc(MAXN, sizeof(StarGNode));
		StarBackwardG = (StarGNode*)calloc(MAXN, sizeof(StarGNode));

		unsigned int nodes_ptr = 0;
		for(unsigned int cur_n=1; cur_n<nameChangeCounter; cur_n++){
			unsigned int size = ForwardGraph[cur_n].nodes.size();
			StarForwardG[cur_n].start = nodes_ptr;
			StarForwardG[cur_n].children = ForwardGraph[cur_n].children;;
			for(unsigned int child=0; child < size; child++, nodes_ptr++)
				Nodes.push_back(ForwardGraph[cur_n].nodes[child]);
		}
		StarForwardG[nameChangeCounter].start = nodes_ptr;

		StarBackwardG[0].start = nodes_ptr;
		for(unsigned int cur_n=1; cur_n<nameChangeCounter; cur_n++){
			unsigned int size = BackwardGraph[cur_n].nodes.size();
			StarBackwardG[cur_n].start = nodes_ptr;
			StarBackwardG[cur_n].children = BackwardGraph[cur_n].children;;
			for(unsigned int child=0; child < size; child++, nodes_ptr++)
				Nodes.push_back(BackwardGraph[cur_n].nodes[child]);
		}
		StarBackwardG[nameChangeCounter].start = nodes_ptr;

		preprocess4();
		preprocess5();

		// free(ForwardGraph);
		// free(BackwardGraph);

		// Creating Threads
		query_pool.setWorkerCount(NUMOFTHREADS-1);
		initThread();
		for(int i=0; i<NUMOFTHREADS-1; i++)
			query_pool.postWork<void>(initThread);
		query_pool.waitAll();

		vector<SP_Params> main_queue;
		main_queue.reserve(20000);

		vector<SP_Params> thread_queue;
		thread_queue.reserve(20000);

		SP_Params params;

		unsigned int totalDeletes = 0;

		cout << "R" << endl << flush;


		cin.clear();
		cin >> c;

		if(averageNeighbors*2*STARFACTOR<maxNeighbors){		// StarGraph with short queries run by single thread
			while(cin >> c) {
				if(c == 'F') {
					// Maybe merge the graph
					if(nameChangeCounter/MERGEFACTOR < totalDeletes){
						totalDeletes = 0;
						merge_stargraph();
					}

					// Give the queries to the threads
					sort(thread_queue.begin(),thread_queue.end(),Mycomparison());
					for(vector<SP_Params>::iterator it=thread_queue.begin(); it!=thread_queue.end(); it++)
						query_pool.postWork<void>([it] {  shortest_path_star(it->a, it->b, it->versionCounter, it->resultsCounter);  });

					// Run your shortest paths
					sort(main_queue.begin(),main_queue.end(),Mycomparison());
					for(vector<SP_Params>::iterator it=main_queue.begin(); it!=main_queue.end(); it++)
						shortest_path_star(it->a, it->b, it->versionCounter, it->resultsCounter);

					thread_queue.clear();
	 				main_queue.clear();
					query_pool.waitAll();
					// Print the results
					for(int i=0; i<resultsCounter; i++)
						cout << results[i] << endl;

					cout << flush;
					cin.clear();

					resultsCounter = 0;

					continue;
				}

				cin >> a >> b;

				if(c == 'Q'){
					if((nameChange[a]<initialNodes && evaluationF[nameChange[a]] && !StarForwardG[nameChange[a]].addition)
						|| (nameChange[b]<initialNodes && evaluationB[nameChange[b]] && !StarBackwardG[nameChange[b]].addition)){
						params.a = a;
						params.b = b;
						params.versionCounter = versionCounter;
						params.resultsCounter = resultsCounter;
						main_queue.push_back(params);
					}else{
						params.a = a;
						params.b = b;
						params.versionCounter = versionCounter;
						params.resultsCounter = resultsCounter;
						thread_queue.push_back(params);
					}

					versionCounter++;
					resultsCounter++;

				}else if(c == 'A'){
					if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
					if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

					add_edge_star(nameChange[a], nameChange[b], versionCounter);

					versionCounter++;

				}else{
					totalDeletes++;
					if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
					if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

					delete_edge_star(nameChange[a], nameChange[b], versionCounter);

					versionCounter++;
				}
			}
		}else{					// StarGraph without short queries run by single thread
			while(cin >> c){
				if(c == 'F') {
					// Give the queries to the threads
					sort(thread_queue.begin(),thread_queue.end(),Mycomparison());
					for(vector<SP_Params>::iterator it=thread_queue.begin(); it!=thread_queue.end(); it++)
						query_pool.postWork<void>([it] {  shortest_path_star(it->a, it->b, it->versionCounter, it->resultsCounter);  });
					thread_queue.clear();

					query_pool.waitAll();
					// Print the results
					for(int i=0; i<resultsCounter; i++)
						cout << results[i] << endl;

					cout << flush;
					cin.clear();

					resultsCounter = 0;

					continue;
				}

				cin >> a >> b;

				if(c == 'Q'){
					params.a = a;
					params.b = b;
					params.versionCounter = versionCounter;
					params.resultsCounter = resultsCounter;
					thread_queue.push_back(params);

					versionCounter++;
					resultsCounter++;

				}else if(c == 'A'){
					if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
					if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

					add_edge_star(nameChange[a], nameChange[b], versionCounter);

					versionCounter++;

				}else{
					if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
					if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

					delete_edge_star(nameChange[a], nameChange[b], versionCounter);

					versionCounter++;
				}
			}
		}
	}else{				// FastGraph without children for predictability
		// Create the new Graph
		FastForwardG = (FastGNode*)calloc(MAXN, sizeof(FastGNode));
		FastBackwardG = (FastGNode*)calloc(MAXN, sizeof(FastGNode));

		unsigned int nodes_ptr = 0;
		for(unsigned int cur_n=1; cur_n<nameChangeCounter; cur_n++){
			unsigned int size = ForwardGraph[cur_n].nodes.size();
			FastForwardG[cur_n].start = nodes_ptr;
			for(unsigned int child=0; child < size; child++, nodes_ptr++)
				Nodes.push_back(ForwardGraph[cur_n].nodes[child]);
		}
		FastForwardG[nameChangeCounter].start = nodes_ptr;

		for(unsigned int cur_n=1; cur_n<nameChangeCounter; cur_n++){
			unsigned int size = BackwardGraph[cur_n].nodes.size();
			FastBackwardG[cur_n].start = nodes_ptr;
			for(unsigned int child=0; child < size; child++, nodes_ptr++)
				Nodes.push_back(BackwardGraph[cur_n].nodes[child]);
		}
		FastBackwardG[nameChangeCounter].start = nodes_ptr;

		preprocess4();
		preprocess5();

		// free(ForwardGraph);
		// free(BackwardGraph);

		// Creating Threads
		query_pool.setWorkerCount(NUMOFTHREADS-1);

		initThread();
		for(int i=0; i<NUMOFTHREADS-1; i++)
			query_pool.postWork<void>(initThread);
		query_pool.waitAll();

		// Main Thread Queue
		std::vector<SP_Params> main_queue;
		main_queue.reserve(20000);

		std::vector<SP_Params> thread_queue;
		thread_queue.reserve(20000);


		SP_Params params;

		cout << "R" << endl << flush;

		cin.clear();
		cin >> c;

		while(cin >> c) {
			if(c == 'F') {
				// Give the queries to the threads
				sort(thread_queue.begin(),thread_queue.end(),Mycomparison());
				for(vector<SP_Params>::iterator it=thread_queue.begin(); it!=thread_queue.end(); it++)
					query_pool.postWork<void>([it] {  shortest_path_fast(it->a, it->b, it->versionCounter, it->resultsCounter);  });
				thread_queue.clear();

				query_pool.waitAll();
				// Print the results
				for(int i=0; i<resultsCounter; i++)
					cout << results[i] << endl;

				cout << flush;
				cin.clear();

				resultsCounter = 0;

				continue;
			}

			cin >> a >> b;

			if(c == 'Q'){
				params.a = a;
				params.b = b;
				params.versionCounter = versionCounter;
				params.resultsCounter = resultsCounter;
				thread_queue.push_back(params);

				versionCounter++;
				resultsCounter++;

			}else if(c == 'A'){
				if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
				if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

				add_edge_fast(nameChange[a], nameChange[b], versionCounter);

				versionCounter++;

			}else{
				if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
				if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

				delete_edge_fast(nameChange[a], nameChange[b], versionCounter);

				versionCounter++;
			}
		}
	}

	return 0;
}
