#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <unistd.h>
#include <algorithm>
#include <unordered_map>
#include <string.h>

#include "threadpool11/threadpool11.hpp"

//#include <fstream>


using namespace std;
using namespace __gnu_cxx;

#define NUMOFTHREADS 24

#define MAXN ((size_t)1<<26)
#define NAMECHANGESIZE ((size_t)1<<28)

#define TOTALNODES 50

typedef unsigned int Node;
typedef unsigned short int Visitor;

typedef struct GraphNode {
   vector<Node> nodes;
   int children;

   GraphNode(){
	   children = 0;
   }
} GraphNode;

typedef struct StarGNode {
	unsigned int start;
	unsigned int addition	: 1,
				 deletion	: 1,
				 children	: 30;
	/*unsigned int children;
	bool addition;
	bool deletion;*/

   StarGNode(){
	   addition = 0;
	   deletion = 0;
	   children = 0;
	   start = 0;
   }
} StarGNode;

typedef struct FastGNode {
	unsigned int addition	: 1,
				 deletion	: 1,
				 start	: 30;
   FastGNode(){
	   addition=false;
	   deletion=false;
	   start = 0;
   }
} FastGNode;



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


//Variables for the node's nameChange
Node *nameChange;
Node *nameReorder;
Node *nameChangeNewPositions;	// This array will hold the positions of the node's new position
Node *nodePositions;			// This array holds the value that is in this position
unsigned int nameChangeCounter = 1;


// Variables for the multiversion
vector<MyPair> *Forward_add,*Forward_del,*Backward_add,*Backward_del;
int results[100000];

int threadCounter = 0;
unordered_map <pthread_t, int> threadIds;
mutex mtx;

void initThread(){
	mtx.lock();
	unordered_map<pthread_t, int>::const_iterator um_it = threadIds.find(pthread_self());
	//cerr << "Thread id " <<  pthread_self() << " -> "<< threadCounter << endl;
	threadIds.emplace(pthread_self(), threadCounter);

	 /*cpu_set_t cpuset;
	 CPU_ZERO(&cpuset);
	 CPU_SET(threadCounter, &cpuset);

	 int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
	 if (rc != 0)
		 std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";*/


	threadCounter++;
	mtx.unlock();

	sleep(1);
}

void shortest_path_star(Node tempa, Node tempb, int localVersion, int resultsCounter) {
	int thread_id = threadIds.find(pthread_self())->second;
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


	Visitor *visited = visited_global[thread_id];
	Node *fFrontQueue = fQueue_global[thread_id][0];
	Node *fNextQueue = fQueue_global[thread_id][1];
	Node *fNextQueue_iter;

	Node *bFrontQueue = bQueue_global[thread_id][0];
	Node *bNextQueue = bQueue_global[thread_id][1];
	Node *bNextQueue_iter;


	if(visitedCounter_global[thread_id] == 65534)
		memset(visited, 0, MAXN);

	visitedCounter_global[thread_id] += 2;
	const unsigned int visitedCounter = visitedCounter_global[thread_id];

	const StarGNode *const FGraph = StarForwardG;
	const StarGNode *const BGraph = StarBackwardG;

	// Initialize Queues and Variables
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

	// The queue saves the currently explore1d nodes
	// The children counters contain the number of children that the already explored nodes have.

	// Main Loop
	while(1){
		//unsigned int iterator = 0;

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
							//iterator++;
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
							//iterator++;
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
								*fNextQueue_iter = child;
								fNextQueue_iter++;
								//iterator++;
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
								//iterator++;
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
							//iterator++;
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
								*bNextQueue_iter = child;
								bNextQueue_iter++;
								//iterator++;
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
								//iterator++;
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
	int thread_id = threadIds.find(pthread_self())->second;
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


	Visitor *visited = visited_global[thread_id];
	Node *fFrontQueue = fQueue_global[thread_id][0];
	Node *fNextQueue = fQueue_global[thread_id][1];

	Node *bFrontQueue = bQueue_global[thread_id][0];
	Node *bNextQueue = bQueue_global[thread_id][1];

	if(visitedCounter_global[thread_id] == 65534)
		memset(visited, 0, MAXN);

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

void preprocess4(){
	int sumNeighbours;
	bool neighbours5;
	evaluationF[0] = false;

	for(unsigned int i=1; i<nameChangeCounter; i++){
		neighbours5 = false;
		sumNeighbours = ForwardGraph[i].children;

		vector<Node> &start1 = ForwardGraph[i].nodes;
		unsigned int size = ForwardGraph[i].nodes.size();
		for(unsigned int j=0; j < size; j++){
			if(ForwardGraph[start1[j]].children){
				neighbours5 = true;
				break;
			}
		}

		if(neighbours5 || sumNeighbours > TOTALNODES)
			evaluationF[i] = false;
		else
			evaluationF[i] = true;
	}
}

void preprocess5(){
	int sumNeighbours;
	bool neighbours5;
	evaluationB[0] = false;

	for(unsigned int i=1; i<nameChangeCounter; i++){
		neighbours5 = false;
		sumNeighbours = BackwardGraph[i].children;

		vector<Node> &start1 = BackwardGraph[i].nodes;
		unsigned int size = BackwardGraph[i].nodes.size();
		for(unsigned int j=0; j < size; j++){
			if(BackwardGraph[start1[j]].children){
				neighbours5 = true;
				break;
			}
		}

		if(neighbours5 || sumNeighbours > TOTALNODES)
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

/*void merge_star(){
	for(unsigned int i=0; i<nameChangeCounter; i++){
		if(StarForwardG[i].deletion && StarForwardG[i].addition){
			unsigned int


		}
	}
}*/

int main() {
	ios_base::sync_with_stdio(false);
	cin.tie(nullptr);

	int versionCounter=0;
	Node a,b;
	char c;
	int resultsCounter = 0;

	//For debugging
	//ifstream input;
	//input.open("/home/thanasis/Desktop/test-harness/input.txt");
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

	unsigned int counter = 0;			// For profiling
	while(cin >> a >> b){
		if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
		if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
		add_edge_final(nameChange[a],nameChange[b]);
		//if(counter == 3232855){		// For profiling
		//	break;
		//}
		counter++;
	}
	unsigned int initialNameCounter = nameChangeCounter;

	Nodes.reserve(counter*2);

	evaluationF = (bool *)calloc(nameChangeCounter,sizeof(bool));
	evaluationB = (bool *)calloc(nameChangeCounter,sizeof(bool));

	preprocess();			// Sorting the biggest Nodes to the front of the lists
	preprocess2();			// Count children


	unsigned int maxNeighbors = getMaxNeighbor();
	//unsigned int averageNeighbors = (unsigned int)(counter)/nameChangeCounter;

	bool star = false;
	if((nameChangeCounter/1000)*maxNeighbors > counter)
		star = true;

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

		free(ForwardGraph);
		free(BackwardGraph);

		// Creating Threads
		threadpool11::Pool query_pool;
		query_pool.setWorkerCount(NUMOFTHREADS-1);

		for(int i=0; i<NUMOFTHREADS-1; i++)
			query_pool.postWork<void>(initThread);
		initThread();
		query_pool.waitAll();

		// Main Thread Queue
		vector<SP_Params> main_queue;
		main_queue.reserve(20000);

		vector<SP_Params> thread_queue;
		thread_queue.reserve(20000);


		SP_Params params;

		cout << "R" << endl << flush;

		//cerr << "Starting Queries!" << endl;

		sleep(2);

		cin.clear();
		cin >> c;

		while(cin >> c) {
			if(c == 'F') {
				//cerr << "Starting Flush!" << endl;

				// Give the queries to the threads
				for(vector<SP_Params>::iterator it=thread_queue.begin(); it!=thread_queue.end(); it++)
					query_pool.postWork<void>([it] {  shortest_path_star(it->a, it->b, it->versionCounter, it->resultsCounter);  });
				thread_queue.clear();

				// Run your shortest paths
				for(vector<SP_Params>::iterator it=main_queue.begin(); it!=main_queue.end(); it++)
					shortest_path_star(it->a, it->b, it->versionCounter, it->resultsCounter);
				main_queue.clear();

				query_pool.waitAll();
				// Print the results
				for(int i=0; i<resultsCounter; i++)
					cout << results[i] << endl;

				cout << flush;
				cin.clear();

				resultsCounter = 0;

				//cerr << "End of Flush!" << endl;

				continue;
			}

			cin >> a >> b;

			if(c == 'Q'){

				//cerr << "Query!" << endl;

				if((nameChange[a]<initialNameCounter && evaluationF[nameChange[a]] && !StarForwardG[nameChange[a]].addition)
					|| (nameChange[b]<initialNameCounter && evaluationB[nameChange[b]] && !StarBackwardG[nameChange[b]].addition)){
					params.a = a;
					params.b = b;
					params.versionCounter = versionCounter;
					params.resultsCounter = resultsCounter;
					main_queue.push_back(params);
					//shortest_path(a, b, versionCounter, resultsCounter);
				}else{
					params.a = a;
					params.b = b;
					params.versionCounter = versionCounter;
					params.resultsCounter = resultsCounter;
					thread_queue.push_back(params);
				}
					//query_pool.postWork<void>([a, b, versionCounter, resultsCounter] {  shortest_path_star(a, b, versionCounter, resultsCounter);  });

				versionCounter++;
				resultsCounter++;

			}else if(c == 'A'){

				//cerr << "Addition!" << endl;

				if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
				if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

				add_edge_star(nameChange[a], nameChange[b], versionCounter);

				versionCounter++;

			}else{

				//cerr << "Deletion!" << endl;

				if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
				if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

				delete_edge_star(nameChange[a], nameChange[b], versionCounter);

				versionCounter++;
			}
		}
	}else{
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

		free(ForwardGraph);
		free(BackwardGraph);

		// Creating Threads
		threadpool11::Pool query_pool;
		query_pool.setWorkerCount(NUMOFTHREADS-1);

		for(int i=0; i<NUMOFTHREADS-1; i++)
			query_pool.postWork<void>(initThread);
		initThread();
		query_pool.waitAll();

		// Main Thread Queue
		vector<SP_Params> main_queue;
		main_queue.reserve(20000);

		vector<SP_Params> thread_queue;
		thread_queue.reserve(20000);


		SP_Params params;

		cout << "R" << endl << flush;

		//cerr << "Starting ... " << endl;

		sleep(2);

		cin.clear();
		cin >> c;

		while(cin >> c) {
			if(c == 'F') {
				//cerr << "Flushing ... " << endl;

				// Give the queries to the threads
				for(vector<SP_Params>::iterator it=thread_queue.begin(); it!=thread_queue.end(); it++)
					query_pool.postWork<void>([it] {  shortest_path_fast(it->a, it->b, it->versionCounter, it->resultsCounter);  });
				thread_queue.clear();

				// Run your shortest paths
				for(vector<SP_Params>::iterator it=main_queue.begin(); it!=main_queue.end(); it++)
					shortest_path_fast(it->a, it->b, it->versionCounter, it->resultsCounter);
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
				//cerr << "Query ... " << endl;

				if((nameChange[a]<initialNameCounter && evaluationF[nameChange[a]] && !FastForwardG[nameChange[a]].addition)
					|| (nameChange[b]<initialNameCounter && evaluationB[nameChange[b]] && !FastBackwardG[nameChange[b]].addition)){
					params.a = a;
					params.b = b;
					params.versionCounter = versionCounter;
					params.resultsCounter = resultsCounter;
					main_queue.push_back(params);
					//shortest_path_fast(a, b, versionCounter, resultsCounter);
				}else{
					params.a = a;
					params.b = b;
					params.versionCounter = versionCounter;
					params.resultsCounter = resultsCounter;
					thread_queue.push_back(params);
				}
					//query_pool.postWork<void>([a, b, versionCounter, resultsCounter] {  shortest_path_fast(a, b, versionCounter, resultsCounter);  });

				versionCounter++;
				resultsCounter++;

			}else if(c == 'A'){
				//cerr << "add ... " << endl;

				if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
				if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

				add_edge_fast(nameChange[a], nameChange[b], versionCounter);

				versionCounter++;

			}else{

				//cerr << "deletion ... " << endl;

				if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
				if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

				delete_edge_fast(nameChange[a], nameChange[b], versionCounter);

				versionCounter++;
			}
		}
	}

	return 0;
}

