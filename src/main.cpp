#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <unistd.h>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <queue>
#include <pthread.h>

#include "threadpool11/threadpool11.hpp"

//#include <fstream>

using namespace std;
using namespace __gnu_cxx;

#define NUMOFTHREADS 3

#define MAXN ((size_t)1<<24)
#define NAMECHANGESIZE ((size_t)1<<26)

#define TOTALNODES 5

/*
abs(ForwardGraph[q1.a].nodes.size() - BackwardGraph[q1.b].nodes.size()) <=  abs(ForwardGraph[q2.a].nodes.size() - BackwardGraph[q2.b].nodes.size());
ForwardGraph[q1.a].nodes.size() > ForwardGraph[q2.a].nodes.size()  &&  BackwardGraph[q1.b].nodes.size() > BackwardGraph[q2.b].nodes.size();
*/
typedef unsigned int Node;
typedef unsigned short int Visitor;

typedef struct GraphNode {
   vector<Node> nodes;
   int children;

   GraphNode(){
	   children = 0;
   }
} GraphNode;

// Try without children
typedef struct GNode {
	unsigned int start;
	unsigned int numOfNodes;
	unsigned int children;
	bool addition;
	bool deletion;

   GNode(){
	   addition=false;
	   deletion=false;
	   children = 0;
	   numOfNodes = 0;
	   start = 0;
   }
} GNode;


struct MyPair{
	Node b;
	int version;
};

struct Query
{
	Node a,b;
	int resultsCounter;
	int localVersion;
};


struct Operation_info{
	char op;
	Node a,b;
};

// Variables for the Queues
GNode *ForwardG, *BackwardG;
vector<Node> Nodes;

Visitor **visited_global;
Node ***fQueue_global;
Node ***bQueue_global;
Node *fQueue_single;
Node *bQueue_single;
Node *visited_single;
unsigned int visitedCounter_single;
Visitor *visitedCounter_global;

// Node **visited_global;
// Node **fQueue_global;
// Node **bQueue_global;


//Variables for the node's nameChange
Node *nameChange;
Node *nameReorder;
Node *nameChangeNewPositions;	// This array will hold the positions of the node's new position
Node *nodePositions;			// This array holds the value that is in this position
unsigned int nameChangeCounter = 1;


// Variables for the multiversion
vector<MyPair> *Forward_add,*Forward_del,*Backward_add,*Backward_del;
int results[100000];
vector<int> rslts;

class Mycomparison
{
  bool reverse;
public:
  bool operator() (const Query& q1, const Query& q2) const
  {
	  int q1smaller = ForwardG[nameChange[q1.a]].children > BackwardG[nameChange[q1.b]].children ? nameChange[q1.a] : nameChange[q1.b];
	  int q2smaller = ForwardG[nameChange[q2.a]].children > BackwardG[nameChange[q2.b]].children ? nameChange[q2.a] : nameChange[q2.b];

   	return q1smaller > q2smaller;
  }
};

priority_queue<Query, std::vector<Query>, Mycomparison> queries;

int threadCounter = 0;
//unordered_map <pthread_t, int> threadIds;
//mutex mtx;
//mutex mtx_q,mtx_work;
//condition_variable noemtpy,nowork;
pthread_mutex_t mtx_q,mtx_work;
int counter_threads=0;
pthread_cond_t noemtpy,nowork;
bool work_done = false;

/*void initThread(){
	mtx.lock();
	unordered_map<pthread_t, int>::const_iterator um_it = threadIds.find(pthread_self());
	//cerr << "Thread id " <<  pthread_self() << " -> "<< threadCounter << endl;
	threadIds.emplace(pthread_self(), threadCounter);

	 cpu_set_t cpuset;
	 CPU_ZERO(&cpuset);
	 CPU_SET(threadCounter, &cpuset);

	 int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
	 if (rc != 0)
		 std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";


	threadCounter++;
	mtx.unlock();

	sleep(1);
}*/


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////single thread



void shortest_path(Node tempa, Node tempb, int localVersion, int resultsCounter) {
	if(nameChange[tempa] == 0 || nameChange[tempb] == 0){
		if(tempa != tempb)
			results[resultsCounter] = -1;
		else
			results[resultsCounter] = 0;
		return;
	}

	Node a=nameChange[tempa];
	Node b=nameChange[tempb];

	if(a == b){
		results[resultsCounter] = 0;
		return;
	}

	Node *visited = visited_single;
	Node *fQueue = fQueue_single;
	Node *bQueue = bQueue_single;

	visitedCounter_single += 2;
	unsigned int &visitedCounter = visitedCounter_single;


	// Initialize Queues and Variables
	unsigned int fChildrenCount = ForwardG[a].numOfNodes;
	unsigned int fCurrentNodes = 1;
	unsigned int bChildrenCount = BackwardG[b].numOfNodes;
	unsigned int bCurrentNodes = 1;

	unsigned int fQueuePointer = 0;
	unsigned int bQueuePointer = 0;
	fQueue[0+fQueuePointer*MAXN] = a;
	bQueue[0+bQueuePointer*MAXN] = b;

	// Initializing Distances and visited Nodes
	unsigned int fGraphDistance = 0;
	unsigned int bGraphDistance = 0;
	visited[a] = visitedCounter;
	visited[b] = visitedCounter+1;

	unsigned int iterator;
	// The queue saves the currently explored nodes
	// The children counters contain the number of children that the already explored nodes have.

	// Main Loop
	while(1){
		iterator = 0;

		if(fChildrenCount <= bChildrenCount){	// Move forward, there are less children there
			fChildrenCount = 0;
			fGraphDistance++;	// Going to the next distance

			// Reading all the children from the nodes in the queue
			for(unsigned int i=0; i < fCurrentNodes; i++){
				Node currentFather = fQueue[i+fQueuePointer*MAXN];
				fChildrenCount += ForwardG[currentFather].children;

				// Check if this currentFather node has any deletions done at him
				if(ForwardG[currentFather].deletion == 0){	// Just check the children nodes

					// Reading the children of the current node in the queue
					Node start = ForwardG[currentFather].start;
					for(unsigned int j=0; j < ForwardG[currentFather].numOfNodes; j++){
						Node currentChild = Nodes[start+j];
						if(visited[currentChild] >= visitedCounter){	// Explored by the other side
							if(visited[currentChild] == visitedCounter+1){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
							continue;
						}
						visited[currentChild] = visitedCounter;
						fQueue[iterator+(fQueuePointer^1)*MAXN] = currentChild;
						iterator++;
					}
				}else{		// Deletion == 1
					// TODO


					// We have to check at every child if it was deleted
					Node start = ForwardG[currentFather].start;
					for(unsigned int j=0; j < ForwardG[currentFather].numOfNodes; j++){
						Node currentChild = Nodes[start+j];

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

						if(visited[currentChild] >= visitedCounter){	// Explored by the other side
							if(visited[currentChild] == visitedCounter+1){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
							continue;
						}
						visited[currentChild] = visitedCounter;
						fQueue[iterator+(fQueuePointer^1)*MAXN] = currentChild;
						iterator++;
					}
				}

				// Check if there were any additions in the current father node
				if(ForwardG[currentFather].addition == 1){
					// For every node added
					if(ForwardG[currentFather].deletion == 0){
						for(std::vector<MyPair> ::iterator it=Forward_add[currentFather].begin(); it!=Forward_add[currentFather].end(); it++){

							// Just add the nodes to the queue
							Node child=it->b;

							// The addition must have appeared before the query
							if(it->version > localVersion)
								continue;

							if(visited[child] >= visitedCounter){	// Explored by the other side
								if(visited[child] == visitedCounter+1){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
								}
								continue;
							}
							visited[child]=visitedCounter;
							fQueue[iterator+(fQueuePointer^1)*MAXN]=child;
							iterator++;
						}
					}else{		// Deletion == 1
						// You have to check if it was deleted first
						for(std::vector<MyPair> ::iterator it=Forward_add[currentFather].begin(); it!=Forward_add[currentFather].end(); it++){

							Node child=it->b;

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

							if(visited[child] >= visitedCounter){	// Explored by the other side
								if(visited[child] == visitedCounter+1){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
								}
								continue;
							}
							visited[child]=visitedCounter;
							fQueue[iterator+(fQueuePointer^1)*MAXN]=child;
							iterator++;
						}
					}
				}
			}

			if(iterator == 0){
				results[resultsCounter] = -1;
				return;
			}
			fCurrentNodes = iterator;
			fQueuePointer = fQueuePointer^1;

		}else{			// bChildrenCounter > fChildrenCounter
			bChildrenCount = 0;
			bGraphDistance++;	// Going to the next distance

			// Reading all the children from the nodes in the queue
			for(unsigned int i=0; i < bCurrentNodes; i++){
				Node currentFather = bQueue[i+bQueuePointer*MAXN];
				bChildrenCount += BackwardG[currentFather].children;

				// Check if this currentFather node has any deletions done at him
				if(BackwardG[currentFather].deletion == 0){	// Just check the children nodes

					// Reading the children of the current node in the queue
					Node start = BackwardG[currentFather].start;
					for(unsigned int j=0; j < BackwardG[currentFather].numOfNodes; j++){
						Node currentChild = Nodes[start+j];
						if(visited[currentChild] >= visitedCounter){	// Explored by the other side
							if(visited[currentChild] == visitedCounter){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
							continue;
						}
						visited[currentChild] = visitedCounter+1;
						bQueue[iterator+(bQueuePointer^1)*MAXN] = currentChild;
						iterator++;
					}
				}else{		// Deletion == 1
					// TODO

					// We have to check at every child if it was deleted
					Node start = BackwardG[currentFather].start;
					for(unsigned int j=0; j < BackwardG[currentFather].numOfNodes; j++){
						Node currentChild = Nodes[start+j];

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

						if(visited[currentChild] >= visitedCounter){	// Explored by the other side
							if(visited[currentChild] == visitedCounter){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
							continue;
						}
						visited[currentChild] = visitedCounter+1;
						bQueue[iterator+(bQueuePointer^1)*MAXN] = currentChild;
						iterator++;
					}
				}

				// Check if there were any additions in the current father node
				if(BackwardG[currentFather].addition == 1){
					// For every node added
					for(std::vector<MyPair> ::iterator it=Backward_add[currentFather].begin(); it!=Backward_add[currentFather].end(); it++){
						if(BackwardG[currentFather].deletion == 0){

							// Just add the nodes to the queue
							Node child=it->b;

							// The addition must have appeared before the query
							if(it->version > localVersion)
								continue;

							if(visited[child] >= visitedCounter){	// Explored by the other side
								if(visited[child] == visitedCounter){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
								}
								continue;
							}
							visited[child]=visitedCounter+1;
							bQueue[iterator+(bQueuePointer^1)*MAXN]=child;
							iterator++;

						}else{		// Deletion == 1
							// You have to check if it was deleted first
							Node child=it->b;

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

							if(visited[child] >= visitedCounter){	// Explored by the other side
								if(visited[child] == visitedCounter){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
								}
								continue;
							}
							visited[child]=visitedCounter+1;
							bQueue[iterator+(bQueuePointer^1)*MAXN]=child;
							iterator++;
						}
					}
				}
			}

			if(iterator == 0){
				results[resultsCounter] = -1;
				return;
			}
			bCurrentNodes = iterator;
			bQueuePointer = bQueuePointer^1;
		}
	}
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// multithread






void shortest_path(int thread_id) {

	 cpu_set_t cpuset;
	 CPU_ZERO(&cpuset);
	 CPU_SET(thread_id, &cpuset);

	 int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
	 if (rc != 0)
		 std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";

	//int thread_id = threadIds.find(pthread_self())->second;
	Visitor *visited = visited_global[thread_id];
	Node **fQueue = fQueue_global[thread_id];
	Node **bQueue = bQueue_global[thread_id];
	Query q;
	Node tempa,tempb,a,b;
	int resultsCounter,localVersion;
	int result;
	unsigned int fChildrenCount;
	unsigned int fCurrentNodes;
	unsigned int bChildrenCount;
	unsigned int bCurrentNodes;
	unsigned int fQueuePointer;
	unsigned int bQueuePointer;
	unsigned int fGraphDistance;
	unsigned int bGraphDistance;

	//cout << "3ekinaw id " << thread_id << endl;

	while(!work_done){

		pthread_mutex_lock (&mtx_q);
		while(queries.size()==0){
			pthread_cond_wait (&noemtpy, &mtx_q);
		}
		q=queries.top()	;
		queries.pop();
		pthread_mutex_unlock (&mtx_q);

		pthread_mutex_lock (&mtx_work);
		counter_threads++;
		pthread_mutex_unlock (&mtx_work);

		tempa=q.a; tempb=q.b; resultsCounter=q.resultsCounter; localVersion=q.localVersion;

		visitedCounter_global[thread_id] += 2;
		const unsigned int visitedCounter = visitedCounter_global[thread_id];

		if(nameChange[tempa] == 0 || nameChange[tempb] == 0){
			if(tempa != tempb)
				result= -1;
			else
				result= 0;
			goto LOOP;
		}

		a = nameChange[tempa];
	 	b = nameChange[tempb];
		// Initialize Queues and Variables
		fChildrenCount = ForwardG[a].numOfNodes;
		fCurrentNodes = 1;
		bChildrenCount = BackwardG[b].numOfNodes;
		bCurrentNodes = 1;

		fQueuePointer = 0;
		bQueuePointer = 0;
		fQueue[fQueuePointer][0] = a;
		bQueue[bQueuePointer][0] = b;

		// Initializing Distances and visited Nodes
		fGraphDistance = 0;
		bGraphDistance = 0;
		visited[a] = visitedCounter;
		visited[b] = visitedCounter+1;

		// The queue saves the currently explore1d nodes
		// The children counters contain the number of children that the already explored nodes have.



		if(a == b){
			result = 0;
			goto LOOP;
		}


		// Main Loop
		while(1){
			unsigned int iterator = 0;

			if(fChildrenCount <= bChildrenCount){	// Move forward, there are less children there
				fChildrenCount = 0;
				fGraphDistance++;	// Going to the next distance

				// Reading all the children from the nodes in the queue
				for(unsigned int n=0; n < fCurrentNodes; n++){
					const Node currentFather = fQueue[fQueuePointer][n];
					fChildrenCount += ForwardG[currentFather].children;

					// Check if this currentFather node has any deletions done at him
					if(ForwardG[currentFather].deletion == 0){	// Just check the children nodes

						// Reading the children of the current node in the queue
						Node start = ForwardG[currentFather].start;
						for(unsigned int i=0; i < ForwardG[currentFather].numOfNodes; i++){
							Node currentChild = Nodes[start+i];
							if(visited[currentChild] >= visitedCounter){	// Explored by the other side
								if(visited[currentChild] == visitedCounter+1){ 	// Found the minimum distance!
									result = fGraphDistance + bGraphDistance;
									goto LOOP;
								}
								continue;
							}
							visited[currentChild] = visitedCounter;
							fQueue[fQueuePointer^1][iterator] = currentChild;
							iterator++;
						}
					}else{		// Deletion == 1
						// We have to check at every child if it was deleted
						Node start = ForwardG[currentFather].start;
						for(unsigned int i=0; i < ForwardG[currentFather].numOfNodes; i++){
							Node currentChild = Nodes[start+i];

							// Check if the node was deleted
							// If it was continue to the next child;

							const unsigned int end = Forward_del[currentFather].size();
							unsigned int j = 0;
							for(; j < end; j++){
								if(Forward_del[currentFather][j].b == currentChild){
									break;
								}
							}
							if(j != end && Forward_del[currentFather][j].version < localVersion){
								continue;
							}

							if(visited[currentChild] >= visitedCounter){	// Explored by the other side
								if(visited[currentChild] == visitedCounter+1){ 	// Found the minimum distance!
									result = fGraphDistance + bGraphDistance;
									goto LOOP;
								}
								continue;
							}
							visited[currentChild] = visitedCounter;
							fQueue[fQueuePointer^1][iterator] = currentChild;
							iterator++;
						}
					}

					// Check if there were any additions in the current father node
					if(ForwardG[currentFather].addition == 1){
						// For every node added
						if(ForwardG[currentFather].deletion == 0){

							const unsigned int end = Forward_add[currentFather].size();
							for(unsigned int i = 0; i < end; i++){
								// Just add the nodes to the queue
								const Node child = Forward_add[currentFather][i].b;

								// The addition must have appeared before the query
								if(Forward_add[currentFather][i].version > localVersion)
									continue;

								if(visited[child] >= visitedCounter){	// Explored by the other side
									if(visited[child] == visitedCounter+1){ 	// Found the minimum distance!
										result = fGraphDistance + bGraphDistance;
										goto LOOP;
									}
									continue;
								}
								visited[child]=visitedCounter;
								fQueue[fQueuePointer^1][iterator] = child;
								iterator++;
							}
						}else{		// Deletion == 1
							// You have to check if it was deleted first
							const unsigned int end = Forward_add[currentFather].size();
							for(unsigned int i = 0; i < end; i++){
								const Node child = Forward_add[currentFather][i].b;

								// The addition must have appeared before the query
								if(Forward_add[currentFather][i].version > localVersion)
									continue;

								const unsigned int end = Forward_del[currentFather].size();
								unsigned int j = 0;
								for(; j < end; j++){
									if(Forward_del[currentFather][j].b == child && Forward_del[currentFather][j].version > Forward_add[currentFather][i].version){		// Delete must have happened after the addition
										break;
									}
								}
								if(j != end && Forward_del[currentFather][j].version < localVersion){
									continue;
								}

								if(visited[child] >= visitedCounter){	// Explored by the other side
									if(visited[child] == visitedCounter+1){ 	// Found the minimum distance!
										result = fGraphDistance + bGraphDistance;
										goto LOOP;
									}
									continue;
								}
								visited[child]=visitedCounter;
								fQueue[fQueuePointer^1][iterator] = child;
								iterator++;
							}
						}
					}
				}

				if(iterator == 0){
					result= -1;
					goto LOOP;
				}
				fCurrentNodes = iterator;
				fQueuePointer = fQueuePointer^1;

			}else{			// bChildrenCounter > fChildrenCounter
				bChildrenCount = 0;
				bGraphDistance++;	// Going to the next distance

				// Reading all the children from the nodes in the queue
				for(unsigned int n=0; n < bCurrentNodes; n++){
					const Node currentFather = bQueue[bQueuePointer][n];
					bChildrenCount += BackwardG[currentFather].children;

					// Check if this currentFather node has any deletions done at him
					if(BackwardG[currentFather].deletion == 0){	// Just check the children nodes

						// Reading the children of the current node in the queue
						Node start = BackwardG[currentFather].start;
						for(unsigned int i=0; i < BackwardG[currentFather].numOfNodes; i++){
							Node currentChild = Nodes[start+i];
							if(visited[currentChild] >= visitedCounter){	// Explored by the other side
								if(visited[currentChild] == visitedCounter){ 	// Found the minimum distance!
									result = fGraphDistance + bGraphDistance;
									goto LOOP;
								}
								continue;
							}
							visited[currentChild] = visitedCounter+1;
							bQueue[bQueuePointer^1][iterator] = currentChild;
							iterator++;
						}
					}else{		// Deletion == 1
						// We have to check at every child if it was deleted
						Node start = BackwardG[currentFather].start;
						for(unsigned int i=0; i < BackwardG[currentFather].numOfNodes; i++){
							Node currentChild = Nodes[start+i];

							// Check if the node was deleted
							// If it was continue to the next child;

							const unsigned int end = Backward_del[currentFather].size();
							unsigned int j = 0;
							for(; j < end; j++){
								if(Backward_del[currentFather][j].b == currentChild){
									break;
								}
							}
							if(j != end && Backward_del[currentFather][j].version < localVersion){
								continue;
							}

							if(visited[currentChild] >= visitedCounter){	// Explored by the other side
								if(visited[currentChild] == visitedCounter){ 	// Found the minimum distance!
									result= fGraphDistance + bGraphDistance;
									goto LOOP;
								}
								continue;
							}
							visited[currentChild] = visitedCounter+1;
							bQueue[bQueuePointer^1][iterator] = currentChild;
							iterator++;
						}
					}

					// Check if there were any additions in the current father node
					if(BackwardG[currentFather].addition == 1){
						// For every node added
						if(BackwardG[currentFather].deletion == 0){

							const unsigned int end = Backward_add[currentFather].size();
							for(unsigned int i = 0; i < end; i++){
								// Just add the nodes to the queue
								const Node child = Backward_add[currentFather][i].b;

								// The addition must have appeared before the query
								if(Backward_add[currentFather][i].version > localVersion)
									continue;

								if(visited[child] >= visitedCounter){	// Explored by the other side
									if(visited[child] == visitedCounter){ 	// Found the minimum distance!
										result = fGraphDistance + bGraphDistance;
										goto LOOP;
									}
									continue;
								}
								visited[child]=visitedCounter+1;
								bQueue[bQueuePointer^1][iterator] = child;
								iterator++;
							}
						}else{		// Deletion == 1
							// You have to check if it was deleted first
							const unsigned int end = Backward_add[currentFather].size();
							for(unsigned int i = 0; i < end; i++){
								const Node child = Backward_add[currentFather][i].b;

								// The addition must have appeared before the query
								if(Backward_add[currentFather][i].version > localVersion)
									continue;

								const unsigned int end = Backward_del[currentFather].size();
								unsigned int j = 0;
								for(; j < end; j++){
									if(Backward_del[currentFather][j].b == child && Backward_del[currentFather][j].version > Backward_add[currentFather][i].version){		// Delete must have happened after the addition
										break;
									}
								}
								if(j != end && Backward_del[currentFather][j].version < localVersion){
									continue;
								}

								if(visited[child] >= visitedCounter){	// Explored by the other side
									if(visited[child] == visitedCounter){ 	// Found the minimum distance!
										result = fGraphDistance + bGraphDistance;
										goto LOOP;
									}
									continue;
								}
								visited[child]=visitedCounter+1;
								bQueue[bQueuePointer^1][iterator] = child;
								iterator++;
							}
						}
					}
				}

				if(iterator == 0){
					result = -1;
					goto LOOP;
				}
				bCurrentNodes = iterator;
				bQueuePointer = bQueuePointer^1;
			}
		}
		LOOP:
		results[resultsCounter]=result;
		pthread_mutex_lock (&mtx_work);
		counter_threads--;
		rslts.push_back(result);
		if(counter_threads==0){
			pthread_cond_signal (&nowork);
		}
		pthread_mutex_unlock (&mtx_work);
		//cout << "Result from " << thread_id << " is " << results[resultsCounter] << endl;
	}
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

void add_edge(Node a, Node b, int versionCounter){
	MyPair fPair, bPair;
	fPair.version=versionCounter;
	fPair.b=b;

	bPair.version=versionCounter;
	bPair.b=a;

	Forward_add[a].push_back(fPair);
	ForwardG[a].addition=true;

	Backward_add[b].push_back(bPair);
	BackwardG[b].addition=true;
}

void delete_edge(Node a, Node b, int versionCounter){
	MyPair fPair, bPair;
	fPair.version=versionCounter;
	fPair.b=b;

	bPair.version=versionCounter;
	bPair.b=a;

	Forward_del[a].push_back(fPair);
	ForwardG[a].deletion=true;

	Backward_del[b].push_back(bPair);
	BackwardG[b].deletion=true;
}

GraphNode *ForwardGraph, *BackwardGraph;

void add_edge_final(Node a, Node b) {
  ForwardGraph[a].nodes.push_back(b);
  BackwardGraph[b].nodes.push_back(a);
}

void delete_edge_final(Node a,Node b) {
	if(!_delete_edge(b,ForwardGraph[a].nodes) || !_delete_edge(a,BackwardGraph[b].nodes) )
		return;
}

// Order the list of vectors, the one with most children first
bool sortF(Node a,Node b) { return (ForwardGraph[a].nodes.size() > ForwardGraph[b].nodes.size()); }
bool sortB(Node a,Node b) { return (BackwardGraph[a].nodes.size() > BackwardGraph[b].nodes.size()); }

void preprocess() {
  for(unsigned int a=1; a < nameChangeCounter; a++) {
    sort(ForwardGraph[a].nodes.begin(),ForwardGraph[a].nodes.end(),sortF);
    sort(BackwardGraph[a].nodes.begin(),BackwardGraph[a].nodes.end(),sortB);
  }
}

// Count the children of a node
void preprocess2(){
	for(unsigned int i=1; i < nameChangeCounter; i++){
		unsigned int start = ForwardG[i].start;
		for(unsigned int j=0; j < ForwardG[i].numOfNodes; j++)
			ForwardG[i].children += ForwardG[Nodes[start+j]].numOfNodes;
	}

	for(unsigned int i=1; i < nameChangeCounter; i++){
		unsigned int start = BackwardG[i].start;
		for(unsigned int j=0; j < BackwardG[i].numOfNodes; j++)
			BackwardG[i].children += BackwardG[Nodes[start+j]].numOfNodes;
	}
}


bool *evaluationF,*evaluationB;

void preprocess4(){
	int sumNeighbours;
	bool level5;
	for(unsigned int i=1; i<nameChangeCounter; i++){
		level5=false;
		sumNeighbours = ForwardG[i].children;

		unsigned int start1 = ForwardG[i].start;
		for(unsigned int j=0; j < ForwardG[i].numOfNodes; j++){

			sumNeighbours += ForwardG[Nodes[start1+j]].children;
			unsigned int start2 = ForwardG[Nodes[start1+j]].start;
			for(unsigned int k=0; k < ForwardG[Nodes[start2]].numOfNodes; k++){
				level5=true;
				sumNeighbours += ForwardG[Nodes[start2+k]].children;
				if(sumNeighbours > TOTALNODES)
					break;
			}
		}

		if(!level5 && sumNeighbours<TOTALNODES)
			evaluationF[i] = true;
		else
			evaluationF[i] = false;
	}
}

void preprocess5(){
	int sumNeighbours;
	bool level5;
	for(unsigned int i=1; i<nameChangeCounter; i++){
		level5=false;
		sumNeighbours = BackwardG[i].children;

		unsigned int start1 = BackwardG[i].start;
		for(unsigned int j=0; j < BackwardG[i].numOfNodes; j++){

			sumNeighbours += BackwardG[Nodes[start1+j]].children;
			unsigned int start2 = BackwardG[Nodes[start1+j]].start;
			for(unsigned int k=0; k < BackwardG[Nodes[start2]].numOfNodes; k++){
				level5=true;
				sumNeighbours += BackwardG[Nodes[start2+k]].children;
				if(sumNeighbours > TOTALNODES)
					break;
			}
		}

		if(!level5 && sumNeighbours<TOTALNODES)
			evaluationF[i] = true;
		else
			evaluationF[i] = false;
	}
}



int main() {
	ios_base::sync_with_stdio(false);
	cin.tie(nullptr);

	//For debugging
	//ifstream input;
	//input.open("/home/thanasis/Desktop/test-harness/input.txt");
	//ofstream output;
	//output.open("/home/thanasis/Desktop/test-harness/output.txt");

	int versionCounter=0;
	Node a,b;
	char c;
	int resultsCounter = 0;
	priority_queue<Query, std::vector<Query>, Mycomparison> singleQ;
	unsigned int multiQuery=0;
	unsigned int initialNameCounter =0;

	// vector<Operation_info> operations(10000);
	// Operation_info info;

	ForwardG = (GNode*)calloc(MAXN, sizeof(GNode));
	BackwardG = (GNode*)calloc(MAXN, sizeof(GNode));
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
	fQueue_single=(Node *)calloc(2*MAXN, sizeof(Node));
	bQueue_single=(Node *)calloc(2*MAXN, sizeof(Node));
	visited_single=(Node*)calloc(MAXN, sizeof(Node));


	for(int i=0; i<NUMOFTHREADS; i++)
		for(int j=0; j<2; j++){
			fQueue_global[i][j] = (Node*)calloc(MAXN, sizeof(Node));
			bQueue_global[i][j] = (Node*)calloc(MAXN, sizeof(Node));
		}

	threadpool11::Pool query_pool;
	query_pool.setWorkerCount(NUMOFTHREADS);

	for(int i=0; i<NUMOFTHREADS; i++){
		query_pool.postWork<void>([i] {shortest_path(i);});
	}
	rslts.reserve(100000);

	int counter = 0;			// For profiling
	while(cin >> a >> b){
		if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
		if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
		add_edge_final(nameChange[a],nameChange[b]);
		//if(counter == 3232855){		// For profiling
		//	break;
		//}
		counter++;
	}

	Nodes.reserve(counter*2);

	preprocess();			// Sorting the biggest Nodes to the front of the lists



	// Create the new Graph
	unsigned int nodes_ptr = 0;
	for(unsigned int cur_n=1; cur_n<nameChangeCounter; cur_n++){
		unsigned int size = ForwardGraph[cur_n].nodes.size();
		ForwardG[cur_n].start = nodes_ptr;
		ForwardG[cur_n].numOfNodes = size;
		for(unsigned int child=0; child < size; child++, nodes_ptr++)
			Nodes.push_back(ForwardGraph[cur_n].nodes[child]);
	}

	for(unsigned int cur_n=1; cur_n<nameChangeCounter; cur_n++){
		unsigned int size = BackwardGraph[cur_n].nodes.size();
		BackwardG[cur_n].start = nodes_ptr;
		BackwardG[cur_n].numOfNodes = size;
		for(unsigned int child=0; child < size; child++, nodes_ptr++)
			Nodes.push_back(BackwardGraph[cur_n].nodes[child]);
	}

	// Delete the old Graph
	free(ForwardGraph);
	free(BackwardGraph);

	initialNameCounter=nameChangeCounter;
	evaluationF = (bool *)calloc(nameChangeCounter,sizeof(bool));
	evaluationB = (bool *)calloc(nameChangeCounter,sizeof(bool));

	preprocess2();		// Count children
	preprocess5();
	preprocess4();

	Query q;

	cout << "R" << endl << flush;

	cin.clear();
	cin >> c;

	unsigned int queries_counter = 0;
	while(cin >> c) {
		if(c == 'F') {
			pthread_mutex_lock (&mtx_work);
			while(rslts.size()!=multiQuery){
				pthread_cond_wait (&nowork , &mtx_work);
			}
			pthread_mutex_unlock (&mtx_work);

			for(int i=0; i<resultsCounter; i++){
				cout << results[i] << endl;
				//sumpaths += results[i];
			}

			cout << flush;
			cin.clear();
			rslts.clear();

			//sumqueries += resultsCounter;
			resultsCounter = 0;
			multiQuery = 0;

			continue;
		}

		cin >> a >> b;

		if(c == 'Q'){
			//query_pool.postWork<void>([a, b, versionCounter, resultsCounter] {  shortest_path(a, b, versionCounter, resultsCounter);  });
			//cerr << "kai edw eimaiste" << endl;
			// q.a=a; q.b=b; q.resultsCounter=resultsCounter; q.localVersion=versionCounter;
			if((evaluationF[nameChange[a]] && nameChange[a] && (nameChange[a]<initialNameCounter) && !ForwardG[nameChange[a]].addition)
					|| (evaluationB[nameChange[b]] && nameChange[b] && (nameChange[b]<initialNameCounter) && !BackwardG[nameChange[b]].addition)){
					queries_counter++;
					shortest_path(a,b,versionCounter,resultsCounter);
			}else{
				q.a=a; q.b=b; q.resultsCounter=resultsCounter; q.localVersion=versionCounter;
				pthread_mutex_lock (&mtx_q);
				queries.push(q);
				if(queries.size()>=1){
					pthread_cond_signal(&noemtpy);
				}
				pthread_mutex_unlock (&mtx_q);
				multiQuery++;
			}

			//cerr << "pairnaeii" << endl;
			versionCounter++;
			resultsCounter++;

		}else if(c == 'A'){
			//cerr << "kai sto add" << endl;
			if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
			if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

			add_edge(nameChange[a], nameChange[b], versionCounter);

			// info.op=c; info.a=nameChange[a]; info.b=nameChange[b];
			// operations.push_back(info);

			versionCounter++;

		}else{
			//cerr << "kai sto delete" << endl;
			if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
			if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

			delete_edge(nameChange[a], nameChange[b], versionCounter);

			// info.op=c; info.a=nameChange[a]; info.b=nameChange[b];
			// operations.push_back(info);

			versionCounter++;
		}
	}
	work_done=true;


	//output.close();

	exit(0);
}


