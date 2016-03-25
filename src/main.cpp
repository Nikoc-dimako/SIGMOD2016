#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <unistd.h>
#include <algorithm>
#include <unordered_map>
#include <mutex>

#include "threadpool11/threadpool11.hpp"

using namespace std;
using namespace __gnu_cxx;

#define NUMOFTHREADS 23

#define MAXN ((size_t)1<<24)
#define FORWARD 0
#define BACKWARD 1
#define ADD 0
#define DELETE 1
#define MAXTH 2
typedef unsigned int Node;


typedef struct GraphNode {
   bool addition;
   bool deletion;
   vector<Node> nodes;

   GraphNode(){	   addition=deletion=false; }
} GraphNode;


struct Operation_info{
	char op;
	Node a,b;
};


struct Child_info{
	Node a;
	vector<Node> children;
};


struct MyPair{
	Node b;
	size_t version;
};

typedef struct SP_params {
   Node a;
   Node b;
   size_t localVersion;
   int resultsCounter;
} SP_params;


// Variables for the Queues
GraphNode *ForwardGraph, *BackwardGraph;
Node **visited_global;
Node **dist_global;
Node **fQueue_global;
Node **bQueue_global;
unsigned int *visitedCounter_global;


//Variables for the node's nameChange
Node *nameChange;
Node *nameReorder;
unsigned int nameChangeCounter = 1;


// Variables for the multiversion
size_t versionCounter=0;
vector<MyPair> *Forward_add,*Forward_del,*Backward_add,*Backward_del;
int results[100000];
SP_params *parameters;
short int threadCounter = 0;
unordered_map <pthread_t, short int> threadIds;
mutex mtx;

void getThreadId(){
	mtx.lock();
	unordered_map<pthread_t, short int>::const_iterator um_it = threadIds.find(pthread_self());
	cerr << "Thread id " <<  pthread_self() << " -> "<< threadCounter << endl;
	threadIds.emplace(pthread_self(), threadCounter);
	threadCounter++;
	mtx.unlock();

	sleep(1);
}


void shortest_path(SP_params parameters) {
	int thread_id = threadIds.find(pthread_self())->second;

	Node a = parameters.a;
	Node b = parameters.b;
	size_t localVersion = parameters.localVersion;
	int resultsCounter = parameters.resultsCounter;

	if(a == b){
		results[resultsCounter] = 0;
		return;
	}

	Node *visited = visited_global[thread_id];
	Node *dist = dist_global[thread_id];
	Node *fQueue = fQueue_global[thread_id];
	Node *bQueue = bQueue_global[thread_id];

	visitedCounter_global[thread_id] += 2;
	unsigned int visitedCounter = visitedCounter_global[thread_id];


	// Initialize Queues and Variables
	unsigned int fChildrenCount = ForwardGraph[a].nodes.size();
	bool fChildrenAdded = ForwardGraph[a].addition;
	unsigned int fCurrentNodes = 1;
	unsigned int bChildrenCount = BackwardGraph[b].nodes.size();
	bool bChildrenAdded = BackwardGraph[b].addition;
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
	dist[a] = 0;
	dist[b] = 0;

	unsigned int iterator;
	Node currentFather;
	Node currentChild;

	// The queue saves the currently explored nodes
	// The children counters contain the number of children that the already explored nodes have.

	// Main Loop
	while(1){
		iterator = 0;

		if(fChildrenCount <= bChildrenCount){	// Move forward, there are less children there
			fChildrenCount = 0;
			fChildrenAdded = false;
			fGraphDistance++;	// Going to the next distance

			// Reading all the children from the nodes in the queue
			for(unsigned int i=0; i < fCurrentNodes; i++){
				currentFather = fQueue[i+fQueuePointer*MAXN];

				// Check if this currentFather node has any deletions done at him
				if(ForwardGraph[currentFather].deletion == 0){	// Just check the children nodes

					// Reading the children of the current node in the queue
					for(unsigned int j=0; j < ForwardGraph[currentFather].nodes.size(); j++){
						currentChild = ForwardGraph[currentFather].nodes[j];
						if(visited[currentChild] >= visitedCounter){	// Explored by the other side
							if(visited[currentChild] == visitedCounter+1){ 	// Found the minimum distance!
								results[resultsCounter] = dist[currentChild] + dist[currentFather] + 1;
								return;
							}
							continue;
						}
						visited[currentChild] = visitedCounter;
						dist[currentChild] = fGraphDistance;
						fChildrenCount += ForwardGraph[currentChild].nodes.size();	// Counting the children of the next stage
						fChildrenAdded = fChildrenAdded || ForwardGraph[currentChild].addition;
						fQueue[iterator+(fQueuePointer^1)*MAXN] = currentChild;
						iterator++;
					}
				}else{		// Deletion == 1
					// We have to check at every child if it was deleted
					for(unsigned int j=0; j < ForwardGraph[currentFather].nodes.size(); j++){
						currentChild = ForwardGraph[currentFather].nodes[j];

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
								results[resultsCounter] = dist[currentChild] + dist[currentFather] + 1;
								return;
							}
							continue;
						}
						visited[currentChild] = visitedCounter;
						dist[currentChild] = fGraphDistance;
						fChildrenCount += ForwardGraph[currentChild].nodes.size();	// Counting the children of the next stage
						fChildrenAdded = fChildrenAdded || ForwardGraph[currentChild].addition;
						fQueue[iterator+(fQueuePointer^1)*MAXN] = currentChild;
						iterator++;
					}
				}

				// Check if there were any additions in the current father node
				if(ForwardGraph[currentFather].addition == 1){
					// For every node added
					if(ForwardGraph[currentFather].deletion == 0){
						for(std::vector<MyPair> ::iterator it=Forward_add[currentFather].begin(); it!=Forward_add[currentFather].end(); it++){

							// Just add the nodes to the queue
							Node child=it->b;

							// The addition must have appeared before the query
							if(it->version > localVersion)
								continue;

							if(visited[child] >= visitedCounter){	// Explored by the other side
								if(visited[child] == visitedCounter+1){ 	// Found the minimum distance!
									results[resultsCounter] = dist[child] + dist[currentFather] + 1;
									return;
								}
								continue;
							}
							visited[child]=visitedCounter;
							dist[child]=fGraphDistance;
							fChildrenCount += ForwardGraph[child].nodes.size();
							fChildrenAdded = fChildrenAdded || ForwardGraph[child].addition;
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
									results[resultsCounter] = dist[child] + dist[currentFather] + 1;
									return;
								}
								continue;
							}
							visited[child]=visitedCounter;
							dist[child]=fGraphDistance;
							fChildrenCount += ForwardGraph[child].nodes.size();
							fChildrenAdded = fChildrenAdded || ForwardGraph[child].addition;
							fQueue[iterator+(fQueuePointer^1)*MAXN]=child;
							iterator++;
						}
					}
				}
			}

			if(fChildrenCount == 0 && !fChildrenAdded){
				results[resultsCounter] = -1;
				return;
			}
			fCurrentNodes = iterator;
			fQueuePointer = fQueuePointer^1;

		}else{			// bChildrenCounter > fChildrenCounter
			bChildrenCount = 0;
			bChildrenAdded = false;
			bGraphDistance++;	// Going to the next distance

			// Reading all the children from the nodes in the queue
			for(unsigned int i=0; i < bCurrentNodes; i++){
				currentFather = bQueue[i+bQueuePointer*MAXN];

				// Check if this currentFather node has any deletions done at him
				if(BackwardGraph[currentFather].deletion == 0){	// Just check the children nodes

					// Reading the children of the current node in the queue
					for(unsigned int j=0; j < BackwardGraph[currentFather].nodes.size(); j++){
						currentChild = BackwardGraph[currentFather].nodes[j];
						if(visited[currentChild] >= visitedCounter){	// Explored by the other side
							if(visited[currentChild] == visitedCounter){ 	// Found the minimum distance!
								results[resultsCounter] = dist[currentChild] + dist[currentFather] + 1;
								return;
							}
							continue;
						}
						visited[currentChild] = visitedCounter+1;
						dist[currentChild] = bGraphDistance;
						bChildrenCount += BackwardGraph[currentChild].nodes.size();	// Counting the children of the next stage
						bChildrenAdded = bChildrenAdded || BackwardGraph[currentChild].addition;
						bQueue[iterator+(bQueuePointer^1)*MAXN] = currentChild;
						iterator++;
					}
				}else{		// Deletion == 1
					// We have to check at every child if it was deleted
					for(unsigned int j=0; j < BackwardGraph[currentFather].nodes.size(); j++){
						currentChild = BackwardGraph[currentFather].nodes[j];

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
								results[resultsCounter] = dist[currentChild] + dist[currentFather] + 1;
								return;
							}
							continue;
						}
						visited[currentChild] = visitedCounter+1;
						dist[currentChild] = bGraphDistance;
						bChildrenCount += BackwardGraph[currentChild].nodes.size();	// Counting the children of the next stage
						bChildrenAdded = bChildrenAdded || BackwardGraph[currentChild].addition;
						bQueue[iterator+(bQueuePointer^1)*MAXN] = currentChild;
						iterator++;
					}
				}

				// Check if there were any additions in the current father node
				if(BackwardGraph[currentFather].addition == 1){
					// For every node added
					for(std::vector<MyPair> ::iterator it=Backward_add[currentFather].begin(); it!=Backward_add[currentFather].end(); it++){
						if(BackwardGraph[currentFather].deletion == 0){

							// Just add the nodes to the queue
							Node child=it->b;

							// The addition must have appeared before the query
							if(it->version > localVersion)
								continue;

							if(visited[child] >= visitedCounter){	// Explored by the other side
								if(visited[child] == visitedCounter){ 	// Found the minimum distance!
									results[resultsCounter] = dist[child] + dist[currentFather] + 1;
									return;
								}
								continue;
							}
							visited[child]=visitedCounter+1;
							dist[child]=bGraphDistance;
							bChildrenCount += BackwardGraph[child].nodes.size();
							bChildrenAdded = bChildrenAdded || BackwardGraph[child].addition;
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
									results[resultsCounter] = dist[child] + dist[currentFather] + 1;
									return;
								}
								continue;
							}
							visited[child]=visitedCounter+1;
							dist[child]=bGraphDistance;
							bChildrenCount += BackwardGraph[child].nodes.size();
							bChildrenAdded = bChildrenAdded || BackwardGraph[child].addition;
							bQueue[iterator+(bQueuePointer^1)*MAXN]=child;
							iterator++;
						}
					}
				}
			}

			if(bChildrenCount == 0 && !bChildrenAdded){
				results[resultsCounter] = -1;
				return;
			}
			bCurrentNodes = iterator;
			bQueuePointer = bQueuePointer^1;
		}
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

void delete_edge(Node a,Node b) {
  if(!_delete_edge(b,ForwardGraph[a].nodes) || !_delete_edge(a,BackwardGraph[b].nodes) )
	  return;
}

bool sortF(Node a,Node b) { return (ForwardGraph[a].nodes.size() > ForwardGraph[b].nodes.size()); }
bool sortB(Node a,Node b) { return (BackwardGraph[a].nodes.size() > BackwardGraph[b].nodes.size()); }

void preprocess() {
  for(unsigned int a=1; a < nameChangeCounter; a++) {
    sort(ForwardGraph[a].nodes.begin(),ForwardGraph[a].nodes.end(),sortF);
    sort(BackwardGraph[a].nodes.begin(),BackwardGraph[a].nodes.end(),sortB);
  }
}

// Order the list of vectors, the one with most children first
void add_edge_final(Node a, Node b) {
  ForwardGraph[a].nodes.push_back(b);
  BackwardGraph[b].nodes.push_back(a);
}

void delete_edge_final(Node a,Node b) {
	if(!_delete_edge(b,ForwardGraph[a].nodes) || !_delete_edge(a,BackwardGraph[b].nodes) )
		return;
}

void multiversion_changing(int graph,int operation,Node a,Node b){
	MyPair pair;

	pair.version=versionCounter;
	pair.b=b;

	if(graph==FORWARD){
		if(operation==ADD){
			Forward_add[a].push_back(pair);
			ForwardGraph[a].addition=true;
		}
		else{
			Forward_del[a].push_back(pair);
			ForwardGraph[a].deletion=true;
		}
	}
	else{
		if(operation==ADD){
			Backward_add[a].push_back(pair);
			BackwardGraph[a].addition=true;
		}
		else{
			Backward_del[a].push_back(pair);
			BackwardGraph[a].deletion=true;
		}
	}
}

int main() {
	ios_base::sync_with_stdio(false);
	cin.tie(nullptr);

	Node a,b;
	char c;
	vector<Operation_info> operations;
	Operation_info info;
	int resultsCounter = 0;

	ForwardGraph = (GraphNode*)calloc(MAXN, sizeof(GraphNode));
	BackwardGraph = (GraphNode*)calloc(MAXN, sizeof(GraphNode));
	nameChange = (Node*)calloc(((size_t)1<<28), sizeof(Node));
	Forward_del = (vector<MyPair>*)calloc(MAXN, sizeof(vector<MyPair>));
	Forward_add = (vector<MyPair>*)calloc(MAXN, sizeof(vector<MyPair>));
	Backward_del = (vector<MyPair>*)calloc(MAXN, sizeof(vector<MyPair>));
	Backward_add = (vector<MyPair>*)calloc(MAXN, sizeof(vector<MyPair>));

	// Creating structures for each thread
	fQueue_global = (Node**)calloc(NUMOFTHREADS, sizeof(Node*));
	bQueue_global = (Node**)calloc(NUMOFTHREADS, sizeof(Node*));
	visited_global = (Node**)calloc(NUMOFTHREADS, sizeof(Node*));
	dist_global = (Node**)calloc(NUMOFTHREADS, sizeof(Node*));
	visitedCounter_global = (unsigned int*)calloc(NUMOFTHREADS, sizeof(unsigned int));
	parameters = (SP_params*)calloc(NUMOFTHREADS, sizeof(SP_params));
	for(int i=0; i<NUMOFTHREADS; i++){
		fQueue_global[i] = (Node*)calloc(MAXN*2, sizeof(Node));
		bQueue_global[i] = (Node*)calloc(MAXN*2, sizeof(Node));
		visited_global[i] = (Node*)calloc(MAXN, sizeof(Node));
		dist_global[i] = (Node*)calloc(MAXN, sizeof(Node));
	}

	while(cin >> a >> b){
		if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
		if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
		add_edge_final(nameChange[a],nameChange[b]);
	}

	SP_params parameters;

	preprocess();			// Sorting the biggest Nodes to the front of the lists

	// Creating Threads
	threadpool11::Pool pool;
	pool.setWorkerCount(NUMOFTHREADS);

	for(int i=0; i<NUMOFTHREADS; i++)
		pool.postWork<void>(getThreadId);
	pool.waitAll();

	cout << "R" << endl << flush;

	cin.clear();
	cin >> c;

	while(cin >> c) {
		if(c == 'F') {
			pool.waitAll();
			// Print the results
			for(int i=0; i<resultsCounter; i++)
				cout << results[i] << endl;

			cout << flush;
			cin.clear();

			/*for(vector<Operation_info>::iterator it=operations.begin(); it!=operations.end(); it++){
				if(it->op=='A'){
					add_edge_final(it->a,it->b);

					Forward_add[a].clear();
					Backward_add[b].clear();

					// here must be changed the flags of the struct ForwardGraph and BackwardGraph
					ForwardGraph[a].addition=false;
					BackwardGraph[b].addition=false;

				}else{
					delete_edge_final(it->a,it->b);

					Forward_del[a].clear();
					Backward_del[b].clear();

					BackwardGraph[b].deletion=false;
					ForwardGraph[a].addition=false;
				}
			}
			operations.clear();*/
			resultsCounter = 0;

			continue;
		}

		cin >> a >> b;

		if(c == 'Q'){
			if(nameChange[a] == 0 || nameChange[b] == 0){
				if(a != b)
					results[resultsCounter] = -1;
				else
					results[resultsCounter] = 0;
				resultsCounter++;
				continue;
			}

			parameters.a = nameChange[a];
			parameters.b = nameChange[b];
			parameters.localVersion = versionCounter;
			parameters.resultsCounter = resultsCounter;

			pool.postWork<void>([parameters] {  shortest_path(parameters);  });
			versionCounter++;
			resultsCounter++;

		}else if(c == 'A'){
			if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
			if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
			//info.op='A'; info.a=nameChange[a]; info.b=nameChange[b];
			//operations.push_back(info);	v2

			multiversion_changing(FORWARD,ADD,nameChange[a],nameChange[b]);
			multiversion_changing(BACKWARD,ADD,nameChange[b],nameChange[a]);

			versionCounter++;

		}else if(c == 'D') {
			if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
			if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
			//info.op='D'; info.a=nameChange[a]; info.b=nameChange[b];
			//operations.push_back(info);	v2

			multiversion_changing(FORWARD,DELETE,nameChange[a],nameChange[b]);
			multiversion_changing(BACKWARD,DELETE,nameChange[b],nameChange[a]);

			versionCounter++;
		}
	}
  return 0;
}