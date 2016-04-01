#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <unistd.h>
#include <algorithm>
#include <unordered_map>

#include "threadpool11/threadpool11.hpp"

using namespace std;
using namespace __gnu_cxx;

#define NUMOFTHREADS 23

#define MAXN ((size_t)1<<26)
#define NAMECHANGESIZE ((size_t)1<<28)

typedef unsigned int Node;
typedef unsigned short int Visitor;

typedef struct GraphNode {
   bool addition;
   bool deletion;
   vector<Node> nodes;
   int children;

   GraphNode(){
	   addition=false;
	   deletion=false;
	   children = 0;
   }
} GraphNode;


struct MyPair{
	Node b;
	int version;
};

struct Operation_info{
	char op;
	Node a,b;
};

// Variables for the Queues
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

	 cpu_set_t cpuset;
	 CPU_ZERO(&cpuset);
	 CPU_SET(threadCounter, &cpuset);

	 int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
	 if (rc != 0)
		 std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";


	threadCounter++;
	mtx.unlock();

	sleep(1);
}

void shortest_path(Node tempa, Node tempb, int localVersion, int resultsCounter) {
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

	int thread_id = threadIds.find(pthread_self())->second;

	Visitor *visited = visited_global[thread_id];
	Node **fQueue = fQueue_global[thread_id];
	Node **bQueue = bQueue_global[thread_id];

	visitedCounter_global[thread_id] += 2;
	const unsigned int visitedCounter = visitedCounter_global[thread_id];


	// Initialize Queues and Variables
	unsigned int fChildrenCount = ForwardGraph[a].nodes.size();
	unsigned int fCurrentNodes = 1;
	unsigned int bChildrenCount = BackwardGraph[b].nodes.size();
	unsigned int bCurrentNodes = 1;

	unsigned int fQueuePointer = 0;
	unsigned int bQueuePointer = 0;
	fQueue[fQueuePointer][0] = a;
	bQueue[bQueuePointer][0] = b;

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

		if(fChildrenCount <= bChildrenCount){	// Move forward, there are less children there
			fChildrenCount = 0;
			fGraphDistance++;	// Going to the next distance

			// Reading all the children from the nodes in the queue
			for(unsigned int n=0; n < fCurrentNodes; n++){
				const unsigned int currentFather = fQueue[fQueuePointer][n];
				fChildrenCount += ForwardGraph[currentFather].children;

				// Check if this currentFather node has any deletions done at him
				if(ForwardGraph[currentFather].deletion == 0){	// Just check the children nodes

					// Reading the children of the current node in the queue
					const unsigned int end = ForwardGraph[currentFather].nodes.size();
					for(unsigned int j=0; j < end; j++){
						const unsigned int currentChild = ForwardGraph[currentFather].nodes[j];
						if(visited[currentChild] >= visitedCounter){	// Explored by the other side
							if(visited[currentChild] == visitedCounter+1){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
							continue;
						}
						visited[currentChild] = visitedCounter;
						fQueue[fQueuePointer^1][iterator] = currentChild;
						iterator++;
					}
				}else{		// Deletion == 1
					// We have to check at every child if it was deleted
					const unsigned int end = ForwardGraph[currentFather].nodes.size();
					for(unsigned int i=0; i < end; i++){
						const unsigned int currentChild = ForwardGraph[currentFather].nodes[i];

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
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
							continue;
						}
						visited[currentChild] = visitedCounter;
						fQueue[fQueuePointer^1][iterator] = currentChild;
						iterator++;
					}
				}

				// Check if there were any additions in the current father node
				if(ForwardGraph[currentFather].addition == 1){
					// For every node added
					if(ForwardGraph[currentFather].deletion == 0){

						const unsigned int end = Forward_add[currentFather].size();
						for(unsigned int i = 0; i < end; i++){
							// Just add the nodes to the queue
							const unsigned int child = Forward_add[currentFather][i].b;

							// The addition must have appeared before the query
							if(Forward_add[currentFather][i].version > localVersion)
								continue;

							if(visited[child] >= visitedCounter){	// Explored by the other side
								if(visited[child] == visitedCounter+1){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
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
							const unsigned int child = Forward_add[currentFather][i].b;

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
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
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
				results[resultsCounter] = -1;
				return;
			}
			fCurrentNodes = iterator;
			fQueuePointer = fQueuePointer^1;

		}else{			// bChildrenCounter > fChildrenCounter
			bChildrenCount = 0;
			bGraphDistance++;	// Going to the next distance

			// Reading all the children from the nodes in the queue
			for(unsigned int n=0; n < bCurrentNodes; n++){
				const unsigned int currentFather = bQueue[bQueuePointer][n];
				bChildrenCount += BackwardGraph[currentFather].children;

				// Check if this currentFather node has any deletions done at him
				if(BackwardGraph[currentFather].deletion == 0){	// Just check the children nodes

					// Reading the children of the current node in the queue
					const unsigned int end = BackwardGraph[currentFather].nodes.size();
					for(unsigned int j=0; j < end; j++){
						const unsigned int currentChild = BackwardGraph[currentFather].nodes[j];
						if(visited[currentChild] >= visitedCounter){	// Explored by the other side
							if(visited[currentChild] == visitedCounter){ 	// Found the minimum distance!
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
							continue;
						}
						visited[currentChild] = visitedCounter+1;
						bQueue[bQueuePointer^1][iterator] = currentChild;
						iterator++;
					}
				}else{		// Deletion == 1
					// We have to check at every child if it was deleted
					const unsigned int end = BackwardGraph[currentFather].nodes.size();
					for(unsigned int i=0; i < end; i++){
						const unsigned int currentChild = BackwardGraph[currentFather].nodes[i];

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
								results[resultsCounter] = fGraphDistance + bGraphDistance;
								return;
							}
							continue;
						}
						visited[currentChild] = visitedCounter+1;
						bQueue[bQueuePointer^1][iterator] = currentChild;
						iterator++;
					}
				}

				// Check if there were any additions in the current father node
				if(BackwardGraph[currentFather].addition == 1){
					// For every node added
					if(BackwardGraph[currentFather].deletion == 0){

						const unsigned int end = Backward_add[currentFather].size();
						for(unsigned int i = 0; i < end; i++){
							// Just add the nodes to the queue
							const unsigned int child = Backward_add[currentFather][i].b;

							// The addition must have appeared before the query
							if(Backward_add[currentFather][i].version > localVersion)
								continue;

							if(visited[child] >= visitedCounter){	// Explored by the other side
								if(visited[child] == visitedCounter){ 	// Found the minimum distance!
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
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
							const unsigned int child = Backward_add[currentFather][i].b;

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
									results[resultsCounter] = fGraphDistance + bGraphDistance;
									return;
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
				results[resultsCounter] = -1;
				return;
			}
			bCurrentNodes = iterator;
			bQueuePointer = bQueuePointer^1;
		}
	}
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

void add_edge(Node a, Node b, int versionCounter){
	MyPair fPair, bPair;
	fPair.version=versionCounter;
	fPair.b=b;

	bPair.version=versionCounter;
	bPair.b=a;

	Forward_add[a].push_back(fPair);
	ForwardGraph[a].addition=true;

	Backward_add[b].push_back(bPair);
	BackwardGraph[b].addition=true;
}

void delete_edge(Node a, Node b, int versionCounter){
	MyPair fPair, bPair;
	fPair.version=versionCounter;
	fPair.b=b;

	bPair.version=versionCounter;
	bPair.b=a;

	Forward_del[a].push_back(fPair);
	ForwardGraph[a].deletion=true;

	Backward_del[b].push_back(bPair);
	BackwardGraph[b].deletion=true;
}

// Count the children of a node
void preprocess2(){
	for(unsigned int i=1; i < nameChangeCounter; i++)
		for(unsigned int j=0; j < ForwardGraph[i].nodes.size(); j++)
			ForwardGraph[i].children += ForwardGraph[ForwardGraph[i].nodes[j]].nodes.size();

	for(unsigned int i=1; i < nameChangeCounter; i++)
		for(unsigned int j=0; j < BackwardGraph[i].nodes.size(); j++)
			BackwardGraph[i].children += BackwardGraph[BackwardGraph[i].nodes[j]].nodes.size();

}

// Initialize vectors
void preprocess3(){
	for(unsigned int i=0; i<nameChangeCounter; i++){
		Forward_del[i].reserve(2);
		Forward_add[i].reserve(2);
		Backward_del[i].reserve(2);
		Backward_add[i].reserve(2);
	}
}

int main() {
	ios_base::sync_with_stdio(false);
	cin.tie(nullptr);

	int versionCounter=0;
	Node a,b;
	char c;
	int resultsCounter = 0;

	vector<Operation_info> operations(10000);
	Operation_info info;

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





	while(cin >> a >> b){
		if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
		if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
		add_edge_final(nameChange[a],nameChange[b]);
	}

	preprocess();			// Sorting the biggest Nodes to the front of the lists
	preprocess2();			// Count children

	// Creating Threads
	//size_t numOfThreads = 1;
	threadpool11::Pool query_pool;
	query_pool.setWorkerCount(NUMOFTHREADS);

	for(int i=0; i<NUMOFTHREADS+1; i++)
		query_pool.postWork<void>(initThread);
	query_pool.waitAll();

	cout << "R" << endl << flush;

	cin.clear();
	cin >> c;

	while(cin >> c) {
		if(c == 'F') {
			query_pool.waitAll();
			// Print the results
			for(int i=0; i<resultsCounter; i++)
				cout << results[i] << endl;

			cout << flush;
			cin.clear();

			resultsCounter = 0;
			/*for(vector<Operation_info>::iterator it=operations.begin(); it!=operations.end(); it++){
				if(it->op=='A'){
					add_edge_final(it->a,it->b);
					if(Forward_add[a].size()){
						Forward_add[a].clear();
					}
					// here must be changed the flags of the struct ForwardGraph and BackwardGraph
					if(Backward_add[b].size()){
						Backward_add[b].clear();
					}
					ForwardGraph[a].addition=false;
					BackwardGraph[b].addition=false;
				}
				else{
					delete_edge_final(it->a,it->b);
					if(Forward_del[a].size()){
						Forward_del[a].clear();
					}
					if(Backward_del[b].size()){
						Backward_del[b].clear();
					}
					BackwardGraph[b].deletion=false;
					ForwardGraph[a].addition=false;
				}
			}
			operations.clear();*/
			continue;
		}

		cin >> a >> b;

		if(c == 'Q'){
			query_pool.postWork<void>([a, b, versionCounter, resultsCounter] {  shortest_path(a, b, versionCounter, resultsCounter);  });
			versionCounter++;
			resultsCounter++;

		}else if(c == 'A'){
			if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
			if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

			add_edge(nameChange[a], nameChange[b], versionCounter);

			info.op=c; info.a=nameChange[a]; info.b=nameChange[b];
			operations.push_back(info);

			versionCounter++;

		}else{
			if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
			if(!nameChange[b]) nameChange[b] = nameChangeCounter++;

			delete_edge(nameChange[a], nameChange[b], versionCounter);

			info.op=c; info.a=nameChange[a]; info.b=nameChange[b];
			operations.push_back(info);

			versionCounter++;
		}
	}

	return 0;
}

