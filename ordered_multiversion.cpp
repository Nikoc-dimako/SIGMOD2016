/*
 * main.cpp
 *
 *  Created on: Mar 6, 2016
 *      Author: thanasis
 */

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <ext/hash_map>
#include <unordered_map>
#include <hash_map>
#include <iostream>
#include <unistd.h>
#include <algorithm>
#include <map>


using namespace std;
using namespace __gnu_cxx;


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

   GraphNode(){ addition=deletion=false; }
} GraphNode;

// Variables for the Queues
GraphNode *ForwardGraph, *BackwardGraph;
Node *visited, *dist;
Node *fQueue;
Node *bQueue;
int fd[MAXTH][2];
unsigned int visitedCounter = 0;


struct Operation_info{
	char op;
	Node a,b;
};


struct Child_info{
	Node a;
	vector<Node> children;
};

struct Graph_info{
	//int id; 8a to steilw san parametro
	unsigned int query,version;
	Node a,b;

};

struct MyPair{
	Node b;
	size_t version;
};




//Variables for the node's nameChange
Node *nameChange;
Node *nameReorder;
unsigned int nameChangeCounter = 1;

size_t versionStart=0,versionCounter=0;
vector<MyPair> *Forward_add,*Forward_del,*Backward_add,*Backward_del;
//std::vector<Node> deletions[2],additions[2];

void construct_the_right_children(int graph,Node a,size_t version,vector<Node>& children){
	vector<MyPair>::iterator itA,itD;
	bool flag;

	if(graph==FORWARD){
		if(ForwardGraph[a].addition && ForwardGraph[a].deletion){
			children=ForwardGraph[a].nodes;
			itA=Forward_add[a].begin();
			itD=Forward_del[a].begin();

			while(itA!=Forward_add[a].end() && itD!=Forward_del[a].end()){
				if(itA->version<version && itD->version<version){
					if(itA->version<itD->version){
						children.push_back(itA->b);
						itA++;
					}
					else{
						std::vector<Node>::iterator it;
						it=find(children.begin(),children.end(),itD->b); // here we need a hash for deletions 
						if(it!=children.end()){
							children.erase(it);
						}
						itD++;
					}
				}
				else{
					break;
				}
			}

			if(itA==Forward_add[a].end() && itD==Forward_del[a].end()){
				return;
			}
			else{
				while(itA!=Forward_add[a].end()){
					if(itA->version<version){
						children.push_back(itA->b);
						itA++;
					}
					else{
						break;
					}
				}
				while(itD!=Forward_del[a].end()){
					if(itD->version<version){
						std::vector<Node>::iterator it;
						it=find(children.begin(),children.end(),itD->b); // here we need a hash for deletions 
						if(it!=children.end()){
							children.erase(it);
						}
						itD++;
					}
					else{
						return;
					}
				}
			}
		}
		else if(ForwardGraph[a].addition){
			children=ForwardGraph[a].nodes;
			itA=Forward_add[a].begin();
			while(itA!=Forward_add[a].end()){
				if(itA->version<version){
					children.push_back(itA->b);
					itA++;
				}
				else{
					break;
				}
			}
		}
		else if(ForwardGraph[a].deletion){
			children=ForwardGraph[a].nodes;
			itD=Forward_del[a].begin();
			while(itD!=Forward_del[a].end()){
				if(itD->version<version){
					std::vector<Node>::iterator it;
					it=find(children.begin(),children.end(),itD->b); // here we need a hash for deletions 
					if(it!=children.end()){
						children.erase(it);
					}
					itD++;
				}
				else{
					return;
				}
			}
		}
		else{
			children=ForwardGraph[a].nodes;
		}
	}
	else{
		if(BackwardGraph[a].addition && BackwardGraph[a].deletion){
			children=BackwardGraph[a].nodes;
			itD=Backward_del[a].begin();
			itA=Backward_add[a].begin();

			while(itA!=Backward_add[a].end() && itD!=Backward_del[a].end()){
				if(itA->version<version && itD->version<version){
					if(itA->version<itD->version){
						children.push_back(itA->b);
						itA++;
					}
					else{
						std::vector<Node>::iterator it;
						it=find(children.begin(),children.end(),itD->b); // here we need a hash for deletions 
						if(it!=children.end()){
							children.erase(it);
						}
						itD++;
					}
				}
				else{
					break;
				}
			}

			if(itA==Backward_add[a].end() && itD==Backward_del[a].end()){
				return ;
			}
			else{
				while(itA!=Backward_add[a].end()){
					if(itA->version<version){
						children.push_back(itA->b);
						itA++;
					}
					else{
						break;
					}
				}
				while(itD!=Backward_del[a].end()){
					if(itD->version<version){
						std::vector<Node>::iterator it;
						it=find(children.begin(),children.end(),itD->b); // here we need a hash for deletions 
						if(it!=children.end()){
							children.erase(it);
						}
						itD++;
					}
					else{
						return;
					}
				}
			}
		}
		else if(BackwardGraph[a].addition){
			children=BackwardGraph[a].nodes;
			itA=Backward_add[a].begin();
			while(itA!=Backward_add[a].end()){
				if(itA->version<version){
					children.push_back(itA->b);
					itA++;
				}
				else{
					break;
				}
			}
		}
		else if(BackwardGraph[a].deletion){
			children=BackwardGraph[a].nodes;
			itD=Backward_del[a].begin();
			while(itD!=Backward_del[a].end()){
				if(itD->version<version){
					std::vector<Node>::iterator it;
					it=find(children.begin(),children.end(),itD->b); // here we need a hash for deletions 
					if(it!=children.end()){
						children.erase(it);
					}
					itD++;
				}
				else{
					return;
				}
			}
		}
		else{
			children=BackwardGraph[a].nodes;
		}
	}

}

int shortest_path(Node a, Node b,size_t localVersion) {
	if(a == b)
		return 0;

	visitedCounter += 2;

	// Initialize Queues and Variables
	unsigned int fChildrenCount = ForwardGraph[a].nodes.size();
	unsigned int fCurrentNodes = 1;
	unsigned int bChildrenCount = BackwardGraph[b].nodes.size();
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
	//cout << "Komple prin th loopa me " << a << " " << b << endl;
	while(1){
		iterator = 0;

		if(fChildrenCount <= bChildrenCount){	// Move forward, there are less children there
			fChildrenCount = 0;
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
							if(visited[currentChild] == visitedCounter+1) 	// Found the minimum distance!
								return dist[currentChild] + dist[currentFather] + 1;
							continue;
						}
						visited[currentChild] = visitedCounter;
						dist[currentChild] = fGraphDistance;
						fChildrenCount += ForwardGraph[currentChild].nodes.size();	// Counting the children of the next stage
						fQueue[iterator+(fQueuePointer^1)*MAXN] = currentChild;
						iterator++;
					}
				}else{
					// We have to check at every child if it was deleted
					for(unsigned int j=0; j < ForwardGraph[currentFather].nodes.size(); j++){
						currentChild = ForwardGraph[currentFather].nodes[j];

						// Check if the node was deleted
						// If it was continue to the next child;
						// TODO

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
							if(visited[currentChild] == visitedCounter+1) 	// Found the minimum distance!
								return dist[currentChild] + dist[currentFather] + 1;
							continue;
						}
						visited[currentChild] = visitedCounter;
						dist[currentChild] = fGraphDistance;
						fChildrenCount += ForwardGraph[currentChild].nodes.size();	// Counting the children of the next stage
						fQueue[iterator+(fQueuePointer^1)*MAXN] = currentChild;
						iterator++;
					}
				}

				// Check if there were any additions in the current father node
				if(ForwardGraph[currentFather].addition == 1){
					// For every node added
					for(std::vector<MyPair> ::iterator it=Forward_add[currentFather].begin(); it!=Forward_add[currentFather].end(); it++){
						if(ForwardGraph[currentFather].deletion == 0){

							// Just add the nodes to the queue
							// TODO
							Node child=it->b;
							visited[child]=visitedCounter;
							dist[child]=fGraphDistance;
							fChildrenCount += ForwardGraph[child].nodes.size();
							fQueue[iterator+(fQueuePointer^1)*MAXN]=child;
							iterator++;
						}
						else{

							// You have to check if it was deleted first
							// TODO
							Node child=it->b;
							std::vector<MyPair>:: iterator itD;
							for(itD=Forward_del[currentFather].begin(); itD!=Forward_del[currentFather].end(); itD++){
								if(itD->b==child){
									break;
								}
							}
							if(itD!=Forward_del[currentFather].end() && itD->version<localVersion){
								continue;
							}
							visited[child]=visitedCounter;
							dist[child]=fGraphDistance;
							fChildrenCount += ForwardGraph[child].nodes.size();
							fQueue[iterator+(fQueuePointer^1)*MAXN]=child;
							iterator++;


						}
					}
				}
			}

			if(fChildrenCount == 0)
				return -1;
			fCurrentNodes = iterator;
			fQueuePointer = fQueuePointer^1;
		}
		else{
			bChildrenCount = 0;
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
							if(visited[currentChild] == visitedCounter) 	// Found the minimum distance!
								return dist[currentChild] + dist[currentFather] + 1;
							continue;
						}
						visited[currentChild] = visitedCounter+1;
						dist[currentChild] = bGraphDistance;
						bChildrenCount += BackwardGraph[currentChild].nodes.size();	// Counting the children of the next stage
						bQueue[iterator+(bQueuePointer^1)*MAXN] = currentChild;
						iterator++;
					}
				}else{
					// We have to check at every child if it was deleted
					for(unsigned int j=0; j < BackwardGraph[currentFather].nodes.size(); j++){
						currentChild = BackwardGraph[currentFather].nodes[j];

						// Check if the node was deleted
						// If it was continue to the next child;
						// TODO

						std::vector<MyPair>:: iterator it;
						for(it=Backward_del[currentFather].begin(); it!=Backward_del[currentFather].end(); it++){
							if(it->b==currentChild){
								break;
							}
						}
						if(it!=Backward_del[currentFather].end() || it->version<localVersion){
							continue;
						}

						if(visited[currentChild] >= visitedCounter){	// Explored by the other side
							if(visited[currentChild] == visitedCounter) 	// Found the minimum distance!
								return dist[currentChild] + dist[currentFather] + 1;
							continue;
						}
						visited[currentChild] = visitedCounter+1;
						dist[currentChild] = bGraphDistance;
						bChildrenCount += BackwardGraph[currentChild].nodes.size();	// Counting the children of the next stage
						bQueue[iterator+(bQueuePointer^1)*MAXN] = currentChild;
						iterator++;
					}
				}

				// Check if there were any additions in the current father node
				if(BackwardGraph[currentFather].addition == 1){
					// For every node added
					//cout << "Mpaine edw sto addition apo BackwardGraph me "  << a << " " << b << endl;
					for(std::vector<MyPair> ::iterator it=Backward_add[currentFather].begin(); it!=Backward_add[currentFather].end(); it++){
						if(BackwardGraph[currentFather].deletion == 0){

							// Just add the nodes to the queue
							// TODO
							Node child=it->b;
							visited[child]=visitedCounter+1;
							dist[child]=bGraphDistance;
							bChildrenCount += BackwardGraph[child].nodes.size();
							bQueue[iterator+(bQueuePointer^1)*MAXN]=child;
							iterator++;
							//cout << "Edw komple apo BackwardGraph me currentFather " << currentFather << endl;
						}
						else{

							// You have to check if it was deleted first
							// TODO
							Node child=it->b;
							std::vector<MyPair>:: iterator itD;
							for(itD=Backward_del[currentFather].begin(); itD!=Backward_del[currentFather].end(); itD++){
								if(itD->b==child){
									break;
								}
							}
							if(itD!=Backward_del[currentFather].end() || itD->version<localVersion){
								continue;
							}
							visited[child]=visitedCounter+1;
							dist[child]=bGraphDistance;
							bChildrenCount += BackwardGraph[child].nodes.size();
							bQueue[iterator+(bQueuePointer^1)*MAXN]=child;
							iterator++;

						}
					}
				}
			}

			if(bChildrenCount == 0)
				return -1;
			bCurrentNodes = iterator;
			bQueuePointer = bQueuePointer^1;
		}
	}
}

// void add_edge(Node a, Node b) {
//   ForwardGraph[a].push_back(b);
//   BackwardGraph[b].push_back(a);
// }

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
	Node a,b;
	char c;
	unsigned int queryCount;
	vector<Operation_info> operations;
	Operation_info info;

	//cerr << "3ekinaei" << endl;
	ForwardGraph = (GraphNode*)calloc(MAXN, sizeof(GraphNode));
	BackwardGraph = (GraphNode*)calloc(MAXN, sizeof(GraphNode));
	fQueue = (Node*)calloc(MAXN*2, sizeof(Node));
	bQueue = (Node*)calloc(MAXN*2, sizeof(Node));
	visited = (Node*)calloc(MAXN, sizeof(Node));
	dist = (Node*)calloc(MAXN, sizeof(Node));
	nameChange = (Node*)calloc(((size_t)1<<28), sizeof(Node));
	Forward_del = (vector<MyPair>*)calloc(MAXN, sizeof(vector<MyPair>));
	Forward_add = (vector<MyPair>*)calloc(MAXN, sizeof(vector<MyPair>));
	Backward_del = (vector<MyPair>*)calloc(MAXN, sizeof(vector<MyPair>));
	Backward_add = (vector<MyPair>*)calloc(MAXN, sizeof(vector<MyPair>));

	//cerr << "Teleiwse thn ana8esh" << endl;
	//cout << "Ola kala" << endl;

	// for(unsigned long long int  i=0; i<MAXN; i++){
	// 	Forward_del[i].reserve(4);
	// 	Forward_add[i].reserve(4);
	// 	Backward_del[i].reserve(4);
	// 	Backward_add[i].reserve(4);
	// 	//cout << "Kala" << endl;
	// }

	while(scanf("%d %d",&a,&b) == 2){
		if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
		if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
		add_edge_final(nameChange[a],nameChange[b]);
	}
	//cerr << "Teleiwse to diavasma" << endl;
	preprocess();			// Sorting the biggest Nodes to the front of the lists
	/*preprocess2();

	for(unsigned int i=0; i < (size_t)1<<28; i++)
		if(nameChange[i])
			nameChange[i] = nameReorder[nameChange[i]];*/

	scanf("%c",&c);
	printf("R\n");
	fflush(stdout);

	while(scanf(" %c",&c) == 1) {
		if(c == 'F') {
			// cerr << "Telos batch" << endl;
			fflush(stdout);
			
				for(vector<Operation_info>::iterator it=operations.begin(); it!=operations.end(); it++){
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
				operations.clear();
				versionStart=versionCounter;
			

			continue;
		}

		scanf("%d %d",&a,&b);


		if(c == 'Q'){
			//cerr << "Arxizei me query" << endl;
			if(nameChange[a] == 0 || nameChange[b] == 0){
				if(a != b)
					printf("-1\n");
				else
					printf("0\n");
				continue;
			}
			// cerr << "Ola kala edw" << endl;
			printf("%d\n",shortest_path(nameChange[a],nameChange[b],versionCounter));
			// cerr << "epestrepse kala" << endl;
		}else if(c == 'A'){
			// cerr << "Arxizei me add" << endl;
			if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
			if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
			info.op='A'; info.a=nameChange[a]; info.b=nameChange[b];
			operations.push_back(info);
			//cerr << "Komple to add" << endl;

			multiversion_changing(FORWARD,ADD,nameChange[a],nameChange[b]);
			ForwardGraph[a].addition=true;
		//	cerr << "Edw omws einai komple?" << endl;
			multiversion_changing(BACKWARD,ADD,nameChange[b],nameChange[a]);
			BackwardGraph[b].addition=true;
			
		//	cerr << "Mallon edw exoume provlima" << endl;
			//additions[FORWARD].push_back(nameChange[a]); additions[BACKWARD].push_back(nameChange[b]);
			//add_edge(nameChange[a],nameChange[b]);
			versionCounter++;


		}else if(c == 'D') {
			// cerr << "Arxizei me delete" << endl;
			if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
			if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
			info.op='D'; info.a=nameChange[a]; info.b=nameChange[b];
			operations.push_back(info);

			multiversion_changing(FORWARD,DELETE,nameChange[a],nameChange[b]);
			ForwardGraph[a].deletion=true;
			multiversion_changing(BACKWARD,DELETE,nameChange[b],nameChange[a]);
			BackwardGraph[b].deletion=true;
			//MyPairdeletions[FORWARD].push_back(nameChange[a]); deletions[BACKWARD].push_back(nameChange[b]);
			versionCounter++;
			//delete_edge(nameChange[a],nameChange[b]);
		}
	}
  return 0;
}
