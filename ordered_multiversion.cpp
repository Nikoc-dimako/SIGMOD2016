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

// Variables for the Queues
vector<Node> *ForwardGraph, *BackwardGraph;
Node *visited, *dist;
Node *fQueue;
Node *bQueue;
int fd[MAXTH][2];
unsigned int visitedCounter = 0;


struct Operation_info{
	char op;
	Node a,b;
};

struct Graph_info{
	//int id; 8a to steilw san parametro
	unsigned int query,version;
	Node a,b;

};

struct MyPair{
	char op;
	Node b;
};



//Variables for the node's nameChange
Node *nameChange;
Node *nameReorder;
unsigned int nameChangeCounter = 1;

size_t versionStart=0,versionCounter=0;
map<Node,map<size_t,MyPair>> multiversion[2];

vector<Node> construct_the_right_children(int graph,Node a,size_t version){
	std::vector<Node> children;
	map<size_t,MyPair> version_of_operations;
	map<Node,map<size_t,MyPair>>::iterator itAD;

	if(graph==FORWARD){
		children=ForwardGraph[a];
	}
	else{
		children=BackwardGraph[a];
	}

	itAD=multiversion[graph].find(a);

	if(itAD!=multiversion[graph].end()){
		version_of_operations=itAD->second;

		map<size_t,MyPair>::iterator it_add_del;

		it_add_del=version_of_operations.begin();

		while(it_add_del!=version_of_operations.end()){
			if(it_add_del->first<version){
				if(it_add_del->second.op=='A'){
					children.push_back(it_add_del->second.b);
				}
				else{
					std::vector<Node>::iterator temp;
					temp=find(children.begin(),children.end(),it_add_del->second.b);
					if(temp!=children.end()){
						children.erase(temp);
					}
				}
				it_add_del++;
			}
			else{
				break;
			}
		}
	}


	return children;

}

int shortest_path(Node a, Node b,size_t localVersion) {
	if(a == b)
		return 0;

	visitedCounter += 2;

	vector<Node> children_F=construct_the_right_children(FORWARD,a,localVersion);
	std::vector<Node> children_B=construct_the_right_children(BACKWARD,b,localVersion);


	// Initialize Queues and Variables
	unsigned int fChildrenCount = children_F.size();
	unsigned int fCurrentNodes = 1;
	unsigned int bChildrenCount = children_B.size();
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
	unsigned int currentFather;
	unsigned int currentChild;

	// The queue saves the currently explored nodes
	// The children counters contain the number of children that the already explored nodes have.

	// Main Loop
	while(1){
		iterator = 0;

		if(fChildrenCount < bChildrenCount){	// Move forward, there are less children there
			fChildrenCount = 0;
			fGraphDistance++;	// Going to the next distance
			for(unsigned int i=0; i < fCurrentNodes; i++){
				currentFather = fQueue[i+fQueuePointer*MAXN];
				children_F=construct_the_right_children(FORWARD,currentFather,localVersion);
				for(unsigned int j=0; j < children_F.size(); j++){
					currentChild = children_F[j];

					if(visited[currentChild] >= visitedCounter){	// Explored by the other side
						if(visited[currentChild] == visitedCounter+1) 	// Found the minimum distance!
							return dist[currentChild] + dist[currentFather] + 1;
						continue;
					}
					visited[currentChild] = visitedCounter;
					dist[currentChild] = fGraphDistance;
					std::vector<Node> tempChild_F=construct_the_right_children(FORWARD,currentChild,localVersion);
					fChildrenCount += tempChild_F.size();	// Counting the children of the next stage
					fQueue[iterator+(fQueuePointer^1)*MAXN] = currentChild;
					iterator++;
				}
			}
			if(fChildrenCount == 0)
				return -1;
			fCurrentNodes = iterator;
			fQueuePointer = fQueuePointer^1;

		}else{
			bChildrenCount = 0;
			bGraphDistance++;	// Going to the next distance
			for(unsigned int i=0; i < bCurrentNodes; i++){
				currentFather = bQueue[i+bQueuePointer*MAXN];
				children_B=construct_the_right_children(BACKWARD,currentFather,localVersion);
				for(unsigned int j=0; j < children_B.size(); j++){
					currentChild = children_B[j];

					if(visited[currentChild] >= visitedCounter){	// Explored by the other side
						if(visited[currentChild] == visitedCounter) 	// Found the minimum distance!
							return dist[currentChild] + dist[currentFather] + 1;
						continue;
					}
					visited[currentChild] = visitedCounter+1;
					dist[currentChild] = bGraphDistance;
					std::vector<Node> tempChild_B=construct_the_right_children(BACKWARD,currentChild,localVersion);
					bChildrenCount += tempChild_B.size();	// Counting the children of the next stage
					bQueue[iterator+(bQueuePointer^1)*MAXN] = currentChild;
					iterator++;
				}
			}
			if(bChildrenCount == 0)
				return -1;
			bCurrentNodes = iterator;
			bQueuePointer = bQueuePointer^1;
		}
	}
}

void add_edge(Node a, Node b) {
  ForwardGraph[a].push_back(b);
  BackwardGraph[b].push_back(a);
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
  if(!_delete_edge(b,ForwardGraph[a]) || !_delete_edge(a,BackwardGraph[b]) )
	  return;
}

bool sortF(Node a,Node b) { return (ForwardGraph[a].size() > ForwardGraph[b].size()); }
bool sortB(Node a,Node b) { return (BackwardGraph[a].size() > BackwardGraph[b].size()); }

void preprocess() {
  for(unsigned int a=1; a < nameChangeCounter; a++) {
    sort(ForwardGraph[a].begin(),ForwardGraph[a].end(),sortF);
    sort(BackwardGraph[a].begin(),BackwardGraph[a].end(),sortB);
  }
}

// Order the list of vectors, the one with most children first
void preprocess2(){
	int maxChildren = 0;
	unsigned int maxChildrenPos;
	unsigned int temp, tempx, tempy;
	nameReorder = (Node*)calloc(nameChangeCounter, sizeof(Node));
	vector<Node> tempNode;

	// Create the swap board
	for(unsigned int i=1; i < nameChangeCounter; i++)
		nameReorder[i] = i;

	// Start sortings
	for(unsigned int i=1; i < nameChangeCounter; i++){
		maxChildren = -1;
		for(unsigned int j=i; j < nameChangeCounter; j++)
			if((int)ForwardGraph[j].size() > maxChildren){
				maxChildren = (int)ForwardGraph[j].size();
				maxChildrenPos = j;
			}

		// Change the position of the max Vector
		tempNode = ForwardGraph[i];
		ForwardGraph[i] = ForwardGraph[maxChildrenPos];
		ForwardGraph[maxChildrenPos] = tempNode;

		// Change the position of the Backward Vectors
		tempNode = BackwardGraph[i];
		BackwardGraph[i] = BackwardGraph[maxChildrenPos];
		BackwardGraph[maxChildrenPos] = tempNode;

		// Swap Board holds the swaps that are made to know which vector Node went where
		// Swap nameswap[nameswap[i]] <-> nameswap[nameswap[maxChildrenPos]]
		tempx = nameReorder[i];
		tempy = nameReorder[maxChildrenPos];
		temp = nameReorder[tempx];
		nameReorder[tempx] = nameReorder[tempy];
		nameReorder[tempy] = temp;

		if(maxChildren < 200)
			break;
	}

	// We have to start swapping the inside of the vectors
	/* Forward Graph */
	for(unsigned int i=1; i < nameChangeCounter; i++)
		for(unsigned int j=0; j < ForwardGraph[i].size(); j++)
			ForwardGraph[i][j] = nameReorder[ForwardGraph[i][j]];

	/* Backward Graph */
	for(unsigned int i=1; i < nameChangeCounter; i++)
		for(unsigned int j=0; j < BackwardGraph[i].size(); j++)
			BackwardGraph[i][j] = nameReorder[BackwardGraph[i][j]];
}

void multiversion_changing(int graph,int operation,Node a,Node b){
	map<Node,map<size_t,MyPair>>::iterator it;
	map<size_t,MyPair> newVersion;
	MyPair pair;

	if(operation){
		pair.op='D';
	}
	else{
		pair.op='A';
	}

	pair.b=b;
	it=multiversion[graph].find(a);

	if(it!=multiversion[graph].end()){
		newVersion=it->second;
		newVersion[versionCounter]=pair;
		multiversion[graph][a]=newVersion;
	}
	else{
		newVersion[versionCounter]=pair;
		multiversion[graph][a]=newVersion;
	}
}

int main() {
	Node a,b;
	char c;
	unsigned int queryCount;
	vector<Operation_info> operations;
	Operation_info info;

	ForwardGraph = (vector<Node>*)calloc(MAXN, sizeof(vector<Node>));
	BackwardGraph = (vector<Node>*)calloc(MAXN, sizeof(vector<Node>));
	fQueue = (Node*)calloc(MAXN*2, sizeof(Node));
	bQueue = (Node*)calloc(MAXN*2, sizeof(Node));
	visited = (Node*)calloc(MAXN, sizeof(Node));
	dist = (Node*)calloc(MAXN, sizeof(Node));
	nameChange = (Node*)calloc(((size_t)1<<28), sizeof(Node));

	while(scanf("%d %d",&a,&b) == 2){
		if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
		if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
		add_edge(nameChange[a],nameChange[b]);
	}

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
			fflush(stdout);
			
				for(vector<Operation_info>::iterator it=operations.begin(); it!=operations.end(); it++){
					if(it->op=='A'){
						add_edge(it->a,it->b);
					}
					else{
						delete_edge(it->a,it->b);
					}
				}
				multiversion[FORWARD].clear();
				multiversion[BACKWARD].clear();
				operations.clear();

				versionStart=versionCounter;
			

			continue;
		}

		scanf("%d %d",&a,&b);


		if(c == 'Q'){
			if(nameChange[a] == 0 || nameChange[b] == 0){
				if(a != b)
					printf("-1\n");
				else
					printf("0\n");
				continue;
			}
			printf("%d\n",shortest_path(nameChange[a],nameChange[b],versionCounter));
		}else if(c == 'A'){
			if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
			if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
			info.op='A'; info.a=nameChange[a]; info.b=nameChange[b];
			operations.push_back(info);

			multiversion_changing(FORWARD,ADD,nameChange[a],nameChange[b]);
			multiversion_changing(BACKWARD,ADD,nameChange[b],nameChange[a]);
			//add_edge(nameChange[a],nameChange[b]);
			versionCounter++;

		}else if(c == 'D') {
			if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
			if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
			info.op='D'; info.a=nameChange[a]; info.b=nameChange[b];
			operations.push_back(info);

			multiversion_changing(FORWARD,DELETE,nameChange[a],nameChange[b]);
			multiversion_changing(BACKWARD,DELETE,nameChange[b],nameChange[a]);
			versionCounter++;
			//delete_edge(nameChange[a],nameChange[b]);
		}
	}
  return 0;
}
