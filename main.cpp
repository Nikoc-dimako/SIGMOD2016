/*
 * main.cpp
 *
 *  Created on: Mar 6, 2016
 *      Author: thanasis
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <unistd.h>

#include <algorithm>


using namespace std;

#define MAXN ((size_t)1<<24)
typedef unsigned int Node;

typedef struct GraphNode {
   bool addition;
   bool deletion;
   vector<Node> nodes;
} GraphNode;

// Variables for the Queues
GraphNode *ForwardGraph, *BackwardGraph;
Node *visited, *dist;
Node *fQueue;
Node *bQueue;
unsigned int visitedCounter = 0;

//Variables for the node's nameChange
Node *nameChange;
unsigned int nameChangeCounter = 1;

int shortest_path(Node a, Node b) {
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
					if(ForwardGraph[currentFather].deletion == 0){

						// Just add the nodes to the queue
						// TODO

					}else{

						// You have to check if it was deleted first
						// TODO

					}
				}
			}

			if(fChildrenCount == 0)
				return -1;
			fCurrentNodes = iterator;
			fQueuePointer = fQueuePointer^1;

		}else{
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
					if(BackwardGraph[currentFather].deletion == 0){

						// Just add the nodes to the queue
						// TODO

					}else{

						// You have to check if it was deleted first
						// TODO

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

void add_edge(Node a, Node b) {
	// TODO
	ForwardGraph[a].nodes.push_back(b);
	BackwardGraph[b].nodes.push_back(a);
}

void delete_edge(Node a,Node b) {
	// TODO
	if(!_delete_edge(b,ForwardGraph[a].nodes) || !_delete_edge(a,BackwardGraph[b].nodes) )
		return;
}


// Sorting
bool sortF(Node a,Node b) { return (ForwardGraph[a].nodes.size()>ForwardGraph[b].nodes.size()); }
bool sortB(Node a,Node b) { return (BackwardGraph[a].nodes.size()>BackwardGraph[b].nodes.size()); }

void preprocess() {
  for(unsigned int a=1; a < nameChangeCounter; a++) {
    sort(ForwardGraph[a].nodes.begin(),ForwardGraph[a].nodes.end(),sortF);
    sort(BackwardGraph[a].nodes.begin(),BackwardGraph[a].nodes.end(),sortB);
  }
}

int main() {
	ios_base::sync_with_stdio(false);
	cin.tie(nullptr);

	Node a,b;
	char c;

	ForwardGraph = (GraphNode*)calloc(MAXN, sizeof(GraphNode));
	BackwardGraph = (GraphNode*)calloc(MAXN, sizeof(GraphNode));
	fQueue = (Node*)calloc(MAXN*2, sizeof(Node));
	bQueue = (Node*)calloc(MAXN*2, sizeof(Node));
	visited = (Node*)calloc(MAXN, sizeof(Node));
	dist = (Node*)calloc(MAXN, sizeof(Node));
	nameChange = (Node*)calloc(((size_t)1<<28), sizeof(Node));

	while(cin >> a >> b){
		if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
		if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
		add_edge(nameChange[a],nameChange[b]);
	}

	preprocess();			// Sorting the biggest Nodes to the front of the lists

	cout << "R" << endl << flush;

	cin.clear();
	cin >> c;

	while(cin >> c) {
		if(c == 'F') {
			cout << flush;
			cin.clear();

			// Add the additions/deletions of this batch to the main graph
			// TODO

			continue;
		}

		cin >> a >> b;

		if(c == 'Q'){
			if(nameChange[a] == 0 || nameChange[b] == 0){
				if(a != b)
					cout << "-1\n";
				else
					cout << "0\n";
				continue;
			}
			cout << shortest_path(nameChange[a],nameChange[b]) << "\n";
		}else if(c == 'A'){
			if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
			if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
			add_edge(nameChange[a],nameChange[b]);
		}
		else if(c == 'D') {
			if(!nameChange[a]) nameChange[a] = nameChangeCounter++;
			if(!nameChange[b]) nameChange[b] = nameChangeCounter++;
			delete_edge(nameChange[a],nameChange[b]);
		}
	}
  return 0;
}
