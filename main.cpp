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

// Variables for the Queues
vector<Node> *ForwardGraph, *BackwardGraph;
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
	unsigned int fChildrenCount = ForwardGraph[a].size();
	unsigned int fCurrentNodes = 1;
	unsigned int bChildrenCount = BackwardGraph[b].size();
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
				for(unsigned int j=0; j < ForwardGraph[currentFather].size(); j++){
					currentChild = ForwardGraph[currentFather][j];

					if(visited[currentChild] >= visitedCounter){	// Explored by the other side
						if(visited[currentChild] == visitedCounter+1) 	// Found the minimum distance!
							return dist[currentChild] + dist[currentFather] + 1;
						continue;
					}
					visited[currentChild] = visitedCounter;
					dist[currentChild] = fGraphDistance;
					fChildrenCount += ForwardGraph[currentChild].size();	// Counting the children of the next stage
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
				for(unsigned int j=0; j < BackwardGraph[currentFather].size(); j++){
					currentChild = BackwardGraph[currentFather][j];

					if(visited[currentChild] >= visitedCounter){	// Explored by the other side
						if(visited[currentChild] == visitedCounter) 	// Found the minimum distance!
							return dist[currentChild] + dist[currentFather] + 1;
						continue;
					}
					visited[currentChild] = visitedCounter+1;
					dist[currentChild] = bGraphDistance;
					bChildrenCount += BackwardGraph[currentChild].size();	// Counting the children of the next stage
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

bool sortF(Node a,Node b) { return (ForwardGraph[a].size()>ForwardGraph[b].size()); }
bool sortB(Node a,Node b) { return (BackwardGraph[a].size()>BackwardGraph[b].size()); }

void preprocess() {
  for(unsigned int a=1; a < nameChangeCounter; a++) {
    sort(ForwardGraph[a].begin(),ForwardGraph[a].end(),sortF);
    sort(BackwardGraph[a].begin(),BackwardGraph[a].end(),sortB);
  }
}

int main() {
	ios_base::sync_with_stdio(false);
	cin.tie(nullptr);

	Node a,b;
	char c;

	ForwardGraph = (vector<Node>*)calloc(MAXN, sizeof(vector<Node>));
	BackwardGraph = (vector<Node>*)calloc(MAXN, sizeof(vector<Node>));
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
