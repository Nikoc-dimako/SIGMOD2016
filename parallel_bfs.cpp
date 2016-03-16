/*
 * main.cpp
 *
 *  Created on: Mar 6, 2016
 *      Author: thanasis
1 2
2 3
3 1
4 1
2 4
S
Q 1 3
A 4 5
Q 1 5
Q 5 1
F
A 5 3
Q 1 3
D 2 3
Q 1 3
F

 */


#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <algorithm>


using namespace std;

#define MAXN ((size_t)1<<22)
#define MAXTH 7
#define READ 0
#define WRITE 1
#define perror2(s,e) fprintf(stderr , "%s: %s\n", s, strerror(e))
typedef unsigned int Node;

vector<Node> *ForwardG, *BackwardG;
Node *Queue[MAXTH], *Distance[MAXTH], *Visited[MAXTH];

class Operations_graph{
  public: 
    char operation;
    Node a,b;
    unsigned int version;
};

class Graph_info{
	public:
		Node a,b;
		vector<Node> *ForwardG, *BackwardG;
		int query;
};

class Thread_info{
	public:
		size_t id;
		int fd[2];
};
vector<Thread_info> list_th;
vector<int> list_res;
int tsk=0;
pthread_mutex_t tsk_mtx=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t fg_bg_mtx=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t res_mtx=PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t thr_ls_mtx=PTHREAD_MUTEX_INITIALIZER;
int err;

int shortest_path(size_t thread,Node a, Node currentNode,vector<Node> *ForwardGraph,vector<Node> *BackwardGraph) {
	// Current Node is used for destination, but also for iterating in order to not create a new one

 // printf("Edw ftanw\n");
  unsigned int VisitedCounter = 0;

	if(a == currentNode)
		return 0;

	VisitedCounter += 2;
	int startOfForward = 0, endOfForward = 0, itForward = 0;
	int startOfBackward = 1, endOfBackward = 1, itBackward = 0;
	Node currentFNode = a, currentBNode = currentNode; /*from = a, to = currentNode;*/
  // Node Queue[MAXN],Distance[MAXN],Visited[MAXTH];
 //  printf("edw edw\n");
	// Node *Queue=(Node *)Queue[thread],*Visited=(Node *)Visited[thread],*Distance=(Node *)Distance[thread];
 //  printf("ekei ekei\n");

  // Queue = (Node*)calloc(MAXN, sizeof(Node));
  // Distance = (Node*)calloc(MAXN, sizeof(Node));
  // Visited = (Node*)calloc(MAXN, sizeof(Node));

	// Initialize Scan Variables
 // printf("Ola kala\n");
	Queue[thread][endOfForward+=2] = a;
	Visited[thread][a] = VisitedCounter;
	Distance[thread][a] = 0;
	Queue[thread][endOfBackward+=2] = currentNode;
	Visited[thread][currentNode] = VisitedCounter+1;
	Distance[thread][currentNode] = 0;
	int minimumDistance = -1;

  // printf("Mexri edw komple\n");

	while(1) { // Main loop
		if(minimumDistance >= 0 && minimumDistance <= Distance[thread][currentFNode]+Distance[thread][currentBNode]+1)
				break;
		if(endOfForward <= endOfBackward) {
			while(!itForward && startOfForward<endOfForward) {
				currentFNode = Queue[thread][startOfForward+=2];
				itForward = ForwardGraph[currentFNode].size();
			}
			if(!itForward)
				break;
			currentNode = ForwardGraph[currentFNode][--itForward];
			if(Visited[thread][currentNode] >= VisitedCounter) {
				if(Visited[thread][currentNode] == VisitedCounter + 1) {
					if(minimumDistance >= 0 && minimumDistance < Distance[thread][currentFNode]+Distance[thread][currentNode]+1)
						break;
					minimumDistance = Distance[thread][currentFNode]+Distance[thread][currentNode]+1;
				}
				continue;
			}
			Queue[thread][endOfForward+=2] = currentNode;
			Visited[thread][currentNode] = VisitedCounter;
			Distance[thread][currentNode] = Distance[thread][currentFNode]+1;
		}else{
			while(!itBackward && startOfBackward<endOfBackward) {
				currentBNode = Queue[thread][startOfBackward+=2];
				itBackward = BackwardGraph[currentBNode].size();
			}
			if(!itBackward)
				break;
			currentNode = BackwardGraph[currentBNode][--itBackward];
			if(Visited[thread][currentNode] >= VisitedCounter) {
				if(Visited[thread][currentNode] == VisitedCounter) {
					if(minimumDistance >= 0 && minimumDistance < Distance[thread][currentBNode]+Distance[thread][currentNode]+1)
					  break;
					minimumDistance = Distance[thread][currentBNode]+Distance[thread][currentNode]+1;
				}
				continue;
			}
			Queue[thread][endOfBackward+=2] = currentNode;
			Visited[thread][currentNode] = VisitedCounter+1;
			Distance[thread][currentNode] = Distance[thread][currentBNode]+1;
		}
	}
  // free(Queue);
  // free(Visited);
  // free(Distance);
	return minimumDistance;
}

void add_edge(Node a, Node b) {
  ForwardG[a].push_back(b);
  BackwardG[b].push_back(a);
}
int _delete_edge(Node a, vector<Node> &E) {
  int sz = E.size(), nsz = 0;
  for(int i=0; i < sz; i++){
    E[nsz] = E[i];
    if(E[i] != a)
		nsz++;
  }
  if(nsz == sz)
	  return 0;
  E.resize(nsz);
  return 1;
}
void delete_edge(Node a,Node b) {
  if(!_delete_edge(b,ForwardG[a]) || !_delete_edge(a,BackwardG[b]) )
	  return;
}

void *thread_computation(void *parametrs){

	Graph_info graph_path;
	Thread_info *thr,thr_inf;
	thr = (Thread_info*) parametrs;
  thr_inf = *thr;
	// close(thr_inf.fd[WRITE]);
	while(1){
		int k=0;
		int temp=0;
		while(temp!=sizeof(Graph_info)){
			k=read(thr_inf.fd[READ],&graph_path,sizeof(Graph_info)-temp);
			if(k<0){
				perror("Read thread");
				exit(-1);
			}
			temp+=k;
		}
    // printf("komple b: %d\n",graph_path.b);
    // printf("Ola komple apo thread: %ld\n",pthread_self());
		int res=shortest_path(thr_inf.id,graph_path.a,graph_path.b,graph_path.ForwardG,graph_path.BackwardG);
    // printf("Ola komple2 apo thread: %ld distance are %d \n",pthread_self(),res);





		if(err=pthread_mutex_lock(&res_mtx)){
	        perror2("pthread_mutex_lock",err);
	        exit(-1);
    }
		list_res.insert(list_res.begin()+graph_path.query,res);
		if(err=pthread_mutex_unlock(&res_mtx)){
	        perror2("pthread_mutex_unlock",err);
	        exit(-1);
    }

		if(err=pthread_mutex_lock(&thr_ls_mtx)){
        	perror2("pthread_mutex_lock",err);
        	exit(-1);
    	}
		list_th.push_back(thr_inf);
		if(err=pthread_mutex_unlock(&thr_ls_mtx)){
        	perror2("pthread_mutex_unlock",err);
        	exit(-1);
    	}
	}

}

int main() {
  Node a,b;
  char c;
  vector<Node> *ForwardGraph[MAXTH], *BackwardGraph[MAXTH];
  int fd[MAXTH][2],queryCount=0;
  pthread_t threads[MAXTH];
  Thread_info thr_temp[MAXTH],thr;
  Graph_info graph;


  ForwardG = (vector<Node>*)calloc(MAXN, sizeof(vector<Node>));
  BackwardG = (vector<Node>*)calloc(MAXN, sizeof(vector<Node>));
  // Queue = (Node*)calloc(MAXN, sizeof(Node));
  // Distance = (Node*)calloc(MAXN, sizeof(Node));
  // Visited = (Node*)calloc(MAXN, sizeof(Node));
  for(int i=0; i<MAXTH; i++){
    Queue[i]=(Node *)calloc(MAXN, sizeof(vector<Node>));
    Distance[i]=(Node *)calloc(MAXN, sizeof(vector<Node>));
    Visited[i]=(Node *)calloc(MAXN, sizeof(vector<Node>));
  }
  while(scanf("%d %d",&a,&b) == 2)
	  add_edge(a,b);
  scanf(" %c",&c);
  for(int i=0; i<MAXTH; i++){
  	pipe(fd[i]);
  	thr_temp[i].id=i;
  	thr_temp[i].fd[READ]=fd[i][READ];
    thr_temp[i].fd[WRITE]=fd[i][WRITE];
  	pthread_create(&threads[i],NULL,thread_computation,(void *) &thr_temp[i]);
  	list_th.push_back(thr_temp[i]);
  }

  printf("R\n");
  fflush(stdout);

  while(scanf(" %c",&c) == 1) {
    if(c == 'F') {
      while(1){
      	if(err=pthread_mutex_lock(&res_mtx)){
        	perror2("pthread_mutex_lock",err);
        	exit(-1);
    	  }
      	if(list_res.size()==queryCount){
      		for(std::vector<int>::iterator it = list_res.begin(); it != list_res.end(); ++it){
      			printf("%d\n",*it); fflush(stdout);
      		}
      		list_res.clear();
      		if(err=pthread_mutex_unlock(&res_mtx)){
     		   perror2("pthread_mutex_unlock",err);
        		exit(-1);
    		}
      		break;
      	}
      	if(err=pthread_mutex_unlock(&res_mtx)){
     		   perror2("pthread_mutex_unlock",err);
        		exit(-1);
    	  }
      }
      queryCount=0;
      
      continue;
    }
    scanf("%d %d",&a,&b);
    if(c=='Q'){
    	queryCount++;
    	while(1){
    		if(err=pthread_mutex_lock(&thr_ls_mtx)){
     		   perror2("pthread_mutex_lock",err);
        		exit(-1);
    		}
    		if(list_th.size()){
    			thr=*list_th.begin();
    			list_th.erase(list_th.begin());
    			if(err=pthread_mutex_unlock(&thr_ls_mtx)){
	     		   perror2("pthread_mutex_unlock",err);
	        		exit(-1);
	    		}
    			break;
    		}
    		if(err=pthread_mutex_unlock(&thr_ls_mtx)){
     		   perror2("pthread_mutex_unlock",err);
        		exit(-1);
    		}
    	}


     
      //  memset(BackwardGraph[thr.id],0,MAXN*sizeof(std::vector<Node>));
      // memset(ForwardGraph[thr.id],0,MAXN*sizeof(std::vector<Node>));
    	// memcpy(ForwardGraph[thr.id],ForwardG,MAXN*sizeof(std::vector<Node>));
     //  memcpy(BackwardGraph[thr.id],BackwardG,MAXN*sizeof(std::vector<Node>));
      // copy(ForwardG,ForwardG+MAXN,ForwardGraph[thr.id]);
      // copy(BackwardG,BackwardG+MAXN,BackwardGraph[thr.id]);
      // printf("main thread \n");
      // printf("Forward 1: %d\n",ForwardGraph[thr.id][1][0]);

    	graph.a=a; graph.b=b; graph.ForwardG=ForwardG; graph.BackwardG=BackwardG; graph.query=queryCount-1;

    	int k=0,temp=0;
    	while(temp!=sizeof(Graph_info)){
    		k=write(thr.fd[WRITE],&graph,sizeof(Graph_info)-temp);
    		if(k<0){
    			perror("Write graph");
    			exit(-1);
    		}
    		temp += k;
    	}

    }
    if(c=='A' || c=='D'){ //printf("%d\n",shortest_path(a,b));

      while(1){
        if(err=pthread_mutex_lock(&thr_ls_mtx)){
             perror2("pthread_mutex_lock",err);
              exit(-1);
        }
        if(list_th.size()==MAXTH){
          if(err=pthread_mutex_unlock(&thr_ls_mtx)){
           perror2("pthread_mutex_unlock",err);
            exit(-1);
          }
          break;
        }
        if(err=pthread_mutex_unlock(&thr_ls_mtx)){
           perror2("pthread_mutex_unlock",err);
            exit(-1);
        }
      }


      if(c=='A') add_edge(a,b);
      if(c=='D') delete_edge(a,b);

    }
  }
  return 0;
}


