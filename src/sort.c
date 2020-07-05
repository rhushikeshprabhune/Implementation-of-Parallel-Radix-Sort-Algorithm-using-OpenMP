#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <omp.h>
#include <math.h> 

#include "sort.h"
#include "edgelist.h"


// Order edges by id of a source vertex, 
// using the Counting Sort
// Complexity: O(E + V)
void countSortEdgesBySource(struct Edge *edges_sorted, struct Edge *edges, int numVertices, int numEdges) {

    
    int i;
    int key;
    int pos;

    // auxiliary arrays, allocated at the start up of the program
    int *vertex_cnt = (int*)malloc(numVertices*sizeof(int)); // needed for Counting Sort

    for(i = 0; i < numVertices; ++i) {
        vertex_cnt[i] = 0;
    }

    // count occurrence of key: id of a source vertex
    for(i = 0; i < numEdges; ++i) {
        key = edges[i].src;
        vertex_cnt[key]++;
    }

    // transform to cumulative sum
    for(i = 1; i < numVertices; ++i) {
        vertex_cnt[i] += vertex_cnt[i - 1];
    }

    // fill-in the sorted array of edges
    for(i = numEdges - 1; i >= 0; --i) {
        key = edges[i].src;
        pos = vertex_cnt[key] - 1;
        edges_sorted[pos] = edges[i];
        vertex_cnt[key]--;
    }


    free(vertex_cnt);

}

//Serial Radix Sort
 void radixSortEdgesBySourceSerial(struct Edge *edges_sorted, struct Edge *edges, int numVertices, int numEdges) {
 int i;
    int key;
    int pos;
    int divisor=1;

    // auxiliary arrays, allocated at the start up of the program
    int *vertex_cnt = (int*)malloc(numVertices*sizeof(int)); // needed for Counting Sort

    int max;
    max = edges[0].src;
    for(i=1;i<numEdges;i++){
       if(edges[i].src > max)
       max=edges[i].src;
    }
       int l = 0;
       int n = max;
       while(n!=0){
        n /= 10;
        ++l;
        }

        for(int j=0;j<l;j++){

    for(i = 0; i < numVertices; ++i) {
        vertex_cnt[i] = 0;
    }

    // count occurrence of key: id of a source vertex
    for(i = 0; i < numEdges; ++i) {
        key = (edges[i].src/divisor)%10;
        vertex_cnt[key]++;
    }


    // transform to cumulative sum
    for(i = 1; i < numVertices; ++i) {
        vertex_cnt[i] += vertex_cnt[i - 1];
    }


    // fill-in the sorted array of edges
    for(i = numEdges - 1; i >= 0; --i) {
        key = (edges[i].src/divisor)%10;
        pos = vertex_cnt[key] - 1;
        edges_sorted[pos] = edges[i];
        vertex_cnt[key]--;
    }

    for(i=0;i<numEdges;i++){
        edges[i]= edges_sorted[i];
         }
        divisor *= 10;
 
 }
 free(vertex_cnt);
 }
#define NUM_THREADS 16
 void radixSortEdgesBySourceParallel(struct Edge *edges_sorted, struct Edge *edges, int numVertices, int numEdges) {
  int i, key;
  int j, l;
  int index;
  int pos;

  omp_set_num_threads(NUM_THREADS); // Request threads
  int Nthreads;
  int max;
  max = edges[0].src;
  for (i = 1; i < numEdges; ++i) {
    key = edges[i].src;
	if (key > max){
	  max = key;
	}
  }
  l = (log10(max)) + 1;

  for (j = 0; j < l; j++){
    int **vertex_cnt = malloc(NUM_THREADS*sizeof(int*));
	for (i=0; i<NUM_THREADS;i++){
    	vertex_cnt[i] = malloc(numVertices*sizeof(int));
	}

	#pragma omp parallel
	{
		int id, i, Nthrds, istart, iend;
		int key, index;
		id = omp_get_thread_num();
		Nthrds = omp_get_num_threads();

		for(i = 0; i < numVertices; ++i) {
		  vertex_cnt[id][i] = 0;
		}
	
		istart = id * numEdges / Nthrds;
		iend = (id+1) * numEdges / Nthrds;
		if (id == Nthrds-1) {iend = numEdges; Nthreads = Nthrds;}
		// count occurrence of key: id of a source vertex
		for(i = istart; i < iend; ++i) {
			key = edges[i].src;
			index = (key/(int)pow(10,j))%10;
			vertex_cnt[id][index]++;
		}
	}


	// transform to cumulative sum
    int id;
    int base = 0;
	for(i = 0; i < numVertices; ++i) {
		for (id = 0; id < Nthreads; ++id){
        	base = base + vertex_cnt[id][i];
            vertex_cnt[id][i] = base;
		}
    }
	#pragma omp parallel
	{
		int id, i, Nthrds, istart, iend;
		int pos;
		int key, index;

		id = omp_get_thread_num();
		Nthrds = omp_get_num_threads();

		istart = id * numEdges / Nthrds;
		iend = (id+1) * numEdges / Nthrds;

		// fill-in the sorted array of edges
		for(i = iend-1; i >= istart; i--) {
			key = edges[i].src;
			index = (key/(int)pow(10,j))%10;
			pos = vertex_cnt[id][index] - 1;
			edges_sorted[pos] = edges[i];
			vertex_cnt[id][index] = vertex_cnt[id][index] - 1;
		}
	}

 	for(i = 0; i < numEdges; i++){
		edges[i] = edges_sorted[i];
 	}

	free(vertex_cnt);
  } 

 }


