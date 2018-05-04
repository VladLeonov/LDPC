#include <stdlib.h>

#include "tanner_graph.h"


graph get_tanner_graph_from_ldpc(ldpc ldpc_object) {
	int n = ldpc_object.n;
	int r = ldpc_object.r;
	int number_of_vertices = n + r;
	int *degree_of_vertices = (int*) malloc(number_of_vertices * sizeof(int));
	int **adjacency_list = (int**) malloc(number_of_vertices * sizeof(int*));
	
	indices_of_nonzero_elements C = ldpc_object.C;
	indices_of_nonzero_elements V = ldpc_object.V;
	int i, j, length;
	
	for (i = 0; i < n; i++) {
		length = C.element_length[i];
		degree_of_vertices[i] = length;
		adjacency_list[i] = (int*) malloc(length * sizeof(int));
		
		for (j = 0; j < length; j++) {
			adjacency_list[i][j] = C.element_data[i][j] + n;
		}
	}
	
	for (i = 0; i < r; i++) {
		length = V.element_length[i];
		degree_of_vertices[n + i] = length;
		adjacency_list[n + i] = (int*) malloc(length * sizeof(int));
		
		for (j = 0; j < length; j++) {
			adjacency_list[n + i][j] = V.element_data[i][j];
		}
	}
	
	graph tanner_graph;
	tanner_graph.number_of_vertices = number_of_vertices;
	tanner_graph.degree_of_vertices = degree_of_vertices;
	tanner_graph.adjacency_list = adjacency_list;
	
	return tanner_graph;
}
