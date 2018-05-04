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


void remove_edge_from_graph(graph graph_object, int vertex_index) {
	if (graph_object.degree_of_vertices[vertex_index] == 0) {
		return;
	}
	
	int other_vertex_index = graph_object.adjacency_list[vertex_index][0];
	
	int i;
	for (i = 0; i < graph_object.degree_of_vertices[vertex_index] - 1; i++) {
		graph_object.adjacency_list[vertex_index][i] = graph_object.adjacency_list[vertex_index][i + 1];
	}
	
	graph_object.degree_of_vertices[vertex_index]--;
	graph_object.adjacency_list[vertex_index] = 
		(int*) realloc(graph_object.adjacency_list[vertex_index], 
		graph_object.degree_of_vertices[vertex_index] * sizeof(int));
		
	for (i = 0; graph_object.adjacency_list[other_vertex_index][i] != vertex_index; i++);
	
	for (; i < graph_object.degree_of_vertices[other_vertex_index] - 1; i++) {
		graph_object.adjacency_list[other_vertex_index][i] = graph_object.adjacency_list[other_vertex_index][i + 1];
	}
	
	graph_object.degree_of_vertices[other_vertex_index]--;
	graph_object.adjacency_list[other_vertex_index] = 
		(int*) realloc(graph_object.adjacency_list[other_vertex_index], 
		graph_object.degree_of_vertices[other_vertex_index] * sizeof(int));
}
