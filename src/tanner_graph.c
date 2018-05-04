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


length_and_number find_shortest_paths_between_vertices(graph graph_object, 
													   int vertex_index1,
													   int vertex_index2) {
	length_and_number shortest_paths = {1, 1}; //TODO don't be a stub
	return shortest_paths;
}


length_and_number find_shortest_cycles_in_graph(graph graph_object) {
	length_and_number shortest_cycles = {INT_MAX, 0};
	length_and_number shortest_paths;
	int vertex_index1, vertex_index2;
	
	for (vertex_index1 = 0; graph_object.degree_of_vertices[vertex_index1] > 0; vertex_index1++) {
		
		vertex_index2 = graph_object.adjacency_list[vertex_index1][0];
		remove_edge_from_graph(graph_object, vertex_index1);
		shortest_paths = find_shortest_paths_between_vertices(graph_object, vertex_index1, vertex_index2);
		
		if (shortest_paths.length < shortest_cycles.length) {
			shortest_cycles = shortest_paths;
		} else if (shortest_paths.length == shortest_cycles.length) {
			shortest_cycles.number += shortest_paths.number;
		}
	}
	
	shortest_cycles.length++;
	return shortest_cycles;
}
