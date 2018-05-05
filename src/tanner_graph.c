#include <stdlib.h>

#include "tanner_graph.h"


#define TRUE !0
#define FALSE 0


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


int perform_iteration_of_graph_colorization(graph graph_object, int* last_vertices, int* num_last_vertices, vertex_color color, vertex_color* colors) {
	int last_vertices_buffer[graph_object.number_of_vertices];
	int num_last_vertices_buffer = 0;
	int current_vertex, neighbor_vertex;
	int number_of_paths = 0;
	int i, j;
	
	for (i = 0; i < *num_last_vertices; i++) {
		current_vertex = last_vertices[i];
		
		for (j = 0; j < graph_object.degree_of_vertices[current_vertex]; j++) {
			neighbor_vertex = graph_object.adjacency_list[current_vertex][j];
			
			if (colors[neighbor_vertex] == WHITE) {
				colors[neighbor_vertex] = color;
				last_vertices_buffer[num_last_vertices_buffer++] = neighbor_vertex;
			} else if (colors[neighbor_vertex] == -color) {
				number_of_paths++;
			}
		}
	}
	
	for (i = 0; i < num_last_vertices_buffer; i++) {
		last_vertices[i] = last_vertices_buffer[i];
	}
	*num_last_vertices = num_last_vertices_buffer;
	
	return number_of_paths;
}


length_and_number find_shortest_paths_between_vertices(graph graph_object, 
													   int vertex_index1,
													   int vertex_index2) {
	int V = graph_object.number_of_vertices;
	vertex_color colors[V];
	
	int i, j;
	for (i = 0; i < V; i++) {
		colors[i] = WHITE;
	}
	colors[vertex_index1] = RED;
	colors[vertex_index2] = BLUE;
	
	int last_red_vertices[V], last_blue_vertices[V];
	last_red_vertices[0] = vertex_index1;
	last_blue_vertices[0] = vertex_index2;
	int num_last_red_vertices = 1, num_last_blue_vertices = 1;
	
	length_and_number shortest_paths = {INT_MAX, 0};
	int path_length = 0;
	int number_of_paths;
	
	while (num_last_red_vertices > 0) {
		path_length++;
		number_of_paths = perform_iteration_of_graph_colorization(graph_object, last_red_vertices, &num_last_red_vertices, RED, colors);
		if (number_of_paths != 0) {
			shortest_paths.length = path_length;
			shortest_paths.number += number_of_paths;
			break;
		}
		
		path_length++;
		number_of_paths = perform_iteration_of_graph_colorization(graph_object, last_blue_vertices, &num_last_blue_vertices, BLUE, colors);
		if (number_of_paths != 0) {
			shortest_paths.length = path_length;
			shortest_paths.number += number_of_paths;
			break;
		}
	}
	
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
