#include "ldpc_generator.h"

#ifndef TANNER_GRAPH_CORRECT_VERSION_H_INCLUDED
#define TANNER_GRAPH_CORRECT_VERSION_H_INCLUDED

typedef struct {
	int number_of_vertices;
    int **adjacency_list;
    int *degree_of_vertices;
} graph;

typedef enum {
    WHITE = 0,
    RED = 1,
    BLUE = -RED
} vertex_color;

graph get_tanner_graph_from_ldpc(ldpc ldpc_object);
void remove_edge_from_graph(graph graph_object, int vertex_index);

#endif // TANNER_GRAPH_CORRECT_VERSION_H_INCLUDED
