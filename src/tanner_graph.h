/**
    LDPC
    tanner_graph.h
    Purpose: Builds a Tanner graph based on the LDPC object 
	and finds the number and length of the shortest cycles in it.

    @author Leonov V.R.
    @version 05.05.18
*/


#include "ldpc_generator.h"


#ifndef TANNER_GRAPH_CORRECT_VERSION_H_INCLUDED
#define TANNER_GRAPH_CORRECT_VERSION_H_INCLUDED

//Simple nondirectional graph.
typedef struct {
	int number_of_vertices;		//Number of vertices in graph.
    int **adjacency_list;		//Adjacency list for graph.
    int *degrees_of_vertices;	//Number of neighboring vertices for each.
} graph;

//Set of colors for coloring the graph.
typedef enum {
    WHITE = 0,		//Not colored.
    RED = 1,		//One color.
    BLUE = -RED		//Other color.
} vertex_color;

//Pair of values for describing paths/cycles of same length in graph.
typedef struct {
	int length;			//Length of paths/cycles.
    int number;			//Number of such paths/cycles.
} length_and_number;

/**
    Builds a Tanner graph based on the LDPC object.

    @param ldpc_object LDPC object.
    @return Tanner graph.
*/
graph get_tanner_graph_from_ldpc(ldpc ldpc_object);

/**
    Find shortest cycles in graph and returns information about them.

    @param graph_object Graph.
    @return Length and number of shortest cycles.
*/
length_and_number find_shortest_cycles_in_graph(graph graph_object);

#endif // TANNER_GRAPH_CORRECT_VERSION_H_INCLUDED
