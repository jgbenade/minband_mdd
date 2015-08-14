/**
 * -------------------------------------------------
 * Independent Set structure - Implementation
 * -------------------------------------------------
 */

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "minband_instance.hpp"
#include "util.hpp"


using namespace std;


//
// Read graph in DIMACS format
//
void Graph::read_dimacs(const char* filename) {

	string buffer;
	char   command;

	ifstream input(filename);

	if (!input.is_open()) {
		cerr << "Error: could not open DIMACS graph file " << filename << endl << endl;
		exit(1);
	}

	int read_edges = 0;
	n_edges = -1;

	int source, target;

	while ( read_edges != n_edges && !input.eof() ) {

		input >> command;

		if (command == 'c') {
			// read comment
			getline(input, buffer);

		} else if (command == 'n') {
			// read weight
			input >> source;
			source--;
			input >> weights[source];

		} else if (command == 'p') {
			// read 'edge' or 'col'
			input >> buffer;

			// read number of vertices and edges
			input >> n_vertices;
			input >> n_edges;

			// allocate adjacent matrix
			adj_m = new bool*[n_vertices];
			for (int i = 0; i < n_vertices; i++) {
				adj_m[i] = new bool[n_vertices];
				memset(adj_m[i], false, sizeof(bool)*n_vertices);
			}

			// allocate adjacent list
			adj_list.resize(n_vertices);

			// initialize weights
			weights = new double[n_vertices];
			for (int i = 0; i < n_vertices; ++i) {
				weights[i] = 1.0;
			}

		} else if (command == 'e') {

			if (input.eof()) {
				break;
			}

			// read edge
			input >> source;
			source--;

			input >> target;
			target--;

			set_adj(source, target);
			set_adj(target, source);

			read_edges++;
		}

	}

	input.close();

	//    // we assume a vertex is adjacent to itself
	//    for( int i = 0; i < n_vertices; i++ ) {
	//        set_adj(i, i);
	//    }

	int count_edges = 0;
	for (int i = 0; i < n_vertices; i++) {
		for (int j = i+1; j < n_vertices; j++) {
			if (is_adj(i,j)) {
				count_edges++;
			}
		}
	}

	n_edges = count_edges;

	for (int i=0; i<n_vertices; i++){
		vertex_neighbourhood.push_back(make_layered_graph(i));
	}

	//print();
	calculate_pairwise_dist();


}

void Graph::calculate_pairwise_dist(){
	/* dist[][] will be the output matrix that will finally have the shortest
	      distances between every pair of vertices */
	int dist[n_vertices][n_vertices], i, j, k;

	/* Initialize the solution matrix same as input graph matrix. Or
	       we can say the initial values of shortest distances are based
	       on shortest paths considering no intermediate vertex. */
	for (i = 0; i < n_vertices; i++)
		for (j = 0; j < n_vertices; j++){
			if (adj_m[i][j] > 0 or i==j)
				dist[i][j] = adj_m[i][j];
			else
				dist[i][j] = n_vertices*n_vertices;
		}

	/* Add all vertices one by one to the set of intermediate vertices.
	      ---> Before start of a iteration, we have shortest distances between all
	      pairs of vertices such that the shortest distances consider only the
	      vertices in set {0, 1, 2, .. k-1} as intermediate vertices.
	      ----> After the end of a iteration, vertex no. k is added to the set of
	      intermediate vertices and the set becomes {0, 1, 2, .. k} */
	for (k = 0; k < n_vertices; k++)
	{
		// Pick all vertices as source one by one
		for (i = 0; i < n_vertices; i++)
		{
			// Pick all vertices as destination for the
			// above picked source
			for (j = 0; j < n_vertices; j++)
			{
				// If vertex k is on the shortest path from
				// i to j, then update the value of dist[i][j]
				if (dist[i][k] + dist[k][j] < dist[i][j])
					dist[i][j] = dist[i][k] + dist[k][j];
			}
		}
	}

	//copy distnace over
	for (int i = 0; i < n_vertices; i++)
	{
		vector<int> d(n_vertices,0);
		for (int j = 0; j < n_vertices; j++)
		{
			d[j] = dist[i][j];
		}
		distance_matrix.push_back(d);
	}

	return;
}




//
// Export to gml
//
void Graph::export_to_gml(const char* output) {

	ofstream file(output);
	file << "graph [\n";

	for (int i = 0; i < n_vertices; i++) {
		file << "node [\n";
		file << "\tid " << i << "\n";
		file << "\tlabel \"" << i << "\"\n";

		file << "\t graphics [ \n";
		file << "\t\t type \"ellipse\"\n";
		file << "\t\t hasFill 0 \n";
		file << "\t\t ] \n";

		file << "\t]\n" << endl;
	}
	int total_edges = 0;
	for (int i = 0; i < n_vertices; ++i) {
		for (int j = i+1; j < n_vertices; ++j) {
			if ( !is_adj(i, j) )
				continue;
			file << "edge [\n";
			file << "\t source " << i << "\n";
			file << "\t target " << j << "\n";
			file << "\t]\n";
			total_edges++;
		}
	}
	file << "\n]";
	file.close();
	cout << "TOTAL EDGES: " << total_edges << endl;
}


//
// Create an isomorphic graph according to a vertex mapping
// Mapping description: mapping[i] = position where vertex i is in new ordering
//
//
Graph::Graph(Graph* graph, vector<int>& mapping)
: n_vertices(graph->n_vertices), n_edges(graph->n_edges)
{
	// allocate adjacent matrix
	adj_m = new bool*[n_vertices];
	for (int i = 0; i < n_vertices; ++i) {
		adj_m[i] = new bool[n_vertices];
		memset(adj_m[i], false, sizeof(bool)*n_vertices);
	}

	// allocate adjacent list
	adj_list.resize(n_vertices);

	// construct graph according to mapping
	for (int i = 0; i < graph->n_vertices; ++i) {
		for (vector<int>::iterator it = graph->adj_list[i].begin();
				it != graph->adj_list[i].end();
				it++)
		{
			set_adj(mapping[i], mapping[*it]);
		}
	}
}


//
// Print graph
//
void Graph::print() {
	cout << "Graph" << endl;
	for (int v = 0; v < n_vertices; ++v) {
		if (adj_list[v].size() != 0) {
			cout << "\t" << v << " --> ";
			for (vector<int>::iterator it = adj_list[v].begin();
					it != adj_list[v].end();
					++it)
			{
				cout << *it << " ";
			}
			cout << endl;
		}
	}
	cout << endl;
}


