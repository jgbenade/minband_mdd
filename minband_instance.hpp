/*
 * --------------------------------------------------------
 * Instance data structures
 * --------------------------------------------------------
 */

#ifndef INSTANCE_HPP_
#define INSTANCE_HPP_

#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
//#include <boost/dynamic_bitset.hpp>

using namespace std;


//
// Graph structure
//
struct Graph {

	int                         n_vertices;         /**< |V| */
	int                         n_edges;            /**< |E| */
	double*						weights;			/**< weight of each vertex */

	bool**                      adj_m;              /**< adjacent matrix */
	vector< vector<int> >       adj_list;           /**< adjacent list */

	vector<vector<vector<int> > > vertex_neighbourhood; /** stores the v_n[v][k] are the vertices
															in the k-th neigbourhood of v**/
	vector<vector<int> >		distance_matrix	;	/**pairwise distances**/


	/** Set two vertices as adjacents */
	void set_adj(int i, int j);

	/** Check if two vertices are adjancent */
	bool is_adj(int i, int j);

	/** Empty constructor */
	Graph();

	/** Create an isomorphic graph according to a vertex mapping */
	Graph(Graph* graph, vector<int>& mapping);

	/** Read graph from a DIMACS format */
	void read_dimacs(const char* filename);

	/** Export to GML format */
	void export_to_gml(const char* output);

	/** Constructor with number of vertices */
	Graph(int num_vertices);

	/** Add edge */
	void add_edge(int i, int j);

	/** Remove edge */
	void remove_edge(int i, int j);

	/** Return degree of a vertex */
	int degree( int v ) { return adj_list[v].size(); }

	/** Print graph */
	void print();

	int dist(int i, int j);

	void calculate_pairwise_dist()	;



	vector<vector<int> > make_layered_graph(int v_);
};


using namespace std;
//
// Independent set instance structure
//
struct MinBandInst {

	Graph*              				graph;             	// independent set graph
	//vector< boost::dynamic_bitset<> >	adj_mask_compl;	 	// complement mask of adjacencies


	/** Read DIMACS independent set instance */
	void read_DIMACS(const char* filename);

	int calculate_Degree_Bound();

	int calculate_Caprara_Bound();

	int calculate_HalfDensity_Bound();

	int 						degree_bound;       /** max degree /2 **/
	int							caprara_bound;		/** Bound by Caprara, Salazar-Gonzalez**/
	int 						half_density_bound; /** 1/2 approx of density bound (Blum et al) **/

	vector<int>					caprara_list;       /**list stores max \lceil (N_k(v)-1)/k\rceil for every v
		 	 	 	 	 	 	 	 	 	 	 	 	 	 Used for computing alpha_1^S efficiently**/

	int CuthillMckee();								/** returns cuthillmckee banwidth**/
};



/*
 * -----------------------------------------------------
 * Inline implementations: Graph
 * -----------------------------------------------------
 */


/**
 * Empty constructor
 */
inline Graph::Graph() : n_vertices(0), n_edges(0), weights(NULL), adj_m(NULL) {
}


/**
 * Constructor with number of vertices
 **/
inline Graph::Graph(int num_vertices)
: n_vertices(num_vertices), n_edges(0), weights(NULL)
{
	adj_m = new bool*[ num_vertices ];
	for (int i = 0; i < num_vertices; ++i) {
		adj_m[i] = new bool[ num_vertices ];
		memset( adj_m[i], false, sizeof(bool) * num_vertices );
	}
	adj_list.resize(num_vertices);
}


/**
 * Check if two vertices are adjacent
 */
inline bool Graph::is_adj(int i, int j) {
/*
	assert(i >= 0);
	assert(j >= 0);
	assert(i < n_vertices);
	assert(j < n_vertices);
*/
	return adj_m[i][j];
}


/**
 * Set two vertices as adjacent
 */
inline void Graph::set_adj(int i, int j) {
/*
	assert(i >= 0);
	assert(j >= 0);
	assert(i < n_vertices);
	assert(j < n_vertices);
*/

	// check if already adjacent
	if (adj_m[i][j])
		return;

	// add to adjacent matrix and list
	adj_m[i][j] = true;
	adj_list[i].push_back(j);
}



/**
 * Add edge
 **/
inline void Graph::add_edge(int i, int j) {
	/*assert(i >= 0);
	assert(j >= 0);
	assert(i < n_vertices);
	assert(j < n_vertices);
*/
	// check if already adjacent
	if (adj_m[i][j])
		return;

	// add to adjacent matrix and list
	adj_m[i][j] = true;
	adj_m[j][i] = true;
	adj_list[i].push_back(j);
	adj_list[j].push_back(i);

	n_edges++;
}

/**
 * Remove edge
 **/
inline void Graph::remove_edge(int i, int j) {
/*
	assert(i >= 0);
	assert(j >= 0);
	assert(i < n_vertices);
	assert(j < n_vertices);
*/

	// check if already adjacent
	if (!adj_m[i][j])
		return;

	// add to adjacent matrix and list
	adj_m[i][j] = false;
	adj_m[j][i] = false;

	for (int v = 0; v < (int)adj_list[i].size(); ++v) {
		if ( adj_list[i][v] == j ) {
			adj_list[i][v] = adj_list[i].back();
			adj_list[i].pop_back();
			break;
		}
	}

	for (int v = 0; v < (int)adj_list[j].size(); ++v) {
		if ( adj_list[j][v] == i ) {
			adj_list[j][v] = adj_list[j].back();
			adj_list[j].pop_back();
			break;
		}
	}
}


inline int Graph::dist(int i, int j){
	return distance_matrix[i][j];
}

inline int MinBandInst::calculate_Degree_Bound(){
	int largest_degree = 0;
	for (int v =0; v < graph->n_vertices; v++){
		largest_degree = std::max((int)graph->adj_list[v].size(), largest_degree);
	}
	degree_bound =  (int)std::ceil(1.*largest_degree/2);
	return (int)std::ceil(1.*largest_degree/2);
}

inline int MinBandInst::calculate_Caprara_Bound(){
	int current_bound = graph->n_vertices;
	caprara_list.clear();

	for (int v = 0; v < graph->n_vertices; v++){
		vector<vector<int> > layered_graph = graph->vertex_neighbourhood[v];
		int cumulative_vertices = 0;
		int internal_max = 0;

		for (vector<vector<int> >::iterator layer = layered_graph.begin()+1; layer!=layered_graph.end(); ++layer){
			cumulative_vertices = cumulative_vertices + (*layer).size();
			int layer_number = layer - layered_graph.begin();

			if (layer_number > 0){
				internal_max = std::max(internal_max, (int)std::ceil((1.*cumulative_vertices-1)/layer_number));
			}
		}

		caprara_list.push_back(internal_max);

		//account for case layer 0
		if (internal_max > 0)
			current_bound = (int)std::min(current_bound, internal_max);
	}

	cout << "Caprara_list: " ;
	for (vector<int>::const_iterator it = caprara_list.begin(); it != caprara_list.end(); ++it)
		cout << *it << ",";

	caprara_bound = current_bound;
	//cout<< "Caprara bound= "<< current_bound<< endl;
	return caprara_bound;
}

inline int MinBandInst::calculate_HalfDensity_Bound(){
	double current_bound = 0;

	for (int v = 0; v<graph->n_vertices; v++){
		// cout << "making layered graph at=" <<v <<endl;
		vector<vector<int> > layered_graph = graph->vertex_neighbourhood[v];
		//cout << "done layered at="<<v<<endl;

		/*for(vector<vector<int> >::const_iterator l = layered_graph.begin(); l!= layered_graph.end(); ++l){
			for (vector<int>::const_iterator ll = (*l).begin(); ll!= (*l).end(); ++ll){
				cout << (*ll);
			} cout << ",";
		}cout << endl;*/

		int cumulative_vertices = 0;
		double internal_max = 0;

		for (vector<vector<int> >::iterator layer = layered_graph.begin()+1; layer != layered_graph.end(); ++layer){

			cumulative_vertices = cumulative_vertices + (*layer).size();
			int layer_number = layer - layered_graph.begin();
			//cout << cumulative_vertices << "\t" << layer_number << endl;

			if (layer_number >0){
				internal_max = std::max(internal_max, (double)((1.*cumulative_vertices-1)/(2*layer_number)) );
			}

		}
		//cout << "nov=" <<cumulative_vertices << "nol=" << layered_graph.size() << endl;
		current_bound = std::max(current_bound, internal_max);
	}

	half_density_bound = (int)std::ceil(current_bound);
	//cout<< "Half-density bound= "<< current_bound<< endl;
	return half_density_bound;
}

inline vector<vector<int> > Graph::make_layered_graph(int v_){
	//assert(v_ < n_vertices);

	vector<int> visited(n_vertices, 0);
	vector<vector<int> > layered_graph;
	vector<int> this_layer;
	this_layer.push_back(v_);
	visited[v_] = 1;

	while (this_layer.size() > 0){

		layered_graph.push_back(this_layer);
		vector<int> new_layer;

		for (vector<int>::const_iterator it = this_layer.begin(); it!=this_layer.end(); ++it){

			for(vector<int>::const_iterator nb = adj_list[(*it)].begin(); nb != adj_list[(*it)].end(); ++nb){

				if (visited[(*nb)] == 0){
					new_layer.push_back((*nb));
					visited[(*nb)] = 1;
				}

			}

		}

		//layered_graph.push_back(new_layer);
		this_layer = new_layer;
	}
	return layered_graph;
}


/*
 * -----------------------------------------------------
 * Inline implementations: Independent Set
 * -----------------------------------------------------
 */


//
// Read DIMACS independent set instance with no costs
//
inline void MinBandInst::read_DIMACS(const char *filename) {

	cout << "\nReading instance " << filename << endl;

	// read graph
	graph = new Graph;
	graph->read_dimacs(filename);

	cout << "\tnumber of vertices: " << graph->n_vertices << endl;
	cout << "\tnumber of edges: " << graph->n_edges << endl;

//	// create complement mask of adjacencies
//	adj_mask_compl.resize(graph->n_vertices);
//	for( int v = 0; v < graph->n_vertices; v++ ) {
//		adj_mask_compl[v].resize(graph->n_vertices, true);
//		for( int w = 0; w < graph->n_vertices; w++ ) {
//			if( graph->is_adj(v,w) ) {
//				adj_mask_compl[v].set(w, false);
//			}
//		}
//		// we assume here a vertex is adjacent to itself
//		adj_mask_compl[v].set(v, false);
//	}

	cout << "\tdone.\n" << endl;

}

inline int MinBandInst::CuthillMckee(){
	int best_solution = graph->n_vertices;

	vector<int> ordering;
	vector<bool> inserted(graph->n_vertices, false);

	//use all vertices as root
	for (int i = 0; i < graph->n_vertices; i++){
		ordering.clear();
		inserted.clear();
		inserted.resize(graph->n_vertices);
		//ordering.push_back(i);
		//inserted[i] = true;

		vector<vector<int> >::const_iterator layer_it =graph->vertex_neighbourhood[i].begin();
		while(layer_it != graph->vertex_neighbourhood[i].end())	{

			//every vertex in this layer must be inserted
			for(int j =0; j < (*layer_it).size(); j++){

				int smallest_deg = graph->n_vertices+1;
				int smallest_vertex = -1;

				//find minimum degree vertex
				for(vector<int>::const_iterator vertex_it = (*layer_it).begin(); vertex_it!= (*layer_it).end(); vertex_it++){

					if ((graph->degree(*vertex_it) < smallest_deg) && !inserted[(*vertex_it)]){
						smallest_deg = graph->degree(*vertex_it);
						smallest_vertex = *vertex_it;
					}

				}

				//insert minimum degree vertex
				if (smallest_vertex > -1){
					ordering.push_back(smallest_vertex);
					inserted[smallest_vertex] = true;
				}
				else cout << "ERROR CuthillMckee" << endl;

			}

			++layer_it;
		}

		/*//ptint ordering
		cout << "\nOrdering: ";
		for (vector<int>::const_iterator it = ordering.begin(); it != ordering.end(); ++it){
			cout <<*it <<",";
		} cout << endl;*/

		//ordering has been created
		// calculate bandwidth of ordering
		int bandwidth = 0;
		for (int u = 0; u < graph->n_vertices -1; u++){
			for (int v = u+1; v < graph->n_vertices; v++){
				if (graph->adj_m[ordering[u]][ordering[v]]){
					if(bandwidth < std::abs(u-v))	bandwidth = std::abs(u-v);
				}
			}
		}
		best_solution = std::min(best_solution, bandwidth);

		//cout<< "CM root "<< i << ": "<<bandwidth;
	}
	return best_solution;
}

#endif /* INSTANCE_HPP_ */
