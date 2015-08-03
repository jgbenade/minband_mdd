/*
 * --------------------------------------------------------
 * Ordering class - implementation
 * --------------------------------------------------------
 */

#include <algorithm>
#include <fstream>
#include <cstdlib>

#include "minband_orderings.hpp"

using namespace std;


// min in state heuristic
int MinInState::vertex_in_layer(BDD* bdd, int layer) {
	return -1;
}

//// randomized min in state heuristic
//int RandomizedMinInState::vertex_in_layer(BDD* bdd, int layer) {
//	return -1;
//}

// fixed ordering
void FixedOrdering::read_ordering(char* filename) {
	ifstream ordering(filename);
	if( !ordering.is_open() ) {
		cout << "ERROR - could not open ordering file " << filename << endl;
		exit(1);
	}

	for( int v = 0; v < inst->graph->n_vertices; v++ ) {
		ordering >> v_in_layer[v];
	}

	ordering.close();
}

// minimum degree ordering
void SpanningTreeOrdering::construct_ordering() {

	v_in_layer.clear();

	vector<bool> selected( inst->graph->n_vertices, false);

	int v = rand() % inst->graph->n_vertices;
	selected[v] = true;
	v_in_layer.push_back(v);

	int previous_size = 1;
	vector<int> nodes_to_visit;

	while ( (int)v_in_layer.size() < inst->graph->n_vertices ) {
		int w=-1;

		for (vector<int>::const_iterator it = inst->graph->adj_list[v].begin(); it != inst->graph->adj_list[v].end(); ++it){

			w = (*it);
			//cout << "selected vertex " << w <<endl;

			if (!selected[w]){
				v_in_layer.push_back(w);
				selected[w] = true;
				nodes_to_visit.push_back(w);
			}
		}

		//if the graph is disconnected
		if (nodes_to_visit.size()>0){
			v = nodes_to_visit.back();
			nodes_to_visit.pop_back();
		} else{
			for (int i=0; i< inst->graph->n_vertices; i++){
				if (!selected[i]){
					v =i;
					v_in_layer.push_back(i);
					selected[i] = true;
				}
			}
		}

		previous_size = (int)v_in_layer.size();
	}

	cout<< endl<< "ordering ";
	for (vector<int>::iterator it = v_in_layer.begin(); it!= v_in_layer.end(); ++it){
		cout << *it;
	}cout<<endl;
}

void DegreeOrdering::construct_ordering(){
	v_in_layer.clear();

	vector<int> count(inst->graph->n_vertices,0);
	vector<bool> inserted(inst->graph->n_vertices, false);

	int vertex = -1;
	int max_count = -1;

	//get highest degree
	for (int i = 0; i< inst->graph->n_vertices; i++){
		if (inst->graph->degree(i)>max_count){
			max_count = inst->graph->degree(i);
			vertex = i;
		}
	}

	while (v_in_layer.size() < inst->graph->n_vertices){
		v_in_layer.push_back(vertex);
			inserted[vertex] = true;
			//update counts
			for (vector<int>::iterator neighbour = inst->graph->adj_list[vertex].begin();
					neighbour != inst->graph->adj_list[vertex].end(); ++neighbour){
				count[*neighbour]++;
			}
			max_count = -1;
			for (int i = 0; i< inst->graph->n_vertices; i++){
				int score = 1*inst->graph->degree(i) + 2* count[i];
				if (score > max_count
						and !inserted[i] ){
					max_count = score;
					vertex = i;
				}
			}
	}

	cout<<  "ordering ";
	for (vector<int>::iterator it = v_in_layer.begin(); it!= v_in_layer.end(); ++it){
		cout << *it;
	}cout<<endl;
}

void FromFrontOrdering::construct_ordering(){
	v_in_layer.clear();

	//get highest degree
	for (int i = 0; i< inst->graph->n_vertices; i++){
		v_in_layer.push_back(i);
	}
/*	for (vector<int>::iterator it = v_in_layer.begin(); it!= v_in_layer.end(); ++it){
		cout << *it;
	}cout<<endl;*/
}

void FrontBackOrdering::construct_ordering(){
	v_in_layer.clear();

	vector<bool> inserted(inst->graph->n_vertices, false);

	//get highest degree
	for (int i = 0; i< inst->graph->n_vertices; i++){
		if (!inserted[i]){
			v_in_layer.push_back(i);
			inserted[i] = true;
		}
		if (!inserted[inst->graph->n_vertices - i - 1]){
			v_in_layer.push_back(inst->graph->n_vertices - i - 1);
			inserted[inst->graph->n_vertices - i - 1] = true;
		}
	}
/*	for (vector<int>::iterator it = v_in_layer.begin(); it!= v_in_layer.end(); ++it){
			cout << *it;
		}cout<<endl;*/
}

// minimum degree ordering
void MinDegreeOrdering::construct_ordering() {

  v_in_layer.clear();

  // compute vertex degree
  vector<int> degree( inst->graph->n_vertices, 0 );
  for ( int i = 0; i < inst->graph->n_vertices; ++i ) {
    for ( int j = i+1; j < inst->graph->n_vertices; ++j ) {
      if ( inst->graph->adj_m[i][j] ) {
	++(degree[i]);
	++(degree[j]);
      }
    }
  }

  vector<bool> selected( inst->graph->n_vertices, false);

  while ( (int)v_in_layer.size() < inst->graph->n_vertices ) {
  
    int min = INF;
    int v = -1;

    for ( int i = 0; i < inst->graph->n_vertices; ++i ) {
      if ( degree[i] > 0 && degree[i] < min && !selected[i] ) {
	min = degree[i];
	v = i;
      }
    }

    if ( v == -1 ) {
      for( int i = 0; i < inst->graph->n_vertices; ++i ) {
	if ( !selected[i] ) {
	  //cout << "\n selected vertex " << i << " --> degree: " << degree[i] << endl;
	  v_in_layer.push_back( i );
	}
      }
    } else {
      selected[v] = true;
      v_in_layer.push_back(v);
      //cout << "\n selected vertex " << v << " --> degree: " << degree[v] << endl;
      for ( int i = 0; i < inst->graph->n_vertices; ++i ) {
	if ( i != v && inst->graph->adj_m[i][v] ) {
	  --(degree[i]);
	}
      }
    }
  }
}


// maximal path decomposition
void MaximalPathDecomp::construct_ordering() {

	int n_maximal_paths = 0;

	v_in_layer.resize(inst->graph->n_vertices);
	vector<bool> visited(inst->graph->n_vertices, false);

	int n = 0;  // number of vertices already considered in the path

	// partial orderings
	vector<int> left;
	vector<int> right;

	while( n < inst->graph->n_vertices ) {
		left.clear();
		right.clear();

		int middle = -1;
		// take first unvisited vertex
		for( int v = 0; v < inst->graph->n_vertices; v++ ) {
			if( !visited[v] ) {
				middle = v;
				break;
			}
		}
		visited[middle] = true;

		// right composition
		int current = middle;
		while( current != -1 ) {
			int next = -1;
			for( int v = 0; v < inst->graph->n_vertices; v++ ) {
				if( !visited[v] && inst->graph->is_adj(current, v) ) {
					next = v;
					right.push_back(next);
					visited[next] = true;
					break;
				}
			}
			current = next;
		}

		// left composition
		current = middle;
		while( current != -1 ) {
			int next = -1;
			for( int v = 0; v < inst->graph->n_vertices; v++ ) {
				if( !visited[v] && inst->graph->is_adj(current, v) ) {
					next = v;
					left.push_back(next);
					visited[next] = true;
					break;
				}
			}
			current = next;
		}

		// compose path from left to right
		for( int i = (int)left.size()-1; i>=0; i-- ) {
			v_in_layer[n++] = left[i];
		}
		v_in_layer[n++] = middle;
		for( int i = 0; i < (int)right.size(); i++ ) {
			v_in_layer[n++] = right[i];
		}
		n_maximal_paths++;
	}

	cout << "\nnumber of maximal paths in decomposition: " << n_maximal_paths << endl << endl;

	// quick check...
	if( n_maximal_paths == 1 ) {
		for( int v = 0; v < inst->graph->n_vertices-1; v++ ) {
			if( !inst->graph->is_adj(v_in_layer[v], v_in_layer[v+1]) ) {
				cout << "ERROR IN MAXIMAL PATH DECOMPOSITION\n";
				exit(1);
			}
		}
	}

}


//void RandomOrdering::construct_ordering() {
//	for( int i = 0; i < inst->graph->n_vertices; i++ ) {
//		v_in_layer[i] = i;
//	}
//
////
////	int a[] = {0,27,28,29,1,42,43,44,3,30,31,32,4,37,36,38,8,35,33,34,6,47,45,46,5,40,39,41,7,50,48,49,2,51,52,53,9,63,
////			64,65,10,78,79,80,12,66,67,68,13,73,72,74,17,71,69,70,15,83,81,82,14,76,75,77,16,86,84,85,11,87,88,89,18,99,
////			100,101,19,114,115,116,21,102,103,104,22,109,108,110,26,107,105,106,24,119,117,118,23,112,111,113,25,122,120,
////			121,20,123,124,125,54,55,56,57,58,59,60,61,62,90,91,92,93,94,95,96,97,98,126,127,128,129,130,131,132,133,134,
////			135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,
////			164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,
////			193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,
////			222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,
////			251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,
////			280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,
////			309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,
////			338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,
////			368,369,370,371,372,373,374,375,376,377 };
////
////	for( int i = 0; i < inst->graph->n_vertices; i++ ) {
////		v_in_layer[i] = a[i];
////	}
//
//
//	random_shuffle(v_in_layer.begin(), v_in_layer.end());
//}


void CutVertexDecomposition::identify_components(vector< vector<int> > &comps, vector<bool> &is_in_graph) {

	vector<int> label(inst->graph->n_vertices, -1);
	int num_comps = -1;

	vector<int> stack;

	vector<bool> visited(inst->graph->n_vertices, false);
	for( int i = 0; i < inst->graph->n_vertices; i++ ) {

		if( is_in_graph[i] && !visited[i]) {

			num_comps++;
			stack.push_back(i);

			while( !stack.empty() ) {

				int v = stack.back();
				stack.pop_back();

				label[v] = num_comps;
				visited[v] = true;

				for( int w = 0; w < inst->graph->n_vertices; w++ ) {
					if( w == v ) continue;
					if( is_in_graph[w] && inst->graph->is_adj(v,w) && !visited[w]) {
						stack.push_back(w);
					}
				}
			}
		}
	}

	comps.clear();
	comps.resize(num_comps+1);

	for( int v = 0; v < inst->graph->n_vertices; v++ ) {
		if( label[v] != -1 ) {
			comps[label[v]].push_back(v);
		}
	}
}

vector<int> CutVertexDecomposition::find_ordering(vector<bool> is_in_graph) {

	int size = 0;
	for( int i = 0; i < inst->graph->n_vertices; i++ ) {
		size += (is_in_graph[i] ? 1 : 0 );
	}

	// find vertex with all components less than half the size of the graph
	vector< vector<int> > comps;
	for( int i = 0; i < inst->graph->n_vertices; i++ ) {
		if( is_in_graph[i] ) {

			// try removing vertex
			is_in_graph[i] = false;
			identify_components(comps, is_in_graph);

			//cout << "componentes when removing vertex " << i << endl;

			bool all_valid = true;
			for( int j = 0; j < (int)comps.size() && all_valid; j++ ) {
				all_valid = ((int)comps[j].size() <= size/2);
				//cout << "\t" << comps[j].size() << endl;
			}

			if( all_valid ) {

				vector<int> ordering;

				// compose ordering for each component separately
				vector<bool> is_in_graph_new(inst->graph->n_vertices, false);
				for( int c = 0; c < (int)comps.size(); c++ ) {

					//cout << "*** Component " << c << endl;

					for( int v = 0; v < (int)comps[c].size(); v++ ) {
						is_in_graph_new[comps[c][v]] = true;
					}

					vector<int> order_bck = find_ordering(is_in_graph_new);

					for( int v = 0; v < (int)comps[c].size(); v++ ) {
						is_in_graph_new[comps[c][v]] = false;
						ordering.push_back(order_bck[v]);
					}
				}
				ordering.push_back(i);
				return ordering;
			}

			// put vertex back again
			is_in_graph[i] = true;
		}
	}
	return( vector<int>(1,-1) );
}




void CutVertexDecomposition::construct_ordering() {

	vector<bool> is_in_graph(inst->graph->n_vertices, true);
	vector<int> ordering = find_ordering(is_in_graph);

	//cout << endl;
	// for( int i = 0; i < (int)ordering.size(); i++ ) {
	//   cout << ordering[i] << " ";
	// }
	// cout << endl;
	//cout << "ordering size: " << ordering.size() << endl;

	for( int i = 0; i < (int)ordering.size(); i++ ) {
		v_in_layer[i] = ordering[i];
	}

}



void CutVertexDecompositionGeneralGraph::identify_components(vector< vector<int> > &comps, vector<bool> &is_in_graph) {

	vector<int> label(inst->graph->n_vertices, -1);
	int num_comps = -1;

	vector<int> stack;

	vector<bool> visited(inst->graph->n_vertices, false);
	for( int i = 0; i < inst->graph->n_vertices; i++ ) {

		if( is_in_graph[i] && !visited[i]) {

			num_comps++;
			stack.push_back(i);

			while( !stack.empty() ) {

				int v = stack.back();
				stack.pop_back();

				label[v] = num_comps;
				visited[v] = true;

				for( int w = 0; w < inst->graph->n_vertices; w++ ) {
					if( w == v ) continue;
					if( is_in_graph[w] && inst->graph->is_adj(v,w) && !visited[w]) {
						stack.push_back(w);
					}
				}
			}
		}
	}

	comps.clear();
	comps.resize(num_comps+1);

	for( int v = 0; v < inst->graph->n_vertices; v++ ) {
		if( label[v] != -1 ) {
			comps[label[v]].push_back(v);
		}
	}
}

vector<int> CutVertexDecompositionGeneralGraph::find_ordering(vector<bool> is_in_graph) {

	int size = 0;
	for( int i = 0; i < inst->graph->n_vertices; i++ ) {
		size += (is_in_graph[i] ? 1 : 0 );
	}

	// find vertex with all components less than half the size of the graph
	vector< vector<int> > comps;
	for( int i = 0; i < inst->graph->n_vertices; i++ ) {
		if( is_in_graph[i] ) {

			// try removing vertex
			is_in_graph[i] = false;
			identify_components(comps, is_in_graph);

			//cout << "componentes when removing vertex " << i << endl;

			bool all_valid = true;
			for( int j = 0; j < (int)comps.size() && all_valid; j++ ) {
				all_valid = ((int)comps[j].size() <= size/2);
				//cout << "\t" << comps[j].size() << endl;
			}

			if( all_valid ) {

				vector<int> ordering;

				// compose ordering for each component separately
				vector<bool> is_in_graph_new(inst->graph->n_vertices, false);
				for( int c = 0; c < (int)comps.size(); c++ ) {

					//cout << "*** Component " << c << endl;

					for( int v = 0; v < (int)comps[c].size(); v++ ) {
						is_in_graph_new[comps[c][v]] = true;
					}

					vector<int> order_bck = find_ordering(is_in_graph_new);

					for( int v = 0; v < (int)comps[c].size(); v++ ) {
						is_in_graph_new[comps[c][v]] = false;
						ordering.push_back(order_bck[v]);
					}
				}
				ordering.push_back(i);
				return ordering;
			}

			// put vertex back again
			is_in_graph[i] = true;
		}
	}
	return( vector<int>(1,-1) );
}




void CutVertexDecompositionGeneralGraph::construct_ordering() {

	vector<bool> is_in_graph(inst->graph->n_vertices, true);
	vector<int> ordering = find_ordering(is_in_graph);

	//cout << endl;
	// for( int i = 0; i < (int)ordering.size(); i++ ) {
	//   cout << ordering[i] << " ";
	// }
	// cout << endl;
	//cout << "ordering size: " << ordering.size() << endl;

	for( int i = 0; i < (int)ordering.size(); i++ ) {
		v_in_layer[i] = ordering[i];
	}

}


void CutVertexDecompositionGeneralGraph::restrict_graph() {
/*
	vector< pair<int,int> > edges;

	vector<int> vertices(inst->graph->n_vertices);
	vector<int> degrees(inst->graph->n_vertices);
	for( int i = 0; i < inst->graph->n_vertices; i++ ) {
		vertices[i] = i;
		degrees[i] = inst->graph->adj_list[i].size()-1;
	}
	IntComparator int_comp(degrees);
	sort(vertices.begin(), vertices.end(),int_comp);

	vector<bool> is_taken(inst->graph->n_vertices, false);
	vector<int> taken;

	cout << "vertices ordered by degree: " << endl;
	for( int i = 0; i < inst->graph->n_vertices; i++ ) {
		cout << vertices[i] << " ";
	}
	cout << endl;

	is_taken[vertices[0]] = true;
	int n_taken = 1;

	while( n_taken != inst->graph->n_vertices ) {
		for( int i = 0; i < inst->graph->n_vertices; i++ ) {

			if( !is_taken[vertices[i]] )
				continue;

			int w = vertices[i];

			bool new_edge = false;

			for( int v = 0; v < inst->graph->n_vertices; v++ ) {
				if( v != w && inst->graph->is_adj(v, w) && !is_taken[v] ) {

					pair<int,int> t;
					t.first = v;
					t.second = w;
					edges.push_back(t);

					is_taken[v] = true;
					n_taken++;

					new_edge = true;
					cout << "took " << v << " due to " << w << endl;
				}
			}

			if( new_edge )
				break;
		}
	}

	cout << endl;
	for( int i = 0; i < (int)edges.size(); i++ ) {
		cout << edges[i].first << "," << edges[i].second << endl;
	}
	cout << endl;

	bool** new_adj = new bool*[inst->graph->n_vertices];
	for( int i = 0; i < inst->graph->n_vertices; i++ ) {
		new_adj[i] = new bool[inst->graph->n_vertices];
	}

	for( int i = 0; i < inst->graph->n_vertices; i++ ) {
		for( int j = i+1; j < inst->graph->n_vertices; j++ ) {
			new_adj[i][j] = false;
			new_adj[j][i] = false;
		}
	}

	for( int i = 0; i < (int)edges.size(); i++ ) {
		new_adj[edges[i].first][edges[i].second] = true;
		new_adj[edges[i].second][edges[i].first] = true;
	}

	original_adj_matrix = inst->graph->adj_m;
	inst->graph->adj_m = new_adj;*/
}


void CutVertexDecompositionGeneralGraph::regenerate_graph() {
	inst->graph->adj_m = original_adj_matrix;
}



