/*
 * --------------------------------------------------------
 * BDD Vertex Ordering class
 * --------------------------------------------------------
 */

#ifndef ORDERING_HPP_
#define ORDERING_HPP_

#include <cassert>
#include <cstdio>
#include <vector>
#include <cstdlib>

#include "mb_bdd.hpp"
#include "minband_instance.hpp"

using namespace std;


struct IntComparator {
	vector<int> &v;
	IntComparator(vector<int> &_v) : v(_v) { }
	bool operator()(int i, int j) {
		return (v[i] > v[j]);
	}
};


using namespace std;

enum OrderType {
  MinState, MaximalPath, CutVertexGen, CutVertex, Fixed,
  MinDegree, RootOrder, LexOrder, SpanningTree
};

// Class representing a general ordering
struct MB_Ordering {

	MinBandInst   *inst;
	char           name[256];
	OrderType	   order_type;

	MB_Ordering(MinBandInst* _inst, OrderType _order_type) : inst(_inst), order_type(_order_type) { }

	virtual ~MB_Ordering() { }

	// returns vertex corresponding to particular layer
	virtual int vertex_in_layer(BDD* bdd, int layer) = 0;
};


// Vertex that it is in the least number of states
struct RootOrdering : MB_Ordering {

	RootOrdering(MinBandInst *_inst) : MB_Ordering(_inst, RootOrder) {
		sprintf(name, "root_ordering");
	}

	int vertex_in_layer(BDD* bdd, int layer) { exit(0); }
};



// Vertex that it is in the least number of states
struct LexOrdering : MB_Ordering {

	LexOrdering(MinBandInst *_inst) : MB_Ordering(_inst, LexOrder) {
		sprintf(name, "lex_ordering");
	}

	int vertex_in_layer(BDD* bdd, int layer) { exit(0); }
};




// Ordering according to root node
struct MinInState : MB_Ordering {

	vector<bool> avail_v;   // available vertices

	MinInState(MinBandInst *_inst) : MB_Ordering(_inst, MinState) {
		avail_v.resize(inst->graph->n_vertices, true);
		sprintf(name, "min_in_state");
	}

	int vertex_in_layer(BDD* bdd, int layer);
};

// Minimum degree ordering
struct SpanningTreeOrdering : MB_Ordering {

  vector<int> v_in_layer;   // vertex at each layer

  SpanningTreeOrdering(MinBandInst *_inst) : MB_Ordering(_inst, SpanningTree) {
    sprintf(name, "spanningtree");
    construct_ordering();
  }

  int vertex_in_layer(BDD* bdd, int layer) {
    assert( layer >= 0 && layer < inst->graph->n_vertices);
    return v_in_layer[layer];
  }

private:
  void construct_ordering();
};

//// Vertex that it is in the least number of states - randomized!!
//struct RandomizedMinInState : IS_Ordering {
//
//	vector<bool> avail_v;   // available vertices
//	double prob;			// probability
//
//	boost::random::mt19937 gen;
//
//	RandomizedMinInState(MinBandInst *_inst, double _prob) : IS_Ordering(_inst, RandMinState), prob(_prob) {
//		avail_v.resize(inst->graph->n_vertices, true);
//		sprintf(name, "rand_min");
//
//		gen.seed(inst->graph->n_vertices + inst->graph->n_edges);
//	}
//
//	int vertex_in_layer(BDD* bdd, int layer);
//};


// Maximal Path Decomposition
struct MaximalPathDecomp : MB_Ordering {

	vector<int> v_in_layer;   // vertex at each layer

	MaximalPathDecomp(MinBandInst *_inst) : MB_Ordering(_inst, MaximalPath) {
		sprintf(name, "maxpath");
		construct_ordering();
	}

	int vertex_in_layer(BDD* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

private:
	void construct_ordering();
};



// Minimum degree ordering
struct MinDegreeOrdering : MB_Ordering {

  vector<int> v_in_layer;   // vertex at each layer

  MinDegreeOrdering(MinBandInst *_inst) : MB_Ordering(_inst, MinDegree) {
    sprintf(name, "mindegree");
    construct_ordering();
  }

  int vertex_in_layer(BDD* bdd, int layer) {
    assert( layer >= 0 && layer < inst->graph->n_vertices);
    return v_in_layer[layer];
  }
  
private:
  void construct_ordering();
};


//// Random ordering
//struct RandomOrdering : IS_Ordering {
//
//	vector<int> v_in_layer;   // vertex at each layer
//	RandomOrdering(MinBandInst *_inst) : IS_Ordering(_inst, Random) {
//		sprintf(name, "random");
//		v_in_layer.resize(inst->graph->n_vertices);
//		construct_ordering();
//	}
//
//	int vertex_in_layer(BDD* bdd, int layer) {
//		assert( layer >= 0 && layer < inst->graph->n_vertices);
//		return v_in_layer[layer];
//	}
//
//private:
//	void construct_ordering();
//};



// FixedOrdering: read from a file
struct FixedOrdering : MB_Ordering {

	vector<int> v_in_layer;   // vertex at each layer
	FixedOrdering(MinBandInst *_inst, char* filename) : MB_Ordering(_inst, Fixed) {

		sprintf(name, "fixed");
		v_in_layer.resize(inst->graph->n_vertices);
		read_ordering(filename);
	}

	int vertex_in_layer(BDD* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

private:
	void read_ordering(char* filename);
};

// Cut vertex decomposition
struct CutVertexDecompositionGeneralGraph : MB_Ordering {

	vector<int> v_in_layer;      // vertex at each layer
	bool** original_adj_matrix;

	CutVertexDecompositionGeneralGraph(MinBandInst *_inst) : MB_Ordering(_inst, CutVertexGen) {
		sprintf(name, "cut-vertex-gen");
		v_in_layer.resize(inst->graph->n_vertices);

		restrict_graph();
		construct_ordering();
		regenerate_graph();
	}

	int vertex_in_layer(BDD* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

private:
	void        restrict_graph();
	void        regenerate_graph();
	void        construct_ordering();
	void        identify_components(vector< vector<int> > &comps, vector<bool> &is_in_graph);
	vector<int> find_ordering(vector<bool> is_in_graph);
};


// Cut vertex decomposition
struct CutVertexDecomposition : MB_Ordering {

	vector<int> v_in_layer;   // vertex at each layer

	CutVertexDecomposition(MinBandInst *_inst) : MB_Ordering(_inst, CutVertex) {
		sprintf(name, "cut-vertex");
		v_in_layer.resize(inst->graph->n_vertices);
		construct_ordering();
	}

	int vertex_in_layer(BDD* bdd, int layer) {
		assert( layer >= 0 && layer < inst->graph->n_vertices);
		return v_in_layer[layer];
	}

private:
	void        construct_ordering();
	void        identify_components(vector< vector<int> > &comps, vector<bool> &is_in_graph);
	vector<int> find_ordering(vector<bool> is_in_graph);
};






#endif
