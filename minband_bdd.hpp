/*
 * ===============================================================
 * MinBandBDD - Definition class
 * ===============================================================
 */

#ifndef MINBANDBDD_HPP_
#define MINBANDBDD_HPP_

#include <map>
#include <queue>
#include <deque>
#include <set>

#include "minband_orderings.hpp"
#include "minband_instance.hpp"
#include "util.hpp"

#define EPS 1e-8

using namespace std;


//
// Parameters
// TODO: include this in a parameter struct
//
#define ROOT_WIDTH	100				// width for initial relaxation



// -------------------------------------------------------------------
// Branching Node
// -------------------------------------------------------------------

struct BranchNode {
	vector<set<int> >		state;				// state to be explored
	int   					cost;				// node cost
	int   					relax_lb;           // upper bound obtained for the BDD where node was created
	vector<int>				vertex;
	//
	// Constructor from a BDD node
	//
	BranchNode(Node* _node, double _cost)
	:	cost(_cost), relax_lb(0)
	{
		const State& state_v = _node->state;
		for (vector<set<int> >::const_iterator it = state_v.begin(); it != state_v.end(); ++it) {
			state.push_back(*it);
		}
	}

	//
	// Constructor from a state vector
	//
	BranchNode(const vector<set<int> >& _state, double _cost)
	:	state(_state), cost(_cost), relax_lb(0)
	{
	}

	//
	// Empty constructor
	//
	BranchNode() {
	}

	void printState();
};

inline void BranchNode::printState(){
	for( std::vector<set<int> >::const_iterator i = state.begin(); i != state.end(); ++i){
		cout << "/";
		for( std::set<int>::const_iterator j = (*i).begin(); j != (*i).end(); ++j){
			std::cout << *j ;
		}

	}
	cout << ":" << cost << ":"<< relax_lb << endl;
}

//
// Branch node comparator by relaxed UB
//
struct BranchNodeComparatorUB {
	//we want to first branch on the lower cost nodes, if equal cost then thosehigher up the tree
	BranchNodeComparatorUB() { }
	bool operator()(const BranchNode* b1, const BranchNode* b2) {
		if (b1->relax_lb == b2->relax_lb) {
			if (b1->state.size() == b2->state.size()) {
				return b1->cost > b2->cost;
			}
			return b1->state.size() > b2->state.size();
		}
		return b1->relax_lb > b2->relax_lb;
	}
};

//
// Branch node comparator for DFS
//
struct BranchNodeComparatorDFS {
	BranchNodeComparatorDFS() { }
	bool operator()(const BranchNode* b1, const BranchNode* b2) {
		//TODO
		return true;
		if (b1->state.size() == b2->state.size()) {
			if (b1->relax_lb == b2->relax_lb) {
				return b1->cost < b2->cost;
			}
			return b1->relax_lb < b2->relax_lb;
		}
		return b1->state.size() > b2->state.size();
	}
};

//comparotor for ILP relaxation triples
struct mytriplecompfv {
	inline bool operator() (const vector<int>& a, const vector<int>& b) {
		return a[1] > b[1];
	}
};
  	 //want to sort with smallest f_v at back
struct mytriplecomplv {
	inline bool operator() (const vector<int>& a, const vector<int>& b) {
		return a[2] > b[2];
	}
};

//
// Queue of branch nodes ordered by relax ub
//
typedef priority_queue<BranchNode*, vector<BranchNode*>, BranchNodeComparatorUB> BranchNodeQueue;


//
// Queue of branch nodes ordered by state size for DFS search
//
typedef priority_queue<BranchNode*, vector<BranchNode*>, BranchNodeComparatorDFS> BranchNodeDFS;


// -------------------------------------------------------------------
// Independent Set BDD: Relaxation / Restriction generator
// -------------------------------------------------------------------

// Independent set solver class
class MinBandBDD {

public:

	//
	// Constructor given a DIMACS instance filename
	//
	// MinBandBDD(char* filename, int _maxWidth)
	// DDX10:
	MinBandBDD(const int _rootWidth,
				const int _ddWidth,
				const string & _instanceName);

	// Generate relaxation. Returns true iff relaxation is exact.
	int generateRelaxation(const int initial_weight);

	// Generate relaxation given a Branch node
	int generateRelaxation(const BranchNode* branch_node);

	// Generate restriction
	int generateRestriction(const int initial_weight);

	// Generate restriction given a Branch node
	int generateRestriction(const BranchNode* branch_node);

	// Update local upper bound, indicating whether it was found locally or remotely
	void updateLocalUB(int _ub, bool isLocal = false);

	// Verify if relaxed BDD is exact
	bool isBDDExact() 				{ return isExact; }

	// Verify if the upper bound has been updated
	bool wasUBUpdated()				{ return isUBUpdated  ; } //we reset ubupdated in genrelaxation

	// Get best lower bound found
	int getBestLB() 		  		{ return best_lb; }
	//int getBestUB()					{ return best_ub; }

	// Get last computed upper bound
	int getUB() 		  			{ return upper_bound; }

	// Add branching nodes to queue
    void addBranchNodesQueue(BranchNodeQueue& queue, unsigned long int &size_pool);

	// Add branching nodes to DFS queue
	void addBranchNodesDFS(BranchNodeDFS& queue);

	// Add branching nodes to DFS queue
    void addBranchNodesDFS(vector<BranchNode*>& queue);

  	int calculateCost(Node* _node);
  	int calculateCost_bounds(Node* _node);
  	int calculateCost_bounds_fast(Node* _node);
  	int calculateCost_caprara(Node* _node);
  	int calculateCost_mu2(Node* _node);
  	int calculateCost_mu2_fast(Node* _node);
  	int calculateCost_ILP(Node* _node);

  	//bool mytriplecomp (const vector<int>& a, const vector<int>& b);

	// Public data members
	MinBandInst* 			inst;				// independent set instance
	vector<int> 			vertex_in_layer;  	// which vertex a layer corresponds to
	Node*					root_node;		  	// relaxation root node

private:

	// ----------------------------------------------
	// Attributes
	// ----------------------------------------------

	BDDNodePool 		node_pool;		  		// list of nodes to explore

	vector<int>			active_vertices;  		// active vertices in the relaxation
	vector<Node*> 		nodes_layer;	  		// nodes in a layer
	State				active_state; 			//the state of the active vertex

	int* 				in_state_counter;  		// min-in-state vertex ordering parameters
	MB_Ordering* 		var_ordering; 			// BDD Vertex ordering

	int 				maxWidth;  				// max allowed width
	int 				best_lb;				// best lower bound
	int					best_ub;				// best upper bound
	bool				isExact;				// if BDD is exact
	State				best_ub_state; 			// the current best heuristic solution

	bool				isUBUpdated;			// if upper bound has been updated

	vector<BranchNode*> branch_nodes; 			// BDD nodes that require branching

	vector<int>			root_ordering;			// relative ordering of each vertex

	int					upper_bound;			// last computed upper bound
	int 				lower_bound;			// last computed lower bound

	int					initialPosPool; 	// pointer to first position in branching pool

    bool variable_width;  // if BDD maximum width is variable

    vector<vector<int> > edges_to_check;         //given the ordering, the mapped corresponding edges.

	// Auxiliary parameters
	CompareNodesCost bn_comparator;


	// ----------------------------------------------
	// Functions
	// ----------------------------------------------

	// Choose vertex
	int choose_next_vertex_min_size_next_layer();

	// Merge two nodes
	void merge(Node* nodeA, Node* nodeB);

	// Merge layer
	void mergeLayer(int layer, vector<Node*> &nodes_layer);

	// Restrict layer
	void restrictLayer(int layer, vector<Node*> &nodes_layer);

	// Add new local node to be explored later
	void addBranchNode(Node* _node);
};


// -------------- Utilities ----------------

struct CompareOrderDescending {
	vector<int> &order;
	CompareOrderDescending(vector<int>& _order) : order(_order) { }
	bool operator()(const int a, const int b) {
		return order[a] > order[b];
	}
};



// ----------------------------------------------------
// Inline implementations
// ----------------------------------------------------

//
// Merge node A into node B. Node A is always deleted.
//
inline void MinBandBDD::merge(Node* nodeA, Node* nodeB) {
	assert( (nodeA->state) == (nodeB->state) );

	if ( (nodeA->exact == nodeB->exact) ) {

		// Both nodes are either relaxed or exact
		nodeB->cost = MIN(nodeA->cost, nodeB->cost);

	} else if (nodeA->exact) {

		// Node A is exact, B is relaxed
		branch_nodes.push_back( new BranchNode(nodeA, nodeA->cost) );
		nodeB->cost = MIN(nodeA->cost, nodeB->cost);

	} else if (nodeB->exact) {

		// Node B is exact, A is relaxed
		branch_nodes.push_back( new BranchNode(nodeB, nodeB->cost));

		nodeB->exact = false;
		nodeB->cost = MIN(nodeA->cost, nodeB->cost);
	}

	delete nodeA;
}

//
// Update lower bound, indicating whether it was found locally or remotely. If locally,
// then mark that bound was updated to send it to other nodes
//
inline void MinBandBDD::updateLocalUB(int _ub, bool isLocal) {
	if (_ub < best_ub) {
		best_ub = _ub;
		if (isLocal) {
			isUBUpdated = true;
		}
	}
}


//
// Add new local node to be explored
//
inline void MinBandBDD::addBranchNode(Node* _node) {
	//cout << "Adding branch node " ; _node->printState();
	branch_nodes.push_back(	new BranchNode(_node, _node->cost) );
}

//
// Add branching nodes to queue
// TODO: apply bounding heuristic
//
inline void MinBandBDD::addBranchNodesQueue(BranchNodeQueue& queue, unsigned long int &size_pool) {

	for (vector<BranchNode*>::iterator st = branch_nodes.begin(); st != branch_nodes.end(); ++st) {
	//todo just check what is relax_lb actually
		//(*st)->printState();
		if ( (*st)->relax_lb >= best_ub ) {
			//cout<< "Not adding" ;		(*st)->printState();
			delete (*st);

		} else {
		  size_pool += (*st)->state.size();
		  queue.push(*st);
		  // cout << "Adding "; 		(*st)->printState();

		}
	}

	branch_nodes.clear();
}


//
// Add branching nodes to DFS queue
// TODO: apply bounding heuristic
//
inline void MinBandBDD::addBranchNodesDFS(BranchNodeDFS& queue) {

	//TODO
	/*for (vector<BranchNode*>::iterator st = branch_nodes.begin(); st != branch_nodes.end(); ++st) {
		if ( (*st)->relax_lb <= best_lb ) {
			delete (*st);
		} else {
			queue.push(*st);
		}
	}
	branch_nodes.clear();*/
}

//
// Add branching nodes to DFS queue
// TODO: apply bounding heuristic
//
inline void MinBandBDD::addBranchNodesDFS(vector<BranchNode*>& queue) {
	for (vector<BranchNode*>::iterator st = branch_nodes.begin(); st != branch_nodes.end(); ++st) {
		if ( (*st)->relax_lb <= best_lb ) {
			delete (*st);
		} else {
			queue.push_back(*st);
		}
	}
	branch_nodes.clear();
}


//
// Generate relaxation given a branch node
//
inline int MinBandBDD::generateRelaxation(const BranchNode* branch_node) {
	active_state = branch_node->state;
	return generateRelaxation(branch_node->cost);
}

//
// Generate restriction given a branch node
//
inline int MinBandBDD::generateRestriction(const BranchNode* branch_node) {
	active_state = branch_node->state;
	return generateRestriction(branch_node->cost);
}


#endif /* INDEPSETBDD_HPP_ */
