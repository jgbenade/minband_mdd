/*
 * ===============================================================
 * MinBandBDD - Implementation class
 *
 * TODO: memory for branching pool can be fixed, since there is
 *       a maximum number of nodes we will add at a time
 * ===============================================================
 */

#include <algorithm>
#include <iostream>
#include <cmath>
#include "minband_bdd.hpp"

using namespace std;
//TODO 1 Merging nodes: DONE
//			Merging algorithm DONE
//			Do states ever match??? yes they do DONE
//TODO 2 Branching order, at least a min spanning tree or something DONE
//TODO 3 Restrictions DONE
//TODO 4 General search part
//todo 5 use lower bounds to speed up cost calc, maybe only fixed nodes (bounds on rest)
//			estimate if you use bounds
//TODO 6 Proper alldiff filtering
//TODO 7 the upperbound cost to prevent branching on min value (take max, min in neighbour domains)
//todo 8 bitset domains...


//
// Constructor given a DIMACS instance filename
//
// MinBandBDD(char* filename, int _maxWidth)
// DDX10:
MinBandBDD::MinBandBDD(const	 int _rootWidth,
			 const int _ddWidth,
			 const string & _instanceName)
	: inst( new MinBandInst ),
	  best_lb(0),
	  isExact(false),
	  isUBUpdated(false),
	  initialPosPool(0),
	  variable_width(false)
{
  //assert( _ddWidth == INF_WIDTH || _ddWidth > 0 );

	inst->read_DIMACS(_instanceName.c_str());
	best_ub = inst->graph->n_vertices;
	//in_state_counter = new int[inst->graph->n_vertices];
	//inst->graph->print();

	/*//testing calculateCosts
	vertex_in_layer.clear();
	vertex_in_layer.push_back(0);
	vertex_in_layer.push_back(3);
	vertex_in_layer.push_back(4);
	vertex_in_layer.push_back(2);
	State s;
	set<int> d1,d2;
	d1.insert(2);
	d2.insert(1);
	d2.insert(3);
	set<int> d3(d2);

	set<int> d4(d1);
	d4.insert(6);
	d4.insert(7);

	set<int> d5(d4); d5.insert(5);
	s.push_back(d4);
	s.push_back(d2);
	s.push_back(d3);
	s.push_back(d1);
	s.push_back(d5);
	Node n1(s,0);
	n1.printState();
	cout << calculateCost(n1);*/

	// ============================================================
	// Generate initial relaxation to derive variable ordering
	// TODO: save initial UB
	// ============================================================

	cout << "[BDD] Root node computation..." << endl;

	if (_rootWidth == -1) {
	  variable_width = true;
	  cout << "\tRoot width: size of state" << endl;
	  maxWidth = INF;
	} else {
	  cout << "\tRoot width: " << _rootWidth << endl;
	  maxWidth = _rootWidth;

	}
	//maxWidth = (_ddWidth == -1? INF : _ddWidth);
	// set min-in-state variable ordering
	//var_ordering = new LexOrdering(inst);
	var_ordering = new SpanningTreeOrdering(inst);

	cout << " ### Degree bound\t" << inst->graph->calculate_Degree_Bound() << endl;
	cout << " ### Half-density bound\t " << inst->graph->calculate_HalfDensity_Bound() << endl;
	cout << " ### Caprara bound   \t " << inst->graph->calculate_Caprara_Bound() << endl;

	// generate relaxation

	active_vertices.resize(inst->graph->n_vertices);
	for (int i = 0; i < inst->graph->n_vertices; ++i) {
		active_vertices[i] = i;
	}

	best_ub = inst->graph->n_vertices;
	upper_bound = inst->graph->n_vertices;

	lower_bound = -1;
	int ub1 = generateRestriction(inst->graph->n_vertices);
	upper_bound = ub1;
	int lb1 = generateRelaxation(-1);
	best_lb = lb1;
	//int ub1 = generateRestriction(inst->graph->n_vertices);

	cout << " ### Lower bound: " << lower_bound << endl;
	cout << " ### Upper bound: " << upper_bound << endl;

	// check if BDD is already exact
	isExact = branch_nodes.empty();
	if (isExact) {
		return;
	}

	// initialize branch node pool
	//update bounds because the higher level search repeatedly copy form the pool
	for (vector<BranchNode*>::iterator st = branch_nodes.begin(); st != branch_nodes.end(); ++st) {
		//(*st)->relax_lb = MIN(ub1, (*st)->cost + (int)(*st)->state.size());
		(*st)->relax_lb = MAX(lb1, (*st)->cost);
	}

	// set new variable ordering
	//	var_ordering = new RootOrdering(inst);
	//root_ordering.resize(inst->graph->n_vertices);
	//for (size_t i = 0; i < vertex_in_layer.size(); ++i) {
	//	root_ordering[vertex_in_layer[i]] = i;
	//}

	// reinitialize width
	if (_ddWidth == -1) {
	  variable_width = true;
	  maxWidth = INF;
	} else {
	  maxWidth = _ddWidth;
	  variable_width = false;
	}
}


//
// Generate MDD relaxation. Returns lower (upper for is) bound.
//
// TODO: prune BDD nodes according to best lower bound (using perhaps some estimate?)
//
int MinBandBDD::generateRelaxation(int initial_lp) {
//TODo
	//return -1;
	// TODO: variable width?
  //	maxWidth = MAX(10, 2000.0 / active_vertices.size());
//	cout << "max width = " << maxWidth << endl;

//	cout << "[BDD - place = " << x10_placeID << "] Node to explore:" << endl;
//	cout << "\tLP = " << initial_lp << endl;
//	cout << "\tState = { ";
//	for (vector<int>::iterator it = active_vertices.begin(); it != active_vertices.end(); ++it) {
//		cout << *it << " ";
//	}
//	cout << "}" << endl << endl;


	// ---------------------------------------------------------------------
	// 1. Initialization
	// ---------------------------------------------------------------------

	if (variable_width){
		maxWidth = INF;
	}

	lower_bound = initial_lp;
	isUBUpdated = false;

  //store vertex in layer data to restore when done with the relaxation.
  	vector<int> orig_vertex_in_layer(vertex_in_layer.begin(), vertex_in_layer.end());

	// Initialize internal structures for vertex ordering
	//vertex_in_layer.clear();
	if (var_ordering->order_type == MinState) {
		for (vector<int>::iterator v = active_vertices.begin(); v != active_vertices.end(); ++v) {
			in_state_counter[*v] = 1;
		}
	} else if ( var_ordering->order_type == RootOrder ) {
		ComparatorAuxIntVectorDescending comp(root_ordering);
		sort( active_vertices.begin(), active_vertices.end(), comp );

	} else if( var_ordering->order_type == LexOrder ) {
		sort( active_vertices.begin(), active_vertices.end() );

	} else if( var_ordering->order_type == SpanningTree ) {
		//do nothing
	} else {
		cout << "Order undefined" << endl;
		exit(0);
	}

	cout<< "Starting relaxation :: ("<< lower_bound << ") :: ["<< maxWidth << "] ::";
	// ---------------------------------------------------------------------
	// 2. Relaxation
	// ---------------------------------------------------------------------

	// Initialize pool for root node

	// create root node state

	//imake sure we have the right number of things in vertex_in_layer
	// we usually have one extra set in the state
	while (vertex_in_layer.size() >= active_state.size()and vertex_in_layer.size() >0)
	{
		vertex_in_layer.pop_back();
	}

	//if this is the very first node active_state is empty so we have to put something in there before we copy it.
	if (active_state.size() < 1 ){
		Domain full_domain;
		for (int i=0;i <= (int)(std::ceil(inst->graph->n_vertices/2)); i++) full_domain.insert(i);
		active_state.push_back(full_domain);
	}

	//copy active_state to root_state
	State root_state;
	for( std::vector<set<int> >::const_iterator i = active_state.begin(); i != active_state.end(); ++i){
		Domain domain;
		for( std::set<int>::const_iterator j = (*i).begin(); j != (*i).end(); ++j){
			//std::cout << *j ;
			domain.insert(*j);
		}
		root_state.push_back(domain);
		//cout << ",";
	}

	//print the initial state just to check
	/*for( std::vector<set<int> >::const_iterator i = active_state.begin(); i != active_state.end(); ++i){
		for( std::set<int>::const_iterator j = (*i).begin(); j != (*i).end(); ++j){
			std::cout << *j ;
		}
		cout << ",";
	}*/

	//cout<< "creating initial root : activestate/vil size " << active_state.size() << vertex_in_layer.size()<<endl;

	// create initial BDD node
	Node* initial_node = new Node(root_state, 0, true);  initial_node->printState();
	node_pool.clear();
	node_pool[ &(root_state) ] = initial_node;
	root_node = initial_node;

	BDDNodePool::iterator node_it, existing_node_it;
	Node* node;

	// relaxation control variables
	int current_vertex;
	int layer = root_state.size() -1 ;
	// so for the first vertex we have a 0, last layer n-1, corresponds to indices in vertex in layer
	const int num_active_vertices = inst->graph->n_vertices - layer;
	//here active vertices are vertices that dont have domains yet

	//Todo check count here
	while ( layer < inst->graph->n_vertices ) {
		//cout << "layer " << layer << " - ";
		// ==================================
		// 1. Vertex and BDD node selection
		// ==================================

		// select next vertex and update active vertex list
		if (var_ordering->order_type == MinState) {
			current_vertex = choose_next_vertex_min_size_next_layer();
		} else if ( var_ordering->order_type == SpanningTree ) {
			current_vertex = var_ordering->vertex_in_layer(new BDD(), layer);
			active_vertices.erase(std::remove(active_vertices.begin(), active_vertices.end(), current_vertex), active_vertices.end());
		}else if ( var_ordering->order_type == RootOrder ) {
			current_vertex = active_vertices.back();
			active_vertices.pop_back();
		} else if (var_ordering->order_type == LexOrder) {
			current_vertex = active_vertices.back();
			active_vertices.pop_back();

		} else {
			exit(0);
		}
		//cout << "Size of vertex_in_layer" << vertex_in_layer.size()<< endl;
		vertex_in_layer.push_back( current_vertex );
		assert( current_vertex != -1 );

		nodes_layer.clear();
		//cout << "Nodes_pool size " << node_pool.size();

		node_it = node_pool.begin();
		int min_cost_pool = inst->graph->n_vertices;
		while (node_it != node_pool.end())	{

			//if (node_it->second->state[current_vertex]) {
			//we need to branch on everything

			if( node_it->second->filterDomains() >= 0 ){
				// First case: state contains current vertex, thus we need to branch on it
				// if var ordering is min-in-state, decrement counter of active vertex list
				//Todo somthing like this, count tedge or something
				/*if (var_ordering->order_type == MinState) {
					for (vector<int>::iterator v = active_vertices.begin(); v != active_vertices.end(); ++v) {
						if ((*node_it->first)[*v]) {
							in_state_counter[*v]--;
						}
					}
				}*/
				//cout << "nodes_layer insert ";
				//node_it->second->printState();
				// move from the pool to current layer list
				nodes_layer.push_back(node_it->second);
				if (min_cost_pool > node_it->second->cost)
					min_cost_pool = node_it->second->cost;
				//node_pool.erase(node_it++);

			} else {
				cout << "filtered infeasible ";
				node_it->second->printState(); cout << endl;

				// Second case: state does not contain vertex; we can maintain it in the pool
				//++node_it;
			}
			node_pool.erase(node_it++);

		}

		lower_bound = MAX(lower_bound, min_cost_pool);

		//PRINT LAYER
		/*cout << "Layer " << layer << " - current vertex: " << current_vertex;
		cout << " - pool size: " << node_pool.size();
		cout << " - before merge: " << nodes_layer.size();
		cout << " - total: " << node_pool.size() + nodes_layer.size();
		cout << endl;*/


		// ==================================
		// 2. Node merging
		// ==================================
		//cout<< "width of layer " << (int)nodes_layer.size() << " max " << maxWidth;
		if( maxWidth != INF && (int)nodes_layer.size() > maxWidth ) {
			mergeLayer(layer, nodes_layer);
		}


		// ==================================
		// 3. Branching
		// ==================================

		Node* branch_node;
		for (vector<Node*>::iterator it = nodes_layer.begin(); it != nodes_layer.end(); ++it) {
			Domain full_domain;
			for (int i=0;i<inst->graph->n_vertices; i++) full_domain.insert(i);

			branch_node = (*it);

			//cout << "Branching on node "  ;
			//branch_node->printState();

			//domains should already be filtered once we get  here
			//remove most recent domain, we are now branching on it
			set<int> branch_domain = branch_node->state.back();
			branch_node->state.pop_back();


			for (set<int>::const_iterator v = branch_domain.begin(); v != branch_domain.end(); ++v){
				Domain domain;
				domain.insert(*v);

				node = new Node(branch_node->state,
									branch_node->cost ,
									branch_node->exact);

				//the one we are branching on
				node->state.push_back(domain);

				//cout << "Filtering..";
				if (node->filterDomains() < 0){
					cout << "infeasible";
					break;
				}
				//cout << "feasible. ";
				//cout << "Calculated cost" << calculateCost(node);
				//TODO filter, check if feasible domains,check costs and update

				//cout << "c1"<< calculateCost(node);
				node->cost = MAX(node->cost, calculateCost(node));
				//cout << "New node "  ;
				//node->printState();
				//cout<< endl;

				//the next variable, will be filtered once we remove from the node_pool
				node->state.push_back(full_domain);


				//rather filter domains after we get currentvertex

				if (node->cost < upper_bound) {

					// Equivalence test: check if node is in list
					existing_node_it = node_pool.find( &(node->state) );
					if (existing_node_it != node_pool.end()) {

						cout <<" Matched, merging";

						// node already exists in the pool: update node match
						merge(node, existing_node_it->second);
						node = existing_node_it->second;

					} else {

						//cout << "inserting" << endl;

						node_pool[ &(node->state)] = node;


						// update eligibility list
						/*if (var_ordering->order_type == MinState) {
							for (vector<int>::iterator v = active_vertices.begin(); v != active_vertices.end(); ++v) {
								if ((node->state)[*v]) {
									in_state_counter[*v]++;
								}
							}
						}*/
					}
				} else {
					delete node;
				}

			}
			//for every value in the current vertex's domain (the last domain) do
			// replace that domain with a singleon, do alldiff filtering
			//create a new node with filtered domains, getcost (best possible cost)
			// if cost < best_ub (from an existing solution)
			// 		add to nodepool, merge if existing or just add
			//for(vector<int>::cons)

		}

		// if the number of nodes in the pool is empty, then we do not need to explore this BDD further
		// since there are no better feasible solutions that can be generated from here
		if (node_pool.empty()) {
			//cout << "deleting nodes from branch_pool";

			while ((int)branch_nodes.size() > initialPosPool) {
				delete branch_nodes.back();
				branch_nodes.pop_back();
			}

			// clean branching pool
			for (vector<BranchNode*>::iterator it = branch_nodes.begin(); it != branch_nodes.end(); ++it) {
				delete (*it);
			}
			branch_nodes.clear();

			// reset internal parameters
			isExact = false;

			// no new lower bound was generated
			vertex_in_layer = orig_vertex_in_layer;
			return 0;
		}

		// go to next layer
		layer++;
	}

	// take info from terminal node
	assert( node_pool.size() > 0 );

	//TODO we have to do the terminal node slightly differently
	/*const Node* terminal = node_pool.begin()->second;

	upper_bound = terminal->cost;
	isExact = terminal->exact;

	delete terminal;*/




	int min_cost_pool = inst->graph->n_vertices;
	int min_exact_cost = inst->graph->n_vertices;
	//cout<< node_pool.size() << " nodes in node_pool after run"<< endl;

	node_it = node_pool.begin();
	bool hasExact = false;
	while (node_it != node_pool.end())	{

			//node_it->second->printState();

			if (node_it->second->exact){
				hasExact = true;
				if (node_it->second->cost < upper_bound	){
					upper_bound = node_it->second->cost;
					//best_ub_state = State(node_it->second->state);
				}
			}

			if (min_cost_pool > node_it->second->cost)
				min_cost_pool = node_it->second->cost;

			node_it++;
	}
	lower_bound = MAX(lower_bound, min_cost_pool);

	//if one of the equally lowest cost nodes are exact then
	//no way we can improve (due to the nature of the cost function)
	//so updatelocalUB, and we can delete the nodes in branch_nodes??
	//if (best_ub == min_cost_pool)
	//isExact = true;

	// if last node is exact, BDD is exact: update lower bound
	if (hasExact) {
		updateLocalUB(upper_bound, true);
	}

	//TODO update branch node pool: we get the lower bound from the node pool, be careful about what you update
	for (vector<BranchNode*>::iterator st = branch_nodes.begin(); st != branch_nodes.end(); ++st) {
		(*st)->relax_lb = MAX(lower_bound, (*st)->cost);
		(*st)->cost = MAX((*st)->cost, min_cost_pool);
			//(*st)->relax_lb = MAX(	lower_bound,(*st)->cost + (*st)->state.size());
	}

	/*cout <<endl<< "(In relax) Branch_nodes size after relaxation "<< branch_nodes.size() << endl;
	for (vector<BranchNode*>::iterator it = branch_nodes.begin(); it != branch_nodes.end(); ++it) {
		(*it)->printState();
	}*/

	/*if (branch_nodes.size() > 0) {
		// update upper bound of the branch nodes
		//cout << "\n[BDD " << x10_placeID << "] Branching nodes: " << endl;
		for (vector<BranchNode*>::iterator it = branch_nodes.begin(); it != branch_nodes.end(); ++it) {
			cout << "\tUB = " << (*it)->relax_ub;
			cout << " - LP = " << (*it)->cost;
			cout << " - State = { ";
			for (vector<int>::iterator v = (*it)->state.begin(); v != (*it)->state.end(); ++v) {
				cout << *v << " ";
			}
			cout << "}" << endl;
		}
	}*/

	//best_lb = MAX(best_lb, lower_bound);

	//vertex_in_layer = orig_vertex_in_layer;
	return lower_bound;
}


//
// Merge nodes in a layer to meet maximum allowed width
//
void MinBandBDD::mergeLayer(int layer, vector<Node*> &nodes_layer) {
	//cout << "Merging layer "<< endl;
	sort(nodes_layer.begin(), nodes_layer.end(), CompareNodesCost());

	// All nodes will be merged to the one with index "width-1", which will
	// be denoted here by central node
	Node* central_node = nodes_layer[maxWidth-1];

	// If central node is exact, we add it to the branch node pool
	// (since it will be maintained in the BDD, we have to make a copy)
	if( central_node->exact ) {
		addBranchNode(central_node);
		central_node->exact = false;
	}


	/*cout << endl << "before merge all nodes sorted"<< endl;
	for (vector<Node*>::iterator node = nodes_layer.begin() ; node != nodes_layer.end(); ++node) {

			(*node)->printState();
		}
	cout<< "finished allnodes"<< endl;*/


	//cout<< endl << "merging :"<< endl;
	// Merge nodes into the central node state (or simply central state)
	// no need to consider cost, central node already has min cost in ordered list

	State* central_state = &( central_node->state );
	for (vector<Node*>::iterator node = nodes_layer.begin() + maxWidth; node != nodes_layer.end(); ++node) {

		//(*node)->printState();
		// update state

		//check if two the nodes you want to merge are on the same layer
		if ((*node)->state.size() != (*central_state).size()){
			cout<< "Invalid Merge operation attempted" << endl;
			assert((*node)->state.size() != (*central_state).size());
		}

		//Take union of domains
		for( int i = 0; i < (*central_state).size(); i++){
			//no need to attempt a merge if ew already have everything
			if ((*central_state)[i].size() != inst->graph->n_vertices)
				(*central_state)[i].insert((*node)->state[i].begin(), (*node)->state[i].end());
		}

		// if any of the remaining nodes to merge is exact, we have to add it to the branch pool as well.
		if ((*node)->exact) {
			addBranchNode((*node));
		}

		// delete node
		delete (*node);
	}
	//cout<< "finished merging"<< endl;

	// resize node layer vector
	nodes_layer.resize(maxWidth);

	/*cout << endl << "begin all nodes after merge"<< endl;
	for (vector<Node*>::iterator node = nodes_layer.begin() ; node != nodes_layer.end(); ++node) {
		(*node)->printState();
	}
	cout<< "finished allnodes after merge"<< endl;*/

	// Check if there are any other nodes which are equivalent to the central node

	Node* node;
	for (int i = 0; i <= maxWidth-2; ++i) {
		node = nodes_layer[i];

		// check if this state already exists in layer nodes
		if( node->state == (*central_state) ) {

			// notice that, if node exists, we do not need to update
			// the costs because the vector is already ordered

			// if existing node is exact, we have to add it to the branch pool
			// (since it will be maintained in the BDD, we have to make a copy)
			if (node->exact) {
				addBranchNode(node);
				node->exact = false;
			}

			// delete the last node (i.e., the central one)
			delete nodes_layer.back();

			// remove it from queue
			nodes_layer.pop_back();
			break;
		}
	}
}


//
// From min-in-state ordering: choose vertex that participates in the least number of states
// and remove it from the active vertex list
//
int MinBandBDD::choose_next_vertex_min_size_next_layer() {
	return -1;
	int sel_index = 0;
	for (size_t i = 1; i < active_vertices.size(); ++i) {

		if (in_state_counter[active_vertices[i]] != 0 && in_state_counter[active_vertices[i]] < in_state_counter[active_vertices[sel_index]]) {
			sel_index = i;

		} else if (in_state_counter[active_vertices[i]] != 0 && (in_state_counter[active_vertices[i]] == in_state_counter[active_vertices[sel_index]]) && active_vertices[i] < active_vertices[sel_index]) {
			// lexicographic tie breaking
			sel_index = i;
		}
	}

	// TODO: check this case!!!
	if (in_state_counter[active_vertices[sel_index]] == 0) {

//		cout << "[BDD " << x10_placeID << "] POOL SIZE = " << node_pool.size() << endl;
//		cout << "State = { ";
//		for (int v = (*(node_pool.begin()->first)).find_first(); v != state_end; v = (*(node_pool.begin()->first)).find_next(v) ) {
//			cout << v << " ";
//		}
//		cout << "}" << endl;
//
//		cout << "Active vertices = { ";
//		for (vector<int>::iterator it = active_vertices.begin(); it != active_vertices.end(); ++it) {
//			cout << *it << " (counter = " << in_state_counter[*it] << ")";
//		}
//		cout << "}" << endl;
//
//		exit(1);
	}

	//assert( in_state_counter[active_vertices[sel_index]] > 0 );

	/*int sel_vertex = active_vertices[sel_index];
	active_vertices[sel_index] = active_vertices.back();
	active_vertices.pop_back();
	return sel_vertex;*/
}


//
// Generate MDD restriction. Returns lower bound.
//
int MinBandBDD::generateRestriction(const int initial_lp) {

    // ---------------------------------------
	// 1. Initialization
	// ---------------------------------------

	// width is the number of vertices
  if (variable_width) {
    maxWidth = active_vertices.size();
  }

	// ---------------------------------------------------------------------
	// 1. Initialization
	// ---------------------------------------------------------------------

	// set width

  //store vertex in layer data to restore when done with the relaxation.
	vector<int> orig_vertex_in_layer(vertex_in_layer.begin(), vertex_in_layer.end());
	// Initialize internal structures for vertex ordering
	//vertex_in_layer.clear();

	if (var_ordering->order_type == MinState) {
		for (vector<int>::iterator v = active_vertices.begin(); v != active_vertices.end(); ++v) {
			in_state_counter[*v] = 1;
		}
	} else if ( var_ordering->order_type == RootOrder ) {
		ComparatorAuxIntVectorDescending comp(root_ordering);
		sort( active_vertices.begin(), active_vertices.end(), comp );

	} else if( var_ordering->order_type == LexOrder ) {
		sort( active_vertices.begin(), active_vertices.end() );

	} else if( var_ordering->order_type == SpanningTree ) {
		//do nothing
	} else {
		cout << "Order undefined" << endl;
		exit(0);
	}

	cout<< "Starting restriction ";
	// ---------------------------------------------------------------------
	// 2. Restriction
	// ---------------------------------------------------------------------

	// Initialize pool for root node

	// create root node state

	//pop vertex in layres, will get again in the loop
	while (vertex_in_layer.size() >= active_state.size() and vertex_in_layer.size() >0)
		{
			vertex_in_layer.pop_back();
		}

	//if this is the very first node active_state is empty so we have to put something in there before we copy it.
	if (active_state.size() < 1 ){
		Domain full_domain;
		for (int i=0; i<= (int)(std::ceil(inst->graph->n_vertices/2)); i++) full_domain.insert(i);
		active_state.push_back(full_domain);
	}

	//copy active_state to root_state
	State root_state;
	for( std::vector<set<int> >::const_iterator i = active_state.begin(); i != active_state.end(); ++i){
		Domain domain;
		for( std::set<int>::const_iterator j = (*i).begin(); j != (*i).end(); ++j){
			//std::cout << *j ;
			domain.insert(*j);
		}
		root_state.push_back(domain);
		//cout << ",";
	}

	/*//print the initial state just to check
	for( std::vector<set<int> >::const_iterator i = active_state.begin(); i != active_state.end(); ++i){
		for( std::set<int>::const_iterator j = (*i).begin(); j != (*i).end(); ++j){
			std::cout << *j ;
		}
		cout << ",";
	}*/

	// create initial BDD node
	Node* initial_node = new Node(root_state, 0, true); initial_node->printState();
	node_pool.clear();
	node_pool[ &(root_state) ] = initial_node;
	root_node = initial_node;

	BDDNodePool::iterator node_it, existing_node_it;
	Node* node;

	// relaxation control variables
	int current_vertex;
	int layer = root_state.size() -1 ;
	// so for the first vertex we have a 0, last layer n-1, corresponds to indices in vertex in layer
	const int num_active_vertices = inst->graph->n_vertices - layer;
	//here active vertices are vertices that dont have domains yet

	//Todo check count here
	while ( layer < inst->graph->n_vertices ) {
		//cout << "layer " << layer << " - ";
		// ==================================
		// 1. Vertex and BDD node selection
		// ==================================

		// select next vertex and update active vertex list
		if (var_ordering->order_type == MinState) {
			current_vertex = choose_next_vertex_min_size_next_layer();
		} else if ( var_ordering->order_type == SpanningTree ) {
			current_vertex = var_ordering->vertex_in_layer(new BDD(), layer);
			active_vertices.erase(std::remove(active_vertices.begin(), active_vertices.end(), current_vertex), active_vertices.end());
		}else if ( var_ordering->order_type == RootOrder ) {
			current_vertex = active_vertices.back();
			active_vertices.pop_back();
		} else if (var_ordering->order_type == LexOrder) {
			current_vertex = active_vertices.back();
			active_vertices.pop_back();

		} else {
			exit(0);
		}

		vertex_in_layer.push_back( current_vertex );
		assert( current_vertex != -1 );

		nodes_layer.clear();
		//cout << "Nodes_pool size " << node_pool.size();

		node_it = node_pool.begin();
		int min_cost_pool = inst->graph->n_vertices;
		while (node_it != node_pool.end())	{

			//this code is moved out of the if to try and filter only once
			nodes_layer.push_back(node_it->second);
			if (min_cost_pool > node_it->second->cost)
				min_cost_pool = node_it->second->cost;

			/* hopefully we can filter one time fewer
			 * if( node_it->second->filterDomains() >= 0 ){

				//node_it->second->printState();
				// move from the pool to current layer list
				nodes_layer.push_back(node_it->second);
				if (min_cost_pool > node_it->second->cost)
					min_cost_pool = node_it->second->cost;

			} else {
				cout << "filtered infeasible ";
				node_it->second->printState(); cout << endl;

			}*/
			node_pool.erase(node_it++);

		}

		//lower_bound = MAX(lower_bound, min_cost_pool);

		//PRINT LAYER
		/*cout << "Layer " << layer << " - current vertex: " << current_vertex;
		cout << " - pool size: " << node_pool.size();
		cout << " - before merge: " << nodes_layer.size();
		cout << " - total: " << node_pool.size() + nodes_layer.size();
		cout << endl;*/


		// ==================================
		// 2. Node merging
		// ==================================
		//cout<< "width of layer " << (int)nodes_layer.size() << " max " << maxWidth;
		if( maxWidth != INF && (int)nodes_layer.size() > maxWidth ) {
			restrictLayer(layer, nodes_layer);
		}


		// ==================================
		// 3. Branching
		// ==================================

		Node* branch_node;
		Domain full_domain;
		//Domain domain;
		for (int i=0;i<inst->graph->n_vertices; i++) full_domain.insert(i);

		for (vector<Node*>::iterator it = nodes_layer.begin(); it != nodes_layer.end(); ++it) {

			branch_node = (*it);

			//cout << "Branching on node "  ;
			//branch_node->printState();

			//domains should already be filtered once we get  here
			//remove most recent domain, we are now branching on it
			set<int> branch_domain = branch_node->state.back();
			//branch_node->state.pop_back();
			branch_node->state.back().clear();


			//add the next domain so long, now you are branching on second to last domain
			int branch_dom_number = branch_node->state.size();
			Domain domain(full_domain);
			branch_node->state.push_back(domain);



			for (set<int>::const_iterator v = branch_domain.begin(); v != branch_domain.end(); ++v){
				//Domain domain(full_domain);
				//domain.insert(*v);

				//insert, create new, remove to insert next
				//branch_node->printState(); cout <<"Inserting into " << branch_dom_number << ", " << *v << "\n";
				branch_node->state[branch_dom_number-1].insert(*v);
				//branch_node->printState(); cout << "After]\n";
				node = new Node(branch_node->state,
						branch_node->cost ,
						branch_node->exact);

				//remove
				branch_node->state[branch_dom_number-1].clear();

				//node->state.push_back(domain);

				//cout << "Filtering..";
				if ( node->filterDomains2() >= 0 ){


					//cout << "feasible. ";
					//cout << "Calculated cost" << calculateCost(node);
					//TODO filter, check if feasible domains,check costs and update

					//cout << "c1"<< calculateCost(node);
					node->cost = MAX(node->cost, calculateCost(node));
					//cout << "New node "  ;
					//node->printState();
					//cout<< endl;

					//the next variable, will be filtered once we remove from the node_pool
					//node->state.push_back(full_domain);


					//rather filter domains after we get currentvertex

					if (node->cost < best_ub) {
						//cout << "Cheap enough: " ; node->printState();
						// Equivalence test: check if node is in list
						existing_node_it = node_pool.find( &(node->state) );
						if (existing_node_it != node_pool.end()) {

							//cout <<" Matched, merging";

							// node already exists in the pool: update node match
							merge(node, existing_node_it->second);
							node = existing_node_it->second;

						} else {

							//cout << "inserting" << endl;
							node_pool[ &(node->state)] = node;

						}
					} else {
						//cout << "Too expensive: "; node->printState();
						delete node;
					}
				}
				else{
					delete node;
				}

			}
			//for every value in the current vertex's domain (the last domain) do
			// replace that domain with a singleon, do alldiff filtering
			//create a new node with filtered domains, getcost (best possible cost)
			// if cost < best_ub (from an existing solution)
			// 		add to nodepool, merge if existing or just add
			//for(vector<int>::cons)

		}

		// go to next layer
		layer++;
	}


	int ub = inst->graph->n_vertices;
	//cout<< node_pool.size() << " nodes in node_pool after run"<< endl;

	node_it = node_pool.begin();
	while (node_it != node_pool.end())	{

		// 	node_it->second->printState();

		if (node_it->second->exact){
			//cout << "exact_restrict_node="<< node_it->second->cost << " ";
			if (node_it->second->cost < ub	){
				ub = node_it->second->cost;
				//best_ub_state = State(node_it->second->state);
			}
		}
		node_it++;
	}
	cout  << "restriction ub  = " 	<< ub << endl;

	best_ub = MIN(best_ub, ub);

	//vertex_in_layer = orig_vertex_in_layer;
	return ub;
 }

//
// Reduce the size of the layer
//
void MinBandBDD::restrictLayer(int layer, vector<Node*> &nodes_layer) {
	sort(nodes_layer.begin(), nodes_layer.end(), CompareNodesCost());
	for( vector<Node*>::iterator node = nodes_layer.begin()+maxWidth;
			node != nodes_layer.end(); ++node)
	{
		delete (*node);
	}
	nodes_layer.resize(maxWidth);
}

int MinBandBDD::calculateCost(Node* _node){

	//TODO different way of calculating cost: only singletons, maybe even only newest vertex

	//cout<< "calculating costs" << endl;
	//cout << "calculating costs on ";
	// _node->printState() ;

	//copy the state (for some reason)
	State _state (_node->state)	;

	// we need to have at least two domains to calculate a cost.
	if (vertex_in_layer.size() < 2 || _node->state.size()<2)
		return 0;


	vector<vector<int> > edges_to_check;
	//TODO assume ordering doesnt change, this can be done once globally
	int vertex_counter = 0;
	for (vector<int>::iterator i = vertex_in_layer.begin(); i != vertex_in_layer.end()-1 ; ++i){
		vector<int>  edges_renamed;
		for (vector<int>::iterator j = i+1; j!=vertex_in_layer.end(); ++j	){
			if (inst->graph->adj_m[*i][*j]){
				//the edge we are pushing is the position of j in vertex_in_layer
				edges_renamed.push_back( j - vertex_in_layer.begin() );
			}
			else{
			}
		}
		edges_to_check.push_back(edges_renamed);
	}


	/*//print
	for( std::vector<int>::const_iterator j = vertex_in_layer.begin(); j != vertex_in_layer.end(); ++j){
		std::cout << *j ;
	}
	cout << endl;*/
	/*cout << "renamed edges to check: ";
	for( std::vector<vector<int> >::const_iterator i = edges_to_check.begin(); i != edges_to_check.end(); ++i){
		//cout << *i <<" - ";
		for( std::vector<int>::const_iterator j = (*i).begin(); j != (*i).end(); ++j){
			std::cout << *j ;
		}
		cout << ";";
	}*/ /*
	cout << endl;
	cout << _state.size()<< vertex_in_layer.size()<< edges_to_check.size(); */

	int largest_smallest_cost = 0;

	//for (int i = 0; i < MIN(_state.size(), vertex_in_layer.size())-1; i++){
	for (int i = 0; i < edges_to_check.size(); i++){
		//cout << "domain1= "<<i<< " - " << edges_to_check[i].size() << " :" << (edges_to_check[i].begin() -edges_to_check[i].end()) ;
		for (vector<int>::const_iterator j = edges_to_check[i].begin(); j != edges_to_check[i].end(); ++j){
			//cout << "enter"<< i<< *j << ";"<<_state[i].size()<< _state[*j].size() << endl;
			int smallest_cost_edge = inst->graph->n_vertices+1;
			bool smallEnough = false;

			//handle case where both domains are everything
			if (_state[i].size() == inst->graph->n_vertices and _state[*j].size() == inst->graph->n_vertices){
				smallest_cost_edge =1;
			}
			else{

				//TODO only upper and lower bounds check distance
				for (set<int>::const_iterator i_val = _state[i].begin(); i_val != _state[i].end() and !smallEnough; ++i_val){
					//cout << "x";
					//test if ew are already had longer on any of the edges, if so this will only get smaller, bail
					//test if we are close enough to i_val that we ca decrease smallest_edge, if not no point in checking
					//for (set<int>::const_iterator j_val = _state[*j].begin(); j_val != _state[*j].end() and !smallEnough and ((*j_val) <= (*i_val) + smallest_cost_edge); ++j_val){
					for (set<int>::const_iterator j_val = _state[*j].begin(); j_val != _state[*j].end() and !smallEnough; ++j_val){
						//int a = (*i_val);
						//int b = (*j_val);
						//cout<< "testing "<<a << b<<endl;

						if ((*i_val) != (*j_val) and (std::abs((*i_val) - (*j_val)) < smallest_cost_edge)){
							smallest_cost_edge = std::abs((*i_val) - (*j_val));
							if (smallest_cost_edge < largest_smallest_cost)
								smallEnough = true;
						}


					}
				}
			}

			//cout << smallest_cost_edge << " ";
			if (smallest_cost_edge > largest_smallest_cost
					and smallest_cost_edge <= inst->graph->n_vertices){
				largest_smallest_cost = smallest_cost_edge;
			}
			//cout << "exit" <<endl;
		}
		//cout << endl;
	}

	//cout << "Done "<< largest_smallest_cost << endl;
	return largest_smallest_cost;

}

