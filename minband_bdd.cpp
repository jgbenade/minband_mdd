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
	//edges_to_check = new vector<vector<int> >;

	//in_node->state_counter = new int[inst->graph->n_vertices];
	//inst->graph->print();


	// ============================================================
	// Generate initial relaxation to derive variable ordering
	// TODO: save initial UB
	// ============================================================

	cout << "[BDD] Root node computation..." << endl;

	lb.resize(inst->graph->n_vertices,-1);
	ub.resize(inst->graph->n_vertices,-1);

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
	var_ordering = new FrontBackOrdering(inst);

	cout << " ### Degree bound\t" << inst->calculate_Degree_Bound() << endl;
	cout << " ### Half-density bound\t " << inst->calculate_HalfDensity_Bound() << endl;
	cout << " ### Half-density bound better\t " << inst->calculate_HalfDensity_Bound_better() << endl;
	cout << " ### Caprara bound   \t " << inst->calculate_Caprara_Bound() << endl;

    best_lb = std::max(inst->caprara_bound, inst->calculate_HalfDensity_Bound_better());
    cout << "best_lb = " << best_lb << "\n";

    // generate relaxation

	active_vertices.resize(inst->graph->n_vertices);
	for (myint i = 0; i < inst->graph->n_vertices; ++i) {
		active_vertices[i] = i;
	}

	best_ub = inst->graph->n_vertices;
	upper_bound = inst->graph->n_vertices;

	lower_bound = -1;
	int ub1 = inst->CuthillMckee();
	cout << "Cuthill-Mckee: " << ub1 <<endl;

	int ub2 = generateRestriction(inst->graph->n_vertices);
	cout << "restriction "<< ub2 << endl;
	upper_bound = MIN(upper_bound,ub2);
	//upper_bound = 11;
	int lb1 = generateFakeRelaxation(-1);
	best_lb = MAX(lb1, best_lb);

	cout << "returned " << lb1 << best_ub<<endl;
	//int ub1 = generateRestriction(inst->graph->n_vertices);

	cout << " ### Lower bound: " << lower_bound << endl;
	cout << " ### Upper bound: " << upper_bound << endl;

	// check if BDD is already exact
	//Not valid for fake relaxation
	/*isExact = branch_nodes.empty();
	if (isExact) {
		return;
	}*/

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
//todo variablewidthd
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
	} else if( var_ordering->order_type == Degree ) {
		//do nothing
	} else if( var_ordering->order_type == FromFront ) {
		//do nothing
	} else if( var_ordering->order_type == FrontBack ) {
		//do nothing
	}
	else {
		cout << "Order undefined" << endl;
		exit(0);
	}

	//cout<< "Starting relaxation :: ("<< lower_bound << ") :: ["<< maxWidth << "] ::";
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
	/*if (active_state.size() < 1 ){
		vector<myint> full_domain;
		for (myint i=0;i <= (int)(std::ceil(inst->graph->n_vertices/2)); i++) full_domain.insert(i);
		active_state.push_back(full_domain);
	}*/

	//initialise all feasible root state //TODO copy active state
	State root_state;
	Domain fulldomain;
	if (active_state.size() < fulldomain.size()	){
		for(int i =0; i<fulldomain.size(); i++){
			Domain domain = ~fulldomain;
			root_state.push_back(domain);
		}
	}
	else{
		root_state = active_state;
	}



	//print the initial state just to check
	/*for( std::vector<set<myint> >::const_iterator i = active_state.begin(); i != active_state.end(); ++i){
		for( std::set<myint>::const_iterator j = (*i).begin(); j != (*i).end(); ++j){
			std::cout << *j ;
		}
		cout << ",";
	}*/

	//cout<< "creating initial root : activestate/vil size " << active_state.size() << vertex_in_layer.size()<<endl;

	// create initial BDD node
	Node* initial_node = new Node(root_state, 0, true); // initial_node->printState();
	node_pool.clear();
	node_pool[ &(root_state) ] = initial_node;
	root_node = initial_node;

	BDDNodePool::iterator node_it, existing_node_it;
	Node* node;
	bool exactNodeInPool = true; // as soon as we have no more exact nodes left we can jump out asd start checking out children

	// relaxation control variables
	int current_vertex;
	//TODO this only allows for root, figure something else out.
	//maybe give branchnode a layer, copy info over when creating rootnode.
	int layer = 0 ;
	// so for the first vertex we have a 0, last layer n-1, corresponds to indices in vertex in layer
	const int num_active_vertices = inst->graph->n_vertices - layer;
	//here active vertices are vertices that dont have domains yet

	//Todo enable/disable early breaks
	while ( layer < inst->graph->n_vertices/* and exactNodeInPool*/ ) {
		//cout << "layer " << layer << " - ";
		exactNodeInPool = false;
		// ==================================
		// 1. Vertex and BDD node selection
		// ==================================

		// select next vertex and update active vertex list
		if (var_ordering->order_type == MinState) {
			current_vertex = choose_next_vertex_min_size_next_layer();
		} else if ( var_ordering->order_type == SpanningTree ) {
			current_vertex = var_ordering->vertex_in_layer(new BDD(), layer);
			active_vertices.erase(std::remove(active_vertices.begin(), active_vertices.end(), current_vertex), active_vertices.end());
		}else if ( var_ordering->order_type == Degree ) {
			BDD* temp = new BDD();
			current_vertex = var_ordering->vertex_in_layer(temp, layer);
			active_vertices.erase(std::remove(active_vertices.begin(), active_vertices.end(), current_vertex), active_vertices.end());
			delete temp;
		}else if ( var_ordering->order_type == FrontBack ) {
			BDD* temp = new BDD();
			current_vertex = var_ordering->vertex_in_layer(temp, layer);
			active_vertices.erase(std::remove(active_vertices.begin(), active_vertices.end(), current_vertex), active_vertices.end());
			delete temp;
		}else if ( var_ordering->order_type == FromFront ) {
			BDD* temp = new BDD();
			current_vertex = var_ordering->vertex_in_layer(temp, layer);
			active_vertices.erase(std::remove(active_vertices.begin(), active_vertices.end(), current_vertex), active_vertices.end());
			delete temp;
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

		for (vector<Node*>::iterator nit = nodes_layer.begin(); nit != nodes_layer.end(); ++nit)	delete *nit;
		nodes_layer.clear();
		//nodes_layer.resize(node_pool.size());
		//cout << "Nodes_pool size " << node_pool.size();

		node_it = node_pool.begin();
		int min_cost_pool = inst->graph->n_vertices;

		while (node_it != node_pool.end())	{
			//nodes_layer[node_pool_counter++] = node_it->second;
			if ((node_it->second->exact)) exactNodeInPool = true;

			nodes_layer.push_back(node_it->second);
			if (min_cost_pool > node_it->second->cost)
				min_cost_pool = node_it->second->cost;

			node_pool.erase(node_it++);
		}

		lower_bound = MAX(lower_bound, min_cost_pool);

		//PRINT LAYER
		cout << "Layer " << layer << " - vertex: " << current_vertex;
		cout << " - pool: " << node_pool.size();
		cout << " - before merge: " << nodes_layer.size();
		cout << " - total: " << node_pool.size() + nodes_layer.size();

		// ==================================
		// 2. Node merging
		// ==================================
		//cout<< "width of layer " << (int)nodes_layer.size() << " max " << maxWidth;
		if( maxWidth != INF && (int)nodes_layer.size() > maxWidth ) {
			mergeLayer(layer, nodes_layer);
		}
		cout << " - bound: " << (*nodes_layer.begin())->cost <<","<<(nodes_layer[nodes_layer.size()-1])->cost ;		if (!exactNodeInPool) cout << "x"; cout << endl;
		// ==================================
		// 3. Branching
		// ==================================

		Node* branch_node;
		Domain domain;

		for (vector<Node*>::iterator it = nodes_layer.begin(); it != nodes_layer.end(); ++it) {

			branch_node = (*it);
			//cout<< "a"<<endl;

			//this only filters the domains that we are about to branch on
			filterBounds2(branch_node->state);

			//cout<<"b"<< (*(--branch_node->state.end())).size() <<endl;
			//cout << "Branching on node "  ;
			//branch_node->printState();

			//domains should already be filtered once we get  here
			//remove most recent domain, we are now branching on it
			// definition can be moved outside, seems to be slower??
			Domain branch_domain = branch_node->state[layer];
			//cout << branch_domain<< endl;
			for (int v = 0; v< inst->graph->n_vertices; v++){
				if (branch_domain[v]){
					//domain.reset();
					//domain.set(v);

					node = new Node(branch_node->state,
							branch_node->cost ,
							branch_node->exact);

					//the one we are branching on
					//cout << domain << "," << node->state[layer];
					//node->state[layer] &= domain;
					node->state[layer].reset();
					node->state[layer].set(v);
					//cout << "," << node->state[layer];

					if (node->filterDomains5(layer) < 0|| filterBounds2(node->state) < 0 ){
						//if (  node->filterDomains5() < 0){
						delete node;
					}
					else{
						int cost=0;
						//	node->cost = MAX(node->cost, calculateCost_caprara_fixed(node));
						if (vertex_in_layer[vertex_in_layer.size()-1]== 0 or vertex_in_layer[vertex_in_layer.size()-1]== inst->graph->n_vertices-1)
							node->cost = MAX(node->cost, inst->caprara_list[v]);
						else
							cost = calculateCost_caprara_pos(node,v);
						node->cost = MAX(node->cost, cost);

						if (node->cost < upper_bound)
							cost = calculateCost_mu2(node)	;
						node->cost = MAX(node->cost,cost);

						/*if (node->exact and node->cost < upper_bound)
							cost = calculateCost_ILP2(node);
						node->cost = MAX(node->cost, cost );*/

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
							}
						} else { // cost > upperbounde
							delete node;
						}
					}
				}
			}
			//for every value in the current vertex's domain (the last domain) do
			// replace that domain with a singleon, do alldiff filtering
			//create a new node with filtered domains, getcost (best possible cost)
			// if cost < best_ub (from an existing solution)
			// 		add to nodepool, merge if existing or just add
			//for(vector<int>::cons)

		}
		//delete branch_node;

		// if the number of nodes in the pool is empty, then we do not need to explore this BDD further
		// since there are no better feasible solutions that can be generated from here
		if (node_pool.empty()) {
			cout << "deleting nodes from branch_pool";

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

	//cout << node_pool.size();
	node_it = node_pool.begin();
	bool hasExact = false;
	int count_low =0;
	int count_exact = 0;


	while (node_it != node_pool.end())	{

			//node_it->second->printState();

			if (node_it->second->exact){
				hasExact = true;
				//for the inexact cost calculations we need to double check
				node_it->second->cost = MAX(node_it->second->cost, calculateCost_bounds(node_it->second));
				//cout << "a" << node_it->second->cost << calculateCost_bounds(node_it->second);
				if (node_it->second->cost < upper_bound	){
					upper_bound = node_it->second->cost;
					//best_ub_node->state = State(node_it->second->state);

					//print ordering
					for (int pos = 0; pos < inst->graph->n_vertices; pos++){
						for  (int v = 0; v< inst->graph->n_vertices; v++)
							if (node_it->second->state[pos][v])	cout << v << ",";
					}cout << ": " << node_it->second->cost <<  endl;
					//
				}
				count_exact++;

			}

			if (node_it->second-> cost < min_cost_pool){
				min_cost_pool = MIN(min_cost_pool, node_it->second->cost);
				count_low = 1;

			}else if (node_it->second-> cost == min_cost_pool)
				count_low++;

			//cout << "n" << node_it->second->cost<< "m"<<min_cost_pool;
			/*if (min_cost_pool > node_it->second->cost){
				min_cost_pool = node_it->second->cost;
			}*/
			delete	 node_it->second;
			node_it++;
	}
	cout <<endl << "at lb: "<< count_low << " / " << node_pool.size() << "\t Exact: " << count_exact << " / " << node_pool.size()<< endl <<branch_nodes.size( ) << endl;
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
	for (vector<Node*>::iterator nit = nodes_layer.begin(); nit != nodes_layer.end(); ++nit)	delete *nit;
	nodes_layer.clear();
	node_pool.clear();

	return lower_bound;
}


//
// Generate MDD relaxation that mimic a partial enumerative procedure, so no merging takes place
// The lowest cost node on every layer is stored in a vector, this is used to compute the lower bounds
//
// TODO: prune BDD nodes according to best lower bound (using perhaps some estimate?)
//
int MinBandBDD::generateFakeRelaxation(int initial_lp) {
//todo variablewidthd
	//	maxWidth = MAX(10, 2000.0 / active_vertices.size());
//	cout << "max width = " << maxWidth << endl;

//	cout << "[BDD - place = " << x10_placeID << "] Node to explore:" << endl;
//	cout << "\tLP = " << initial_lp << endl;
//	cout << "\tState = { ";
//	for (vector<int>::iterator it = active_vertices.begin(); it != active_vertices.end(); ++it) {
//		cout << *it << " ";
//	}
//	cout << "}" << endl << endl;

	//lowest cost on every layer among the nodes that would have been merged (now just ignored)
	vector<int> lowest_cost;
	int lowest_unbranched = inst->graph->n_vertices;

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
	} else if( var_ordering->order_type == Degree ) {
		//do nothing
	} else if( var_ordering->order_type == FromFront ) {
		//do nothing
	} else if( var_ordering->order_type == FrontBack ) {
		//do nothing
	}
	else {
		cout << "Order undefined" << endl;
		exit(0);
	}

	//cout<< "Starting relaxation :: ("<< lower_bound << ") :: ["<< maxWidth << "] ::";
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
	/*if (active_state.size() < 1 ){
		vector<myint> full_domain;
		for (myint i=0;i <= (int)(std::ceil(inst->graph->n_vertices/2)); i++) full_domain.insert(i);
		active_state.push_back(full_domain);
	}*/

	//initialise all feasible root state //TODO copy active state
	State root_state;
	Domain fulldomain;
	if (active_state.size() < fulldomain.size()	){
		for(int i =0; i<fulldomain.size(); i++){
			Domain domain = ~fulldomain;
			root_state.push_back(domain);
		}
	}
	else{
		root_state = active_state;
	}

	// create initial BDD node
	Node* initial_node = new Node(root_state, 0, true); // initial_node->printState();
	node_pool.clear();
	node_pool[ &(root_state) ] = initial_node;
	root_node = initial_node;

	BDDNodePool::iterator node_it, existing_node_it;
	Node* node;
	bool exactNodeInPool = true; // as soon as we have no more exact nodes left we can jump out asd start checking out children

	// relaxation control variables
	int current_vertex;
	//maybe give branchnode a layer, copy info over when creating rootnode.
	int layer = 0 ;
	// so for the first vertex we have a 0, last layer n-1, corresponds to indices in vertex in layer
	const int num_active_vertices = inst->graph->n_vertices - layer;
	//here active vertices are vertices that dont have domains yet

	//Todo enable/disable early breaks
	while ( layer < inst->graph->n_vertices/* and exactNodeInPool*/ ) {
		//cout << "layer " << layer << " - ";
		exactNodeInPool = false;
		// ==================================
		// 1. Vertex and BDD node selection
		// ==================================

		// select next vertex and update active vertex list
		if (var_ordering->order_type == MinState) {
			current_vertex = choose_next_vertex_min_size_next_layer();
		} else if ( var_ordering->order_type == SpanningTree ) {
			current_vertex = var_ordering->vertex_in_layer(new BDD(), layer);
			active_vertices.erase(std::remove(active_vertices.begin(), active_vertices.end(), current_vertex), active_vertices.end());
		}else if ( var_ordering->order_type == Degree ) {
			BDD* temp = new BDD();
			current_vertex = var_ordering->vertex_in_layer(temp, layer);
			active_vertices.erase(std::remove(active_vertices.begin(), active_vertices.end(), current_vertex), active_vertices.end());
			delete temp;
		}else if ( var_ordering->order_type == FrontBack ) {
			BDD* temp = new BDD();
			current_vertex = var_ordering->vertex_in_layer(temp, layer);
			active_vertices.erase(std::remove(active_vertices.begin(), active_vertices.end(), current_vertex), active_vertices.end());
			delete temp;
		}else if ( var_ordering->order_type == FromFront ) {
			BDD* temp = new BDD();
			current_vertex = var_ordering->vertex_in_layer(temp, layer);
			active_vertices.erase(std::remove(active_vertices.begin(), active_vertices.end(), current_vertex), active_vertices.end());
			delete temp;
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

		for (vector<Node*>::iterator nit = nodes_layer.begin(); nit != nodes_layer.end(); ++nit)	delete *nit;
		nodes_layer.clear();

		//nodes_layer.resize(node_pool.size());
		//cout << "Nodes_pool size " << node_pool.size();

		node_it = node_pool.begin();
		int min_cost_pool = inst->graph->n_vertices;

		while (node_it != node_pool.end())	{
			//nodes_layer[node_pool_counter++] = node_it->second;
			if ((node_it->second->exact)) exactNodeInPool = true;

			nodes_layer.push_back(node_it->second);
			if (min_cost_pool > node_it->second->cost)
				min_cost_pool = node_it->second->cost;

			node_pool.erase(node_it++);
		}

		//here we look at all the nodes on the new layer, and the
		//lowest cost among all higher nodes that were not branched on
		lower_bound = MAX(lower_bound, MIN(min_cost_pool, lowest_unbranched) );


		//PRINT LAYER
		cout << "Layer " << layer << " - vertex: " << current_vertex;
		cout << " - pool: " << node_pool.size();
		cout << " - before merge: " << nodes_layer.size();
		cout << " - total: " << node_pool.size() + nodes_layer.size();

		// ==================================
		// 2. Node merging
		// ==================================
		//cout<< "width of layer " << (int)nodes_layer.size() << " max " << maxWidth;

		//instead of merging we just sort the nodes here.
		if( maxWidth != INF && (int)nodes_layer.size() > maxWidth ) {
			//mergeLayer(layer, nodes_layer);
			sort(nodes_layer.begin(), nodes_layer.end(), CompareNodesCostDelta());
			lowest_cost.push_back(nodes_layer[maxWidth]->cost);
			lowest_unbranched = MIN(lowest_unbranched, nodes_layer[maxWidth]->cost);
			isExact = false;
		}
		cout << " - bound: " << (*nodes_layer.begin())->cost <<","<<(nodes_layer[nodes_layer.size()-1])->cost ;		if (!exactNodeInPool) cout << "x"; cout << " " << lowest_unbranched;  cout << endl;
		// ==================================
		// 3. Branching
		// ==================================

		lower_bound = MAX(lower_bound, MIN(min_cost_pool, lowest_unbranched));


		Node* branch_node;
		Domain domain;

		//only on the first few, the rest would have been merged and are now ignored
		for (vector<Node*>::iterator it = nodes_layer.begin(); it != nodes_layer.end() && it != nodes_layer.begin()+ maxWidth ; ++it) {

			branch_node = (*it);


			//this only filters the domains that we are about to branch on
			filterBounds2(branch_node->state);

			//domains should already be filtered once we get  here
			//remove most recent domain, we are now branching on it
			// definition can be moved outside, seems to be slower??
			Domain branch_domain = branch_node->state[layer];
			//cout << branch_domain<< endl;
			for (int v = 0; v< inst->graph->n_vertices; v++){
				if (branch_domain[v]){
					//domain.reset();
					//domain.set(v);

					node = new Node(branch_node->state,
							branch_node->cost ,
							branch_node->exact);

					//the one we are branching on
					//cout << domain << "," << node->state[layer];
					//node->state[layer] &= domain;
					node->state[layer].reset();
					node->state[layer].set(v);
					//cout << "," << node->state[layer];

					if (node->filterDomains5(layer) < 0|| filterBounds2(node->state) < 0 ){
						//if (  node->filterDomains5() < 0){
						delete node;
					}
					else{
						int cost=0;
						//	node->cost = MAX(node->cost, calculateCost_caprara_fixed(node));
						if (vertex_in_layer[vertex_in_layer.size()-1]== 0 or vertex_in_layer[vertex_in_layer.size()-1]== inst->graph->n_vertices-1)
							node->cost = MAX(node->cost, inst->caprara_list[v]);
						else
							cost = calculateCost_caprara_pos(node,v);
						node->cost = MAX(node->cost, cost);

						if (node->cost < upper_bound)
							cost = calculateCost_mu2(node)	;
						node->cost = MAX(node->cost,cost);

						//cgen cost
						/*if (node->cost < upper_bound)
							cost = calculateCost_caprara_gen(node)	;
						node->cost = MAX(node->cost,cost);*/

						//ilp cost
						/*
						if (node->exact and node->cost < upper_bound)
							cost = calculateCost_ILP2(node);
						node->cost = MAX(node->cost, cost );*/

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
							}
						} else { // cost > upperbounde
							delete node;
						}
					}
				}
			}
		}
		// if the number of nodes in the pool is empty, then we do not need to explore this BDD further
		// since there are no better feasible solutions that can be generated from here
		if (node_pool.empty()) {
			cout << "deleting nodes from branch_pool";

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

			if (lowest_cost.size() > 0 ){
				cout << endl << lowest_unbranched << endl;
				return lowest_unbranched	;
			}
			else{
				// Exact, never ignored nodes, return upper bound
				return upper_bound;
			}
		}

		// go to next layer
		layer++;
	}

	// take info from terminal node
	assert( node_pool.size() > 0 );

	int min_cost_pool = inst->graph->n_vertices;
	int min_exact_cost = inst->graph->n_vertices;
	//cout<< node_pool.size() << " nodes in node_pool after run"<< endl;

	//cout << node_pool.size();
	node_it = node_pool.begin();
	bool hasExact = false;
	int count_low =0;
	int count_exact = 0;


	while (node_it != node_pool.end())	{

			//node_it->second->printState();

			if (node_it->second->exact){
				hasExact = true;
				//for the inexact cost calculations we need to double check
				node_it->second->cost = MAX(node_it->second->cost, calculateCost_bounds(node_it->second));
				//cout << "a" << node_it->second->cost << calculateCost_bounds(node_it->second);
				if (node_it->second->cost < upper_bound	){
					upper_bound = node_it->second->cost;
					//best_ub_node->state = State(node_it->second->state);

					//print ordering
					for (int pos = 0; pos < inst->graph->n_vertices; pos++){
						for  (int v = 0; v< inst->graph->n_vertices; v++)
							if (node_it->second->state[pos][v])	cout << v << ",";
					}cout << ": " << node_it->second->cost <<  endl;
					//
				}
				count_exact++;
			}

			if (node_it->second-> cost < min_cost_pool){
				min_cost_pool = MIN(min_cost_pool, node_it->second->cost);
				count_low = 1;

			}else if (node_it->second-> cost == min_cost_pool)
				count_low++;

			delete	 node_it->second;
			node_it++;
	}
	cout <<endl << "at lb: "<< count_low << " / " << node_pool.size() << "\t Exact: " << count_exact << " / " << node_pool.size()<< endl <<branch_nodes.size( ) << endl;
	lower_bound = MAX(lower_bound, MIN(min_cost_pool, lowest_unbranched));

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

	//best_lb = MAX(best_lb, lower_bound);
	//vertex_in_layer = orig_vertex_in_layer;
	for (vector<Node*>::iterator nit = nodes_layer.begin(); nit != nodes_layer.end(); ++nit)	delete *nit;
	nodes_layer.clear();
	node_pool.clear();

	return lower_bound;
}


//
// Merge nodes in a layer to meet maximum allowed width
//
void MinBandBDD::mergeLayer(int layer, vector<Node*> &nodes_layer) {
	//cout << "Merging layer "<< endl;
	sort(nodes_layer.begin(), nodes_layer.end(), CompareNodesCostDelta());
	/*for (vector<Node*>::const_iterator node_it =  nodes_layer.begin(); node_it != nodes_layer.end() ; ++node_it){
		(*node_it)->printState();
	} cout << endl;*/

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
		for( int i = 0; i < inst->graph->n_vertices; i++){
			(*central_state)[i] |= (*node)->state[i];
		}
		if ((*node)->exact) {
			addBranchNode((*node));
		}
		delete (*node);
	}
	//cout<< "finished merging"<< endl;

	// resize node layer vector
	nodes_layer.resize(maxWidth);

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
	} else if( var_ordering->order_type == FromFront ) {
		//do nothing
	} else if( var_ordering->order_type == FrontBack ) {
		//do nothing
	} else {
		cout << "Order undefined" << endl;
		exit(0);
	}

	//cout<< "Starting restriction ";
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
/*	if (active_state.size() < 1 ){
		vector<int> full_domain;
		for (unsigned char i=0; i<= (int)(std::ceil(inst->graph->n_vertices/2)); i++) full_domain.insert(i);
		active_state.push_back(full_domain);
	}*/

	//initialise all feasible root state //TODO copy active state
		State root_state;
		Domain fulldomain;
		for(int i =0; i<fulldomain.size(); i++){
			Domain domain = ~fulldomain;
			root_state.push_back(domain);
		}

	/*//print the initial state just to check
	for( std::vector<set<myint> >::const_iterator i = active_state.begin(); i != active_state.end(); ++i){
		for( std::set<myint>::const_iterator j = (*i).begin(); j != (*i).end(); ++j){
			std::cout << *j ;
		}
		cout << ",";
	}*/

	// create initial BDD node
	Node* initial_node = new Node(root_state, 0, true);
	node_pool.clear();
	node_pool[ &(root_state) ] = initial_node;
	root_node = initial_node;

	BDDNodePool::iterator node_it, existing_node_it;
	Node* node;

	// relaxation control variables
	int current_vertex;
	//int layer = root_state.size() -1 ;
	int layer =0; //todo
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

		}else if ( var_ordering->order_type == FrontBack ) {
			BDD* temp = new BDD();
			current_vertex = var_ordering->vertex_in_layer(temp, layer);
			active_vertices.erase(std::remove(active_vertices.begin(), active_vertices.end(), current_vertex), active_vertices.end());
			delete temp;
		}else if ( var_ordering->order_type == FromFront ) {
			BDD* temp = new BDD();
			current_vertex = var_ordering->vertex_in_layer(temp, layer);
			active_vertices.erase(std::remove(active_vertices.begin(), active_vertices.end(), current_vertex), active_vertices.end());
			delete temp;
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

		for (vector<Node*>::iterator it = nodes_layer.begin(); it != nodes_layer.end(); ++it) {

			branch_node = (*it);

			//cout << "Branching on node "  ;
			//branch_node->printState();

			//domains should already be filtered once we get  here
			//remove most recent domain, we are now branching on it
			Domain branch_domain = branch_node->state[layer];
			Domain domain;

			for (int v=0; v<inst->graph->n_vertices; v++){
				if (branch_domain[v]){
					domain.reset();
					domain.set(v);
					//domain.insert(*v);

					node = new Node(branch_node->state,
							branch_node->cost ,
							branch_node->exact);

					//remove
					node->state[layer] &= domain;

					//node->state.push_back(domain);

					//cout << "Filtering..";
					if (node->filterDomains5(layer) >= 0 &&  filterBounds2(node->state) >=0   ){
						int cost=0;
						if (vertex_in_layer[vertex_in_layer.size()-1]== 0 or vertex_in_layer[vertex_in_layer.size()-1]== inst->graph->n_vertices-1)
							node->cost = MAX(node->cost, inst->caprara_list[v]);
						else
							cost =  calculateCost_caprara_pos(node, v);
						node->cost = MAX(node->cost,cost);

						cost =  calculateCost_mu2(node);
						node->cost = MAX(node->cost,cost);

						if (node->cost < upper_bound)
							cost = calculateCost_ILP2(node);
						node->cost = MAX(node->cost,cost );


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
			}

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

				//
				//best_ub_node->state = State(node_it->second->state);
			}
		}
		node_it++;
	}

	best_ub = MIN(best_ub, ub);
	cout  << "restriction ub  = " 	<< best_ub << "," << lower_bound << endl;

	//edges_to_check.clear();

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


int MinBandBDD::calculateCost_caprara_gen(Node* _node){
	int outside_max = -1;
	vector<bool> false_vector(inst->graph->n_vertices, false);
	//todo smarter initialisation/

	//for every position
	for (int i=0; i<inst->graph->n_vertices/2+.1; i++){
		//cout << "x " << i << " x " << inst->graph->n_vertices -1 -i<< endl;
		int position_min = inst->graph->n_vertices;

		//for evrey vertex allowed in that position
		for (int v = 0; v < inst->graph->n_vertices; v++){
			if (_node->state[i][v]){
				vector<vector<int> > layered_graph = inst->graph->vertex_neighbourhood[v];
				int cumulative_vertices = 0;
				int internal_max = 0;

				for (vector<vector<int> >::iterator layer = layered_graph.begin()+1; layer!=layered_graph.end(); ++layer){
					int disjoint_domains =0;

					cumulative_vertices = cumulative_vertices + (*layer).size();
					int layer_number = layer - layered_graph.begin();

					//calculate current range
					int low_range, high_range;
					if(cumulative_vertices > 2*i){
						low_range = 0;
						high_range = cumulative_vertices ;
					}
					else{
						low_range = i - cumulative_vertices/2;
						high_range = i + cumulative_vertices/2;
					}
					//count vertices whose domains do not intersect range
					bool intersect = false;
					for (int vert = 0; vert < inst->graph->n_vertices; vert++){
						if (inst->graph->dist(vert,v) <= layer_number){
							intersect = false;
							for (int poss = low_range; poss< high_range; poss++)
								if (_node->state[poss][vert]) intersect = true;
							if (!intersect) disjoint_domains++;
						}
					}

					if (layer_number > 0){
						internal_max = std::max(internal_max,
								(int)std::ceil((1.*(cumulative_vertices+disjoint_domains > 2*i
										? cumulative_vertices - i + disjoint_domains
										: (cumulative_vertices+disjoint_domains)/2) - 1 ) /layer_number));
					}
				}
				//account for case layer 0
				if (internal_max > 0)
					position_min = (int)std::min(position_min, internal_max);
			}
		}
		//cout << position_min << endl;

		int position_min2 = inst->graph->n_vertices;


		for (int v = 0; v < inst->graph->n_vertices; v++){
			if (_node->state[inst->graph->n_vertices-i-1][v]){
				vector<vector<int> > layered_graph = inst->graph->vertex_neighbourhood[v];
				int cumulative_vertices = 0;
				int internal_max = 0;

				for (vector<vector<int> >::iterator layer = layered_graph.begin()+1; layer!=layered_graph.end(); ++layer){
					int disjoint_domains = 0;

					cumulative_vertices = cumulative_vertices + (*layer).size();
					int layer_number = layer - layered_graph.begin();


					//calculate current range
					int low_range, high_range;
					if(cumulative_vertices > 2*i){
						low_range = inst->graph->n_vertices - cumulative_vertices;
						high_range = inst->graph->n_vertices ;
					}
					else{
						low_range = i - cumulative_vertices/2;
						high_range = i + cumulative_vertices/2;
					}

					//count vertices whose domains do not intersect range
					bool intersect = false;
					for (int vert = 0; vert < inst->graph->n_vertices; vert++){
						if (inst->graph->dist(vert,v) <= layer_number){
							intersect = false;
							for (int poss = low_range; poss< high_range; poss++)
								if (_node->state[poss][vert]) intersect = true;
							if (!intersect) disjoint_domains++;
						}
					}

					if (layer_number > 0){
						internal_max = std::max(internal_max,
								(int)std::ceil((1.*(cumulative_vertices+disjoint_domains > 2*i
										? cumulative_vertices - i + disjoint_domains
										: (cumulative_vertices+disjoint_domains)/2) - 1 ) /layer_number));
					}
				}
				//account for case layer 0
				if (internal_max > 0)
					position_min2 = (int)std::min(position_min2, internal_max);
			}
		}
		//cout << position_min2 << endl;


		outside_max = (int)std::max(outside_max,
				std::max(position_min, position_min2));

	}
	return outside_max;

}
int MinBandBDD::calculateCost_caprara_pos(Node* _node, int _vertex){
	int pos = vertex_in_layer[vertex_in_layer.size()-1];
	int cumulative_vertices = 0;
	int internal_max = 0;

	for (vector<vector<int> >::iterator layer = inst->graph->vertex_neighbourhood[_vertex].begin()+1; layer!=inst->graph->vertex_neighbourhood[_vertex].end(); ++layer){
		int disjoint_domains = 0;

		cumulative_vertices = cumulative_vertices + (*layer).size();
		int layer_number = layer - inst->graph->vertex_neighbourhood[_vertex].begin();

		if (layer_number > 0){
			internal_max = std::max(internal_max,
					(int)std::ceil((1.*(cumulative_vertices > 2*pos
							? cumulative_vertices - pos
									: (cumulative_vertices)/2) - 1 ) /layer_number));
		}
	}

	return internal_max;

}

int MinBandBDD::calculateCost_caprara_fixed(Node* _node){
	int bound1 = inst->graph->n_vertices;;
	int boundn = inst->graph->n_vertices;
	for (int i = 0; i< inst->graph->n_vertices; i++){
		if (_node->state[0][i])
			bound1 = min(bound1, inst->caprara_list[i]);
		if (_node->state[inst->graph->n_vertices-1][i])
			boundn = min(boundn, inst->caprara_list[i]);
	}
	return max(bound1, boundn);
}
/*
int MinBandBDD::calculateCost_caprara(Node* _node){
	vector<bool> possible_start(inst->graph->n_vertices, true);

	for(vector<Domain>::const_iterator state_it = _node->state.begin(); state_it!= _node->state.end()-1; ++ state_it){
		if ((*state_it).find(1) == (*state_it).end()	){
			possible_start[vertex_in_layer[state_it - _node->state.begin()]] = false;
		}
		else if ((*state_it).size() == 1){ //if the node is fixed there
			return inst->caprara_list[vertex_in_layer[state_it - _node->state.begin()]] ;
		}
	}

	int caprara_min = inst->graph->n_vertices;

	//for every vertex that can still be in position 1
	for (int i = 0; i<inst->graph->n_vertices; i++){
		if (possible_start[i])
			caprara_min = std::min(caprara_min, inst->caprara_list[i]);
	}

	return caprara_min;
}*/

int MinBandBDD::calculateCost_bounds(Node* _node){
	/* Only use the largest nd smalles t element in every domain to calc
	 * costs
	 */
	int largest_smallest_cost = 0;

//	vector<int> lb(inst->graph->n_vertices,-1);
	//vector<int> ub(inst->graph->n_vertices,-1);
	for (int i = 0; i<inst->graph->n_vertices; i++){
		lb[i]  = -1; ub[i] = -1;
		for (int pos = 0; pos<inst->graph->n_vertices; pos++){
			if(_node->state[pos][i]){
				if (lb[i]<0)
					lb[i] = pos;
				ub[i] = pos;
			}
		}
		if (lb[i] == -1 or ub[i] == -1) return -1;
	}

	for (int vi = 0; vi < inst->graph->n_vertices-1; vi++){
		for (int vj = vi+1; vj< inst->graph->n_vertices; vj++){

			if (inst->graph->adj_m[vi][vj]){
				int smallest_cost_edge = inst->graph->n_vertices+1;

				//contained, dist 1
				if (lb[vi] <= lb[vj] and ub[vj] <= ub[vi]) smallest_cost_edge = 1;
				if (lb[vj] <= lb[vi] and ub[vi] <= ub[vj]) smallest_cost_edge = 1;
				//overlapping, dist 1
				if (lb[vi] <= lb[vj] and ub[vi] <= ub[vj] and ub[vi] >= lb[vj]) smallest_cost_edge = 1;
				if (lb[vj] <= lb[vi] and ub[vj] <= ub[vi] and ub[vj] >= lb[vi]) smallest_cost_edge = 1;
				//disjoint, calc dist
				if (ub[vi] <= lb[vj]) smallest_cost_edge = std::abs(ub[vi] - lb[vj]);
				if (ub[vj] <= lb[vi]) smallest_cost_edge = std::abs(ub[vj] - lb[vi]);

				if (smallest_cost_edge > largest_smallest_cost
						and smallest_cost_edge <= inst->graph->n_vertices){

					largest_smallest_cost = smallest_cost_edge;
					if (largest_smallest_cost > upper_bound)
							return largest_smallest_cost;
				}
			}

		}
	}

	return largest_smallest_cost;

}

/*
int MinBandBDD::calculateCost_bounds_fast(Node* _node){
	 Only use the largest nd smalles t element in every domain to calc
	 * costs


	//what happens to vertex_in_layer between iterations? reset?
	while (edges_to_check.size() > vertex_in_layer.size())
		edges_to_check.pop_back();


	//see if we are on a new level, map the ne wedges if so
	//this assumes that all nodes have the same ordering
	//cout << "Starting edges to check update " << vertex_in_layer.size() << " " << edges_to_check.size() << endl;
	vector<int>  edges_renamed;

	while(vertex_in_layer.size() > edges_to_check.size()){
		edges_renamed.clear();
		std::vector<int>::iterator from = vertex_in_layer.begin()+edges_to_check.size();

		for (vector<int>::iterator to = vertex_in_layer.begin(); to != from; ++to){
			if (inst->graph->adj_m[*to][*from] ){
				edges_renamed.push_back( to - vertex_in_layer.begin() );
				//cout << "addindg edge (" << edges_to_check.size()<< "," << to-vertex_in_layer.begin() <<") ";
			}
		}
		edges_to_check.push_back(edges_renamed);
	}

	// we need to have at least two domains to calculate a cost.
	if (vertex_in_layer.size() < 2 || _node->state.size()<2)
		return 0;

	int largest_smallest_cost = 0;

	if (_node->state.size()% 10 == 0 ){
		//for (int i = 0; i < MIN(_node->state.size(), vertex_in_layer.size())-1; i++){
		for (int i = 0; i < edges_to_check.size(); i++){
			//cout << "domain1= "<<i<< " - " << edges_to_check[i].size() << " :" << (edges_to_check[i].begin() -edges_to_check[i].end()) ;
			for (vector<int>::const_iterator j = edges_to_check[i].begin(); j != edges_to_check[i].end(); ++j){
				//cout << "enter"<< i<< *j << ";"<<_node->state[i].size()<< _node->state[*j].size() << endl;
				int smallest_cost_edge = inst->graph->n_vertices+1;

				//handle case where both domains are everything
				if (_node->state[i].size() == inst->graph->n_vertices or _node->state[*j].size() == inst->graph->n_vertices){
					smallest_cost_edge = 1;
				}
				else{
					int lb1 = *(_node->state[i].begin());
					int ub1 = *(--_node->state[i].end());

					int lb2 = *(_node->state[*j].begin());
					int ub[vi] = *(--_node->state[*j].end());

					//contained, dist 1
					if (lb1 <= lb2 and ub2 <= ub1) smallest_cost_edge = 1;
					if (lb2 <= lb1 and ub1 <= ub2) smallest_cost_edge = 1;
					//overlapping, dist 1
					if (lb1 <= lb2 and ub1 <= ub2 and ub1 >= lb2) smallest_cost_edge = 1;
					if (lb2 <= lb1 and ub2 <= ub1 and ub2 >= lb1) smallest_cost_edge = 1;
					//disjoint, calc dist
					if (ub1 <= lb2) smallest_cost_edge = std::abs(ub1 - lb2);
					if (ub2 <= lb1) smallest_cost_edge = std::abs(ub2 - lb1);

				}

				//cout << smallest_cost_edge << " ";
				if (smallest_cost_edge > largest_smallest_cost
						and smallest_cost_edge <= inst->graph->n_vertices){
					largest_smallest_cost = smallest_cost_edge;

					if (largest_smallest_cost > upper_bound)
										return largest_smallest_cost;
				}
				//cout << "exit" <<endl;
			}
			//cout << endl;
		}
	}
	else{

		int i = edges_to_check.size()-1;
		for (vector<int>::const_iterator j = edges_to_check[i].begin(); j != edges_to_check[i].end(); ++j){
			//cout << "enter"<< i<< *j << ";"<<_node->state[i].size()<< _node->state[*j].size() << endl;
			int smallest_cost_edge = inst->graph->n_vertices+1;

			//handle case where both domains are everything
			if (_node->state[i].size() == inst->graph->n_vertices or _node->state[*j].size() == inst->graph->n_vertices){
				smallest_cost_edge = 1;
			}
			else{
				int lb1 = *(_node->state[i].begin());
				int ub1 = *(--_node->state[i].end());

				int lb2 = *(_node->state[*j].begin());
				int ub2 = *(--_node->state[*j].end());

				//contained, dist 1
				if (lb1 <= lb2 and ub2 <= ub1) smallest_cost_edge = 1;
				if (lb2 <= lb1 and ub1 <= ub2) smallest_cost_edge = 1;
				//overlapping, dist 1
				if (lb1 <= lb2 and ub1 <= ub2 and ub1 >= lb2) smallest_cost_edge = 1;
				if (lb2 <= lb1 and ub2 <= ub1 and ub2 >= lb1) smallest_cost_edge = 1;
				//disjoint, calc dist
				if (ub1 <= lb2) smallest_cost_edge = std::abs(ub1 - lb2);
				if (ub2 <= lb1) smallest_cost_edge = std::abs(ub2 - lb1);

			}

			//cout << smallest_cost_edge << " ";
			if (smallest_cost_edge > largest_smallest_cost
					and smallest_cost_edge <= inst->graph->n_vertices){
				largest_smallest_cost = smallest_cost_edge;

				if (largest_smallest_cost > upper_bound)
					return largest_smallest_cost;
			}
			//cout << "exit" <<endl;
		}
	}

	//cout << "Done "<< largest_smallest_cost << endl;
	return largest_smallest_cost;
}
*/


int MinBandBDD::calculateCost_mu2(Node* _node){
	/* We assume we have a dummy domain that we dont care about,
	 * the domain for the next layer has already been added but
	 * there is nothing in vertex_in_layer yet. -2, -1 in for loops.
	 */

	//vector<int> lb(inst->graph->n_vertices,-1);
	//vector<int> ub(inst->graph->n_vertices,-1);
	for (int i = 0; i<inst->graph->n_vertices; i++){
		lb[i]  = -1; ub[i] = -1;

		for (int pos = 0; pos<inst->graph->n_vertices; pos++){
			if(_node->state[pos][i]){
				if (lb[i]<0)
					lb[i] = pos;
				ub[i] = pos;
			}
		}
		if (lb[i] == -1 or ub[i] == -1) return -1;
	}

	int largest_guaranteed_cost = 0;

	for (int vi = 0; vi < inst->graph->n_vertices-1; vi++){
		for (int vj = vi+1; vj< inst->graph->n_vertices; vj++){

			if (inst->graph->adj_m[vi][vj]){
				int smallest_cost_edge = inst->graph->n_vertices+1;

				//contained, dist 1
				if (lb[vi] <= lb[vj] and ub[vj] <= ub[vi]) smallest_cost_edge = 1;
				if (lb[vj] <= lb[vi] and ub[vi] <= ub[vj]) smallest_cost_edge = 1;
				//overlapping, dist 1
				if (lb[vi] <= lb[vj] and ub[vi] <= ub[vj] and ub[vi] >= lb[vj]) smallest_cost_edge = 1;
				if (lb[vj] <= lb[vi] and ub[vj] <= ub[vi] and ub[vj] >= lb[vi]) smallest_cost_edge = 1;
				//disjoint, calc dist
				if (ub[vi] <= lb[vj]) smallest_cost_edge = std::abs(ub[vi] - lb[vj]);
				if (ub[vj] <= lb[vi]) smallest_cost_edge = std::abs(ub[vj] - lb[vi]);

				largest_guaranteed_cost = MAX(largest_guaranteed_cost,
						(int)std::ceil(smallest_cost_edge/(1.*inst->graph->dist(vi,vj))));

				if (largest_guaranteed_cost > upper_bound)
					return largest_guaranteed_cost;
			}
		}
	}

	//cout << "Done "<< largest_smallest_cost << endl;
	return largest_guaranteed_cost;

}
/*

int MinBandBDD::calculateCost_mu2_fast(Node* _node){
	// we need to have at least two domains to calculate a cost.
	if (vertex_in_layer.size() < 2 || _node->state.size()<2)
		return 0;

	int largest_guaranteed_cost = 0;

	if (_node->state.size() % 10 == 0){
		for (int i = 0; i < _node->state.size() -2; i++){
			//cout << "domain1= "<<i<< " - " << edges_to_check[i].size() << " :" << (edges_to_check[i].begin() -edges_to_check[i].end()) ;
			for (int j = i+1; j != _node->state.size()-1; ++j){

				//cout << "enter"<< i<< *j << ";"<<_node->state[i].size()<< _node->state[*j].size() << endl;
				int smallest_cost_edge = inst->graph->n_vertices+1;

				//handle case where both domains are everything
				if (_node->state[i].size() == inst->graph->n_vertices or _node->state[j].size() == inst->graph->n_vertices){
					smallest_cost_edge = 1;
				}
				else{
					int lb1 = *(_node->state[i].begin());
					int ub1 = *(--_node->state[i].end());

					int lb2 = *(_node->state[j].begin());
					int ub2 = *(--_node->state[j].end());

					//contained, dist 1
					if (lb1 <= lb2 and ub2 <= ub1) smallest_cost_edge = 1;
					if (lb2 <= lb1 and ub1 <= ub2) smallest_cost_edge = 1;
					//overlapping, dist 1
					if (lb1 <= lb2 and ub1 <= ub2 and ub1 >= lb2) smallest_cost_edge = 1;
					if (lb2 <= lb1 and ub2 <= ub1 and ub2 >= lb1) smallest_cost_edge = 1;
					//disjoint, calc dist
					if (ub1 <= lb2) smallest_cost_edge = std::abs(ub1 - lb2);
					if (ub2 <= lb1) smallest_cost_edge = std::abs(ub2 - lb1);

				}

				largest_guaranteed_cost = MAX(largest_guaranteed_cost,
						(int)std::ceil(smallest_cost_edge/(1.*inst->graph->dist(vertex_in_layer[i], vertex_in_layer[j]))));

				if (largest_guaranteed_cost > upper_bound)
								return largest_guaranteed_cost;
				//cout << smallest_cost_edge << " ";

				//cout << "exit" <<endl;
			}
			//cout << endl;
		}

	}
	else{
		// vil is more general than state.size (do we hav the next domain already or not?)
		int i = vertex_in_layer.size() - 1;
		//cout << "domain1= "<<i<< " - " << edges_to_check[i].size() << " :" << (edges_to_check[i].begin() -edges_to_check[i].end()) ;
		for (int j = 0; j <i; ++j){

			//cout << "enter"<< i<< *j << ";"<<_node->state[i].size()<< _node->state[*j].size() << endl;
			int smallest_cost_edge = inst->graph->n_vertices+1;

			//handle case where both domains are everything
			if (_node->state[i].size() == inst->graph->n_vertices or _node->state[j].size() == inst->graph->n_vertices){
				smallest_cost_edge = 1;
			}
			else{
				int lb1 = *(_node->state[i].begin());
				int ub1 = *(--_node->state[i].end());

				int lb2 = *(_node->state[j].begin());
				int ub2 = *(--_node->state[j].end());

				//contained, dist 1
				if (lb1 <= lb2 and ub2 <= ub1) smallest_cost_edge = 1;
				if (lb2 <= lb1 and ub1 <= ub2) smallest_cost_edge = 1;
				//overlapping, dist 1
				if (lb1 <= lb2 and ub1 <= ub2 and ub1 >= lb2) smallest_cost_edge = 1;
				if (lb2 <= lb1 and ub2 <= ub1 and ub2 >= lb1) smallest_cost_edge = 1;
				//disjoint, calc dist
				if (ub1 <= lb2) smallest_cost_edge = std::abs(ub1 - lb2);
				if (ub2 <= lb1) smallest_cost_edge = std::abs(ub2 - lb1);

			}

			largest_guaranteed_cost = MAX(largest_guaranteed_cost,
					(int)std::ceil(smallest_cost_edge/(1.*inst->graph->dist(vertex_in_layer[i], vertex_in_layer[j]))));

			//cout << smallest_cost_edge << " ";

			//cout << "exit" <<endl;
		}
		//cout << endl;
	}	//cout << "Done "<< largest_smallest_cost << endl;
	return largest_guaranteed_cost;

}
*/
/*
int MinBandBDD::calculateCost_mu1(Node* _node){
	// we need to have at least two domains to calculate a cost.
	if (vertex_in_layer.size() < 2 || _node->state.size()<2)
		return 1;

	int largest_guaranteed_cost = 0;
	set<myint>::const_iterator i_val,j_val;

	//for (int i = 0; i < MIN(_node->state.size(), vertex_in_layer.size())-1; i++){
	for (int i = 0; i < _node->state.size() -2; i++){
		//cout << "domain1= "<<i<< " - " << edges_to_check[i].size() << " :" << (edges_to_check[i].begin() -edges_to_check[i].end()) ;
		for (int j = i+1; j != _node->state.size()-1; ++j){

			//cout << "enter"<< i<< *j << ";"<<_node->state[i].size()<< _node->state[*j].size() << endl;
			int smallest_cost_edge = inst->graph->n_vertices+1;

			for (i_val = _node->state[i].begin(); i_val != _node->state[i].end(); ++i_val){
				for (j_val = _node->state[j].begin(); j_val != _node->state[j].end() ; ++j_val){

					//todo only need bounds with all of other range

					if ((*i_val) != (*j_val) and (std::abs((*i_val) - (*j_val)) < smallest_cost_edge)){
						smallest_cost_edge = std::abs((*i_val) - (*j_val));
					}

				}
			}

			largest_guaranteed_cost = MAX(largest_guaranteed_cost,
								(int)std::ceil(smallest_cost_edge/(1.*inst->graph->dist(vertex_in_layer[i], vertex_in_layer[j]))));
			if (largest_guaranteed_cost > upper_bound)
				return largest_guaranteed_cost;
			//todo return as soon as bigger than upperbound

			//cout << "exit" <<endl;
		}
		//cout << endl;
	}

	//cout << "Done "<< largest_smallest_cost << endl;
	return largest_guaranteed_cost;

}*/

int MinBandBDD::calculateCost_ILP(Node* _node){


	int low_phi =  max(lower_bound, inst->caprara_bound) ;
	low_phi =  max(low_phi, inst->half_density_bound_b) ;

	int hi_phi = inst->graph->n_vertices ;
    int phi = (int) low_phi + (hi_phi-low_phi)/2;

	vector<int> ordering(inst->graph->n_vertices,-1);
	vector<int> fixed_pos;
	vector<int> fixed_v;

	//vector<int> lb(inst->graph->n_vertices,-1);
	//vector<int> ub(inst->graph->n_vertices,-1);
	for (int i = 0; i<inst->graph->n_vertices; i++){
		lb[i]  = -1; ub[i] = -1;

		for (int pos = 0; pos<inst->graph->n_vertices; pos++){
			if(_node->state[pos][i]){
				if (lb[i]<0)
					lb[i] = pos;
				ub[i] = pos;
			}
		}
		if (lb[i] == -1 or ub[i] == -1) return -1;
		if (lb[i] == ub[i]){
			fixed_pos.push_back(lb[i]);
			ordering[lb[i]] = i;
			fixed_v.push_back(i);
		}
	}

	vector<vector<int> > data_vertices ;

	while (low_phi - hi_phi < 0){
		data_vertices.clear();

		//compute all the l_v, f_v
		for (int i =0;i< inst->graph->n_vertices; i++){
			vector<int> trip(3,-1);
			trip[0] = i;
			trip[1] = inst->graph->n_vertices; //fv
			trip[2] = -1; //lv

			//calculate lv,fv
			for(int f = 0; f< fixed_pos.size(); f++){
				int h = inst->graph->dist(i,fixed_v[f]);
				//min {h*phi +i }
				trip[1] = std::min(trip[1], - h*phi + fixed_pos[f] );
				trip[2] = std::max(trip[2], h*phi + fixed_pos[f] );
			}

			//within bounds
			trip[1] = std::max(trip[1], 0);
			trip[2] = std::min(trip[2], inst->graph->n_vertices);

			//if we have a domain for this vertex, update bounds
			trip[1] = std::max(trip[1], lb[i]);
			trip[2] = std::min(trip[2], ub[i]);

			//cout << trip[0] <<"," <<trip[1]<< "," << trip[2]<<endl;
			data_vertices.push_back(trip);
		}

		//sort in order of f_v, want to use popback
		std::sort(data_vertices.begin(), data_vertices.end(), mytriplecompfv());

		//print data

		//vector<vector<int> > feasible_vertices;
		std::priority_queue<vector<int>, vector<vector<int> >,  mytriplecomplv> feasible_vertices;

		int pos=0;
		bool feasible_flag = true;

		while(pos < inst->graph->n_vertices and feasible_flag){


			while(data_vertices.size() >0 and data_vertices.back()[1] <= pos ){
				feasible_vertices.push(data_vertices.back());
				data_vertices.pop_back();
			}
			if (feasible_vertices.size() == 0){
				//cout << "size ";
				feasible_flag = false;
			}
			else{

				vector<int> insert = feasible_vertices.top();
				feasible_vertices.pop();

				//test last position violated? insert if not
				if (insert[2] < pos){
					//cout << "last";
					//cout << pos << " d "<<insert[0] <<"," <<insert[1]<< "," << insert[2]<<endl;
					feasible_flag = false;
				}
				else
					ordering[pos] = insert[0];
			}

			++pos;
		}

		//cout << phi << "-"<< feasible_flag << "...";
		//update binary search bounds
		if (feasible_flag){
			hi_phi = phi-1;
			phi = (int) low_phi + (hi_phi-low_phi)/2;
		}
		else{
			low_phi = phi+1;
			phi = (int) low_phi + (hi_phi-low_phi)/2;

			if (phi > upper_bound)
				return phi;
		}

	}

	//cout << phi << "="; _node->printState();

	return phi;

}

int MinBandBDD::calculateCost_ILP2(Node* _node){


	int low_phi =  max(lower_bound, inst->caprara_bound) ;
	low_phi =  max(low_phi, inst->half_density_bound_b) ;

	int hi_phi = upper_bound ;
    int phi = (int) low_phi + (hi_phi-low_phi)/2;

	vector<int> ordering(inst->graph->n_vertices,-1);
	vector<int> fixed_pos;
	vector<int> fixed_v;

	//vector<int> lb(inst->graph->n_vertices,-1);
	//vector<int> ub(inst->graph->n_vertices,-1);
	for (int i = 0; i<inst->graph->n_vertices; i++){
		lb[i]  = inst->graph->n_vertices; ub[i] = -1;

		//Only works for frontback ordering
		int pos = 0;
		for (int pos=0; pos< inst->graph->n_vertices; pos=pos+2){
			if(_node->state[pos][i]){
				lb[i] = std::min(lb[i], pos/2);
				ub[i]  = std::max(ub[i], pos/2);
			}
		}
		for (int pos=1; pos< inst->graph->n_vertices; pos=pos+2){
			if(_node->state[pos][i]){
				lb[i] = std::min(lb[i],inst->graph->n_vertices - 1 - (pos-1)/2);
				ub[i]  = std::max(ub[i],inst->graph->n_vertices - 1 - (pos-1)/2);
			}
		}/*
		for (int pos = 0; pos<inst->graph->n_vertices; pos++){
			if(_node->state[pos][i]){
				if (lb[i]<0)
					lb[i] = pos;
				ub[i] = pos;
			}
		}*/
		if (lb[i] == inst->graph->n_vertices or ub[i] == -1) return -1;
		if (lb[i] == ub[i]){
			fixed_pos.push_back(lb[i]);
			ordering[lb[i]] = i;
			fixed_v.push_back(i);
		}
	}
	//cout << "a" << endl;
	//v, s_v, l_v (f_v for t_vect)
	vector<vector<int> > s_vect(inst->graph->n_vertices, vector<int>(4,-1));
	vector<vector<int> > t_vect(inst->graph->n_vertices, vector<int>(4,-1));

	vector<int> Left;
	vector<int> Right;

	int ii = 0;
	//cout << "Left ";
	while (ii < inst->graph->n_vertices and ordering[ii] !=-1 ){
		//cout << ordering[ii] <<" "	;
		Left.push_back(ordering[ii++]);
	}
	ii = inst->graph->n_vertices-1;
	//cout << "Right ";
	while (ii >= 0 and ordering[ii] != -1 ){
		//cout << ordering[ii] << " ";
		Right.push_back(ordering[ii--]);
	}
	//cout << " -- "<<  Left.size()<< "," << Right.size() << endl;

	for (int i =0; i< inst->graph->n_vertices; i++){
		s_vect[i][0] = i;
		t_vect[i][0] = i;
		s_vect[i][1] = inst->graph->n_vertices;
		t_vect[i][1] = inst->graph->n_vertices;
		for(int j = 0; j< Left.size() ; j++){
			//cout << "d("<< i<< "," <<Left[j] <<")=" <<  inst->graph->dist(i,Left[j])<<" ";
			s_vect[i][1] = MIN(s_vect[i][1], inst->graph->dist(i,Left[j]));
		}
		for(int j = 0; j< Right.size() ; j++)
			t_vect[i][1] = MIN(t_vect[i][1], inst->graph->dist(i,Right[j]));
		//ifthe vertex is in L or R, then flag =0 means fixed.
		if (t_vect[i][1]== 0 or s_vect[i][1] == 0)
		{
			t_vect[i][3] = 0;
			s_vect[i][3] = 0;
		}

		//cout  << i << "[" << s_vect[i][1] << "," << t_vect[i][1] << "] " ;
	}
	//cout  << endl;


	vector<vector<int> > data_vertices ;

	while (low_phi - hi_phi < 0){
		// sort, we will soon update s/t_vect by direct indexing,
		// this is only valid if it is in order.
		std::sort(s_vect.begin(), s_vect.end(), mycompinc());
		std::sort(t_vect.begin(), t_vect.end(), mycompinc());

		//cout << "phi = " << phi  << endl;
		data_vertices.clear();
		//cout << phi << endl;
		//compute all the l_v, f_v
		for (int i =0;i< inst->graph->n_vertices; i++){
			vector<int> trip(3,-1);
			// (v, fv, lv) based on phi and distances and fixed vertices.
			trip[0] = i;
			trip[1] = -1; //fv
			trip[2] = inst->graph->n_vertices; //lv

			//calculate lv,fv
			for(int f = 0; f< fixed_pos.size(); f++){
				int h = inst->graph->dist(i,fixed_v[f]);
				//min {h*phi +i }
				trip[1] = std::max(trip[1], - h*phi + fixed_pos[f] );
				trip[2] = std::min(trip[2], h*phi + fixed_pos[f] );
			}

			//within bounds
			trip[1] = std::max(trip[1], 0);
			trip[2] = std::min(trip[2], inst->graph->n_vertices);

			//if we have a domain for this vertex, update bounds
			trip[1] = std::max(trip[1], lb[i]);
			trip[2] = std::min(trip[2], ub[i]);
			if (lb[i] == ub[i]){
				trip[1] = lb[i];
				trip[2] = ub[i];
			}

			s_vect[i][2] = trip[2];
			t_vect[i][2] = trip[1];

			//cout << trip[0] <<"(" <<trip[1]<< "," << trip[2]<<")"<<" ";
			data_vertices.push_back(trip);
		}
		//cout <<"L" << endl;
		if (Left.size()> 0){
			//inner LP for computing l_v
			//for (int i = 0; i< inst->graph->n_vertices; i++) cout << s_vect[i][0] << "-"<< s_vect[i][1]<< "-" <<s_vect[i][2]<<"-" <<s_vect[i][3] << ","; cout << endl;
			std::sort(s_vect.begin(), s_vect.end(), mydoublecompincinc()); //increasing order of s_v, decresing in l_v
			//for (int i = 0; i< inst->graph->n_vertices; i++) cout << s_vect[i][0] << "-"<< s_vect[i][1]<< "-" <<s_vect[i][2]<<"-" <<s_vect[i][3] << ","; cout << endl;
			for (int i = 0; i< inst->graph->n_vertices; i++){
				//cout << s_vect[i][1] << "(" << s_vect[i][2] << ")";
				if (s_vect[i][1] > 0 and s_vect[i][3] != 0 ){

					if (s_vect[i][1] > s_vect[i-1][1])
						std::sort(s_vect.begin(), s_vect.begin()+i, mydoublecompincinc());

					//cout << s_vect[i][0] <<"~";
					int j = i-1;
					int pos = inst->graph->n_vertices;
					while (j>=0 and s_vect[j][1] >= s_vect[i][1] - 1){
						//onl for things in N_1^L
						if ( s_vect[j][1]  < s_vect[i][1] and  inst->graph->dist(s_vect[j][0], s_vect[i][0]) == 1){
						//	cout << s_vect[j][0]<<"{"<<s_vect[j][2]<<"}";
							if (pos == inst->graph->n_vertices) // first one goes in l_v, next goes min
								pos = s_vect[j][2];
							else
								pos = MIN(pos-1, s_vect[j][2]);
							//cout <<"=" << pos;
						}

						j--;
					}

					s_vect[i][2] = MIN(s_vect[i][2], pos+phi);
					//cout << ":"<< s_vect[i][2]<< " ";
					data_vertices[s_vect[i][0]][2] = MIN(data_vertices[s_vect[i][0]][2], pos+phi);
				}
			}
		}
		//cout << endl << "R ";
		//inner lp computing f_v
		if (Right.size() > 0){
			std::sort(t_vect.begin(), t_vect.end(), mydoublecompincdec());
			//for (int i = 0; i< inst->graph->n_vertices; i++) cout << t_vect[i][0] << "-"<< t_vect[i][1]<< "-" <<t_vect[i][2]<<"-"<< t_vect[i][3]<< ","; cout << endl;
			for (int i = 0; i< inst->graph->n_vertices; i++){
				//cout << t_vect[i][1];
				if (t_vect[i][1] > 0 and t_vect[i][3] != 0 ){

					if (t_vect[i][1] > t_vect[i-1][1])
							std::sort(t_vect.begin(), t_vect.begin()+i, mydoublecompincdec());
					//cout << t_vect[i][0] <<"~";

					int j = i-1;
					int pos =-1;
					while (j>=0 and t_vect[j][1] >= t_vect[i][1] - 1){
						//onl for things in N_1^L (so at dist 1 )
						if ( t_vect[j][1]  < t_vect[i][1] and  inst->graph->dist(t_vect[j][0], t_vect[i][0]) == 1  ){
						//	cout << t_vect[j][0]<<"{"<<t_vect[j][2]<<"}";

							if (pos == -1) // first one goes in f_v, next goes min
								pos = t_vect[j][2];
							else
								pos = MAX(pos+1, t_vect[j][2]);
							//cout <<"=" << pos;

						}
						j--;
					}
					t_vect[i][2] =MAX(t_vect[i][2], pos-phi);
					//cout << ":"<< t_vect[i][2]<< " ";

					data_vertices[t_vect[i][0]][1] = MAX(data_vertices[t_vect[i][0]][1],pos-phi);
				}
			}
		}
		//cout<< endl;
		//for (int i = 0; i< inst->graph->n_vertices; i++) cout << data_vertices[i][0]<<"(" <<data_vertices[i][1] << "," << data_vertices[i][2] << ") ";cout <<endl;

		//sort in order of f_v, want to use popback, so decreasing order
		std::sort(data_vertices.begin(), data_vertices.end(), mytriplecompfv());
//		cout << "after sort: ";
		///for (int i = 0; i< inst->graph->n_vertices; i++) cout << data_vertices[i][0]<<"(" <<data_vertices[i][1] << "," << data_vertices[i][2] << ") ";cout <<endl;

		//print data

		//vector<vector<int> > feasible_vertices;
		std::priority_queue<vector<int>, vector<vector<int> >,  mytriplecomplv> feasible_vertices;

		int pos=0;
		bool feasible_flag = true;

		while(pos < inst->graph->n_vertices and feasible_flag){

			// if f_v <= pos, potential vertex
			while(data_vertices.size() >0 and data_vertices.back()[1] <= pos ){
				feasible_vertices.push(data_vertices.back());
				data_vertices.pop_back();
			}
			//cout << pos << feasible_vertices.size() << endl;
			if (feasible_vertices.size() == 0){
				//cout << "size ";
				feasible_flag = false;
			}
			else{

				vector<int> insert = feasible_vertices.top();
				feasible_vertices.pop();
				//cout << insert[2] << "-" ;
				//if	 (feasible_vertices.size() > 0) cout << feasible_vertices.top()[2] << endl;

				//test last position violated? insert if not
				if (insert[2] < pos){
					//cout << "last";
//					cout << pos << " d "<<insert[0] <<"," <<insert[1]<< "," << insert[2]<<endl;
					feasible_flag = false;
				}
				else
					ordering[pos] = insert[0];
			}

			++pos;
		}

//		cout << low_phi << "," << phi << "," << hi_phi  << "-"<< feasible_flag << "..."<<endl;
		//update binary search bounds
		if (feasible_flag){
			hi_phi = phi;
			phi = (int) low_phi + (hi_phi-low_phi)/2;

			//calculate the upper bound from the orderoing we just got
			int val = 0;
			for (int i=0; i < inst->graph->n_vertices-1; i++)
				for (int j=i+1; j< inst->graph->n_vertices; j++)
					if (inst->graph->dist(ordering[i], ordering[j]) == 1)
						val = MAX(val, j-i );

			if (val < upper_bound) {
				cout << "setting upperbound = " << val << endl;
				/*for (int pos = 0; pos < inst->graph->n_vertices; pos++)
					cout << ordering[pos] << ","; cout << ": " << val << endl;*/
				upper_bound = val;
			}

		}
		else{
			low_phi = phi+1;
			phi = (int) low_phi + (hi_phi-low_phi)/2;

			if (phi > upper_bound){
				//cout << "returning, val=" << phi << endl;
				return phi;
			}
		}
		//cout << low_phi << "," << phi << "," << hi_phi <<endl;
	}

	//cout << phi << "="; _node->printState();


	//cout << "returning, val=" << phi << endl;

	return phi;

}


/*
int MinBandBDD::calculateCost(Node* _node){

	//TODO different way of calculating cost: only singletons, maybe even only newest vertex

	//cout<< "calculating costs" << endl;
	//cout << "calculating costs on ";
	// _node->printState() ;

	//what happens to vertex_in_layer between iterations? reset?
	while (edges_to_check.size() > vertex_in_layer.size())
		edges_to_check.pop_back();


	//see if we are on a new level, map the ne wedges if so
	//this assumes that all nodes have the same ordering
	//cout << "Starting edges to check update " << vertex_in_layer.size() << " " << edges_to_check.size() << endl;
	std::vector<int>::iterator from,to;
	while(vertex_in_layer.size() > edges_to_check.size()){
			vector<int>  edges_renamed;
			from = vertex_in_layer.begin()+edges_to_check.size();

			for (to = vertex_in_layer.begin(); to != from; ++to){
				if (inst->graph->adj_m[*to][*from] ){
					edges_renamed.push_back( to - vertex_in_layer.begin() );
					//cout << "addindg edge (" << edges_to_check.size()<< "," << to-vertex_in_layer.begin() <<") ";
				}
			}
			edges_to_check.push_back(edges_renamed);
		}
	//cout <<"Finishing edges to check update" << endl;

	// we need to have at least two domains to calculate a cost.
	if (vertex_in_layer.size() < 2 || _node->state.size()<2)
		return 0;

	int largest_smallest_cost = 0;
	vector<int>::const_iterator j;
	set<myint>::const_iterator i_val,j_val;

	//for (int i = 0; i < MIN(_node->state.size(), vertex_in_layer.size())-1; i++){
	for (int i = 0; i < edges_to_check.size(); i++){
		//cout << "domain1= "<<i<< " - " << edges_to_check[i].size() << " :" << (edges_to_check[i].begin() -edges_to_check[i].end()) ;
		for (j = edges_to_check[i].begin(); j != edges_to_check[i].end(); ++j){
			//cout << "enter"<< i<< *j << ";"<<_node->state[i].size()<< _node->state[*j].size() << endl;
			int smallest_cost_edge = inst->graph->n_vertices+1;
			bool smallEnough = false;

			//handle case where both domains are everything
			if (_node->state[i].size() == inst->graph->n_vertices or _node->state[*j].size() == inst->graph->n_vertices){
				smallest_cost_edge = 1;
			}
			else{

				for ( i_val = _node->state[i].begin(); i_val != _node->state[i].end() and !smallEnough; ++i_val){

					//test if ew are already had longer on any of the edges, if so this will only get smaller, bail
					//test if we are close enough to i_val that we ca decrease smallest_edge, if not no point in checking
					//for (set<myint>::const_iterator j_val = _node->state[*j].begin(); j_val != _node->state[*j].end() and !smallEnough and ((*j_val) <= (*i_val) + smallest_cost_edge); ++j_val){
					for ( j_val = _node->state[*j].begin(); j_val != _node->state[*j].end() and !smallEnough; ++j_val){
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

				if (largest_smallest_cost > upper_bound)
					return largest_smallest_cost;
			}
			//cout << "exit" <<endl;
		}
		//cout << endl;
	}

	//cout << "Done "<< largest_smallest_cost << endl;
	return largest_smallest_cost;

}

double MinBandBDD::calculateCost_bounds_delta(Node* _node){


	//what happens to vertex_in_layer between iterations? reset?
	while (edges_to_check.size() > vertex_in_layer.size())
		edges_to_check.pop_back();


	//see if we are on a new level, map the ne wedges if so
	//this assumes that all nodes have the same ordering

	//this may get called at any stage, add all edges we missed
	std::vector<int>::iterator from,to;
	vector<int>  edges_renamed;
	while(vertex_in_layer.size() > edges_to_check.size()){
		edges_renamed.clear();
		from = vertex_in_layer.begin() + edges_to_check.size();

		for (to = vertex_in_layer.begin(); to != from; ++to){
			if (inst->graph->adj_m[*to][*from] ){
				edges_renamed.push_back( to - vertex_in_layer.begin() );
			}
		}
		edges_to_check.push_back(edges_renamed);
	}

	// we need to have at least two domains to calculate a cost.
	if (vertex_in_layer.size() < 2 || _node->state.size()<2)
		return 0;

	int largest_smallest_cost = 0;
	vector<int> cost_counts(upper_bound+1,0);

	//for (int i = 0; i < MIN(_node->state.size(), vertex_in_layer.size())-1; i++){
	for (int i = 0; i < edges_to_check.size(); i++){
		//cout << "domain1= "<<i<< " - " << edges_to_check[i].size() << " :" << (edges_to_check[i].begin() -edges_to_check[i].end()) ;
		for (vector<int>::const_iterator j = edges_to_check[i].begin(); j != edges_to_check[i].end(); ++j){
			//cout << "enter"<< i<< *j << ";"<<_node->state[i].size()<< _node->state[*j].size() << endl;
			int smallest_cost_edge = inst->graph->n_vertices+1;

			//handle case where both domains are everything
			if (_node->state[i].size() == inst->graph->n_vertices or _node->state[*j].size() == inst->graph->n_vertices){
				smallest_cost_edge = 1;
			}
			else{
				int lb1 = *(_node->state[i].begin());
				int ub1 = *(--_node->state[i].end());

				int lb2 = *(_node->state[*j].begin());
				int ub2 = *(--_node->state[*j].end());

				//contained, dist 1
				if (lb1 <= lb2 and ub2 <= ub1) smallest_cost_edge = 1;
				if (lb2 <= lb1 and ub1 <= ub2) smallest_cost_edge = 1;
				//overlapping, dist 1
				if (lb1 <= lb2 and ub1 <= ub2 and ub1 >= lb2) smallest_cost_edge = 1;
				if (lb2 <= lb1 and ub2 <= ub1 and ub2 >= lb1) smallest_cost_edge = 1;
				//disjoint, calc dist
				if (ub1 <= lb2) smallest_cost_edge = std::abs(ub1 - lb2);
				if (ub2 <= lb1) smallest_cost_edge = std::abs(ub2 - lb1);

			}

			if (smallest_cost_edge < upper_bound)
				cost_counts[smallest_cost_edge]++;
			//cout << smallest_cost_edge << " ";
			if (smallest_cost_edge > largest_smallest_cost
					and smallest_cost_edge <= inst->graph->n_vertices){
				largest_smallest_cost = smallest_cost_edge;

				if (largest_smallest_cost >= upper_bound)
						return largest_smallest_cost;
			}
			//cout << "exit" <<endl;
		}
		//cout << endl;
	}
	//calculate the additional terms of delta
	vector<int> products(largest_smallest_cost+1,0);
	products[largest_smallest_cost] = inst->graph->n_vertices+1;
	for(int i=largest_smallest_cost-1; i>=0; i--)
		products[i] = products[i+1]*(inst->graph->n_vertices + largest_smallest_cost - i +1 );

	double sum =0;
	for (int i=1; i<=largest_smallest_cost; i++)
		sum = sum + ((1.*cost_counts[i])/products[i]);

	_node->cost_delta = sum;
	//cout << largest_smallest_cost + sum<< endl;
	return largest_smallest_cost;

}
*/
int MinBandBDD::filterBounds2(State& state){
	//only filters the layer we are branching on
	//cout << "fb";

	//vector<int> lb(inst->graph->n_vertices,-1);
	//vector<int> ub(inst->graph->n_vertices,-1);

	for (int i = 0; i<inst->graph->n_vertices; i++){
		lb[i]  = -1; ub[i] = -1;

		for (int pos = 0; pos<inst->graph->n_vertices; pos++){
			if(state[pos][i]){
				if (lb[i]<0)
					lb[i] = pos;
				ub[i] = pos;
			}
		}
		if (lb[i] == -1 or ub[i] == -1) return -1;
	}
	for (int i = 0; i<inst->graph->n_vertices; i++){
			for (int pos = 0; pos<inst->graph->n_vertices; pos++){
				if(state[pos][i]){
					if (lb[i]<0)
						lb[i] = pos;
					ub[i] = pos;
				}
			}
			if (lb[i] == -1 or ub[i] == -1) return -1;
		}

	for (int i = 0; i<inst->graph->n_vertices; i++){
		for (int j = i; j<i+1; j++){
			int jump = (upper_bound) * inst->graph->dist(i,j);

			//prune state[.][i]
			if (lb[i] < lb[j] - jump){
				for (int pos=lb[i]; pos<lb[j] - jump; pos++){
						state[pos].reset(i);
				}
				lb[i] = lb[j] - jump;
			}
			if (ub[i] > ub[j] + jump){
				for (int pos=ub[j] + jump; pos<inst->graph->n_vertices; pos++){
					state[pos].reset(i);
				}
				ub[i] = ub[j]+jump;
			}

			//prune state[.][i]
			if (lb[j] < lb[i] - jump){
				for (int pos=lb[j]; pos<lb[i] - jump; pos++){
					state[pos].reset(i);
				}
				lb[j] = lb[i] - jump;
			}
			if (ub[j] > ub[i] + jump){
				for (int pos=ub[i] + jump; pos<inst->graph->n_vertices; pos++){
					state[pos].reset(i);
				}
				ub[j] = ub[i]+jump;
			}
		}
	}

	//test if any row/column is empty
	for (int i=0; i <inst->graph->n_vertices; i++){
		if (state[i].none())
			return -1;
		/*int count = 0;
		for (int j=0; j<inst->graph->n_vertices; j++)
			if (state[j][i] )
				count++;
		if (count == 0)
			return -1;*/
	}


	return 0;
}

/*int MinBandBDD::filterBounds(State& state){
	//cout << "fb";

	vector<int> lb(inst->graph->n_vertices,-1);
	vector<int> ub(inst->graph->n_vertices,-1);

	for (int i = 0; i<inst->graph->n_vertices; i++){
		for (int pos = 0; pos<inst->graph->n_vertices; pos++){
			if(state[pos][i]){
				if (lb[i]<0)
					lb[i] = pos;
				ub[i] = pos;
			}
		}
		if (lb[i] == -1 or ub[i] == -1) return -1;
	}

	for (int i = 0; i<inst->graph->n_vertices; i++){
		for (int j = i; j<i+1; j++){
			int jump = (upper_bound-1) * inst->graph->dist(i,j);

			//prune state[.][i]
			if (lb[i] < lb[j] - jump){
				for (int pos=lb[i]; pos<lb[j] - jump; pos++){
					state[pos].reset(i);
				}
				lb[i] = lb[j] - jump;
			}
			if (ub[i] > ub[j] + jump){
				for (int pos=ub[j] + jump; pos<inst->graph->n_vertices; pos++){
					state[pos].reset(i);
				}
				ub[i] = ub[j]+jump;
			}

			//prune state[.][i]
			if (lb[j] < lb[i] - jump){
				for (int pos=lb[j]; pos<lb[i] - jump; pos++){
					state[pos].reset(i);
				}
				lb[j] = lb[i] - jump;
			}
			if (ub[j] > ub[i] + jump){
				for (int pos=ub[i] + jump; pos<inst->graph->n_vertices; pos++){
					state[pos].reset(i);
				}
				ub[j] = ub[i]+jump;
			}

		}
	}

	//test if any row/column is empty
	for (int i=0; i <inst->graph->n_vertices; i++){
		if (state[i].none())
			return -1;
		int count = 0;
		for (int j=0; j<inst->graph->n_vertices; j++)
			if (state[j][i] )
				count++;
		if (count == 0)
			return -1;
	}


	return 0;
}*/
