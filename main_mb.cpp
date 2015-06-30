//============================================================================
// Name        : BDD_MinBand.cpp
//============================================================================

#include <ctime>
#include <deque>
#include <iostream>
#include <vector>
#include "minband_instance.hpp"
#include "minband_bdd.hpp"

#define MAX_NODES_MEMORY        1700000000
#define MIN_NODES_DFS		10000
#define TIME_LIMIT		1800

using namespace std;

//
// Buffer of branching nodes
// (Global to preserve memory)
//
vector<BranchNode*> branchNodeBuffer;


//
// Store branching nodes in file. Returns maximum upper bound of stored node
// TODO: do more efficiently: currently it is removing all nodes from queue and adding
// them back
//
inline int storeNodesInFile(const int num_file, const int size, BranchNodeQueue& queue) {
/*
	assert(size < queue.size());
	int maxUB = 0;

	char filename[256];
	sprintf(filename, "branches/input_%d.dat", num_file);

	// create buffer
	branchNodeBuffer.clear();
	for (size_t i = 0; i < queue.size(); ++i) {
		branchNodeBuffer.push_back( queue.top() );
		queue.pop();
	}

	// store node in file
	ofstream states(filename);
	states << size << endl;
	for (int i = 0; i < size; ++i) {

		// get node
		BranchNode* branchNode = branchNodeBuffer.back();
		branchNodeBuffer.pop_back();
		maxUB = MAX(maxUB, branchNode->relax_ub);

		// add to file
		states << branchNode->relax_ub;
		states << " " << branchNode->cost;
		states << " " << branchNode->state.size();
		for (vector<int>::iterator it = branchNode->state.begin(); it != branchNode->state.end(); ++it) {
			states << " " << *it << " ";
		}
		states << endl;

		// delete from memory
		delete branchNode;
	}
	states.close();

	// add remaining nodes back to queue
	for (vector<BranchNode*>::iterator it = branchNodeBuffer.begin(); it != branchNodeBuffer.end(); ++it) {
		queue.push(*it);
	}
	branchNodeBuffer.clear();

	return maxUB;*/
}


//
// Read branching nodes from file
//
inline void readNodesInFile(const int num_file, BranchNodeQueue& queue) {
	return;
	/*char filename[256];
	sprintf(filename, "branches/input_%d.dat", num_file);
	ifstream states(filename);

	int size;
	states >> size;
	for (int i = 0; i < size; ++i) {
		BranchNode* branchNode = new BranchNode;
		states >> branchNode->relax_lb;
		states >> branchNode->cost;

		int state_size;
		states >> state_size;
		branchNode->state.resize(state_size);
		for (int j = 0; j < state_size; ++j) {
			states >> branchNode->state[j];
		}

		queue.push(branchNode);
	}
	states.close();
	std::remove(filename);*/
}


// -------------------------------------------------------
// Global variables
// -------------------------------------------------------
int 					global_lb;
int 					global_ub;
long int 				nodes_explored;
clock_t 				init;
clock_t 				maxTime;
vector<BranchNode*>		        dfsBranchNodes;
unsigned long int                       size_pool;


//
// Perform DFS search on nodes in the queue
//
void performDFS(MinBandBDD& minband_bdd, BranchNodeQueue& queue) {
	//todo
	cout << "\n\nswitching to DFS..." << endl << endl;
	clock_t current;

	while (size_pool > MIN_NODES_DFS) {

		// update bound
		global_ub = MIN(global_ub, queue.top()->relax_lb);

		// take best node
		dfsBranchNodes.push_back( queue.top() );
		queue.pop();
		size_pool -= dfsBranchNodes.back()->state.size();

		// explore nodes to reduce size of the queue
		while (!dfsBranchNodes.empty()) {

			// take branch node from the top
			BranchNode* branch_node = dfsBranchNodes.back();
			dfsBranchNodes.pop_back();

			// update statistics
			nodes_explored++;

			if (nodes_explored % 100 == 0 ) {
			  //cout << "**Explored = " << nodes_explored;
				cout << "\t\tTo Explore = " << (queue.size() + dfsBranchNodes.size());
				cout << "\t\tLB = " << global_lb;
				cout << "\t\tUB = " << global_ub;
				//cout << "\t\tgap = " << ((double)(global_ub - global_lb) / global_ub)*100.0;
				cout << "\t\tpool = " << size_pool << endl;
				cout << endl;
			}

			// explores node if not pruned due to global lower bound
			if (branch_node->relax_lb > global_lb) {

				// explore node
				int ub = minband_bdd.generateRelaxation(branch_node);

				// update global lower bound if BDD is exact
				if (minband_bdd.isBDDExact()) {
					global_lb = std::max(global_lb, ub);
				} else {
					// Primal heuristic **************************************
					int lb = minband_bdd.generateRestriction(branch_node);
					minband_bdd.updateLocalUB(lb, true);
					global_lb = std::max(global_lb, lb);
					// *******************************************************
				}

				// add open nodes to pool
				minband_bdd.addBranchNodesDFS(dfsBranchNodes);
			}

			// erase branch node
			delete branch_node;

			// check time
			current = clock();
			if (current - init >= maxTime) {
				for (vector<BranchNode*>::iterator it = dfsBranchNodes.begin(); it != dfsBranchNodes.end(); ++it) {
					queue.push(*it);
				}
				return;
			}
		}

		// if bounds match, optimum has been found
		if (global_lb >= global_ub) {
			global_ub = global_lb;
			return;
		}
	}


	cout << "\n\nswitching back to priority queue..." << endl << endl;
}


void printSet(set<int>& _set){
	for (set<int>::const_iterator i = _set.begin(); i!=_set.end(); ++i){
		cout << *i;
	}
	cout << ',';
}

void testFilter(){
	set<int> d1,d2;
	d1.insert(1);
	d2.insert(2);
	d2.insert(3);

	set<int> d3(d2);

	set<int> d4(d2);
	d4.insert(1);
	d4.insert(4);

	set<int> d5(d4); d5.insert(5);
	State state;
	state.push_back(d4);
	state.push_back(d2);
	state.push_back(d3);
	state.push_back(d1);
	state.push_back(d5);

	Node node(state, 0);

	cout << "State before  filtering" << endl;
	node.printState();
	cout << "State after filtering: " << node.filterDomains2() << endl;
	node.printState();

	node.state.push_back(d5);
	cout << "State before  filtering" << endl;
	node.printState();
	cout << "State after filtering: "<<		node.filterDomains2() << endl;
	node.printState();

	State state2;
	d1.insert(2);
	state2.push_back(d1);
	state2.push_back(d2);
	d4.erase(1); d4.erase(2);
	Domain d6;
	d6.insert(3); d6.insert(4);
	Domain d7 ;
	d7.insert(1);
	d7.insert(4);
	state2.push_back(d6);
	state2.push_back(d7);
	state2.push_back(d5);

	Node node2(state2, 0);

	cout << "State before  filtering" << endl;
	node2.printState();
	cout << "State after filtering: " << node2.filterDomains2() << endl;
	node2.printState();


}
void testStates(){
	set<int> d1,d2;
	d1.insert(1);
	d2.insert(2);
	d2.insert(3);

	set<int> d3(d2);

	set<int> d4(d2);
	d4.insert(1);
	d4.insert(4);

	set<int> d5(d4); d5.insert(5);
	State state;
	state.push_back(d4);
	state.push_back(d2);
	state.push_back(d3);
	state.push_back(d1);
	state.push_back(d5);
	State state2(state);


	Node node(state, 0);
	Node node2(state, 0);

	cout << "State compare "<< (state == state2)<< endl ;
	state[1].insert(9);
	cout << "State compare "<< (state == state2)<< endl ;
	state2[1].insert(8);
	cout << "State compare "<< (state == state2)<< endl ;

	state.clear();
	state2.clear();
	cout << "State compare "<< (state == state2)<< endl ;
	set<int> a;
	set<int> b;
	state.push_back(a);
	state2.push_back(b);
	cout << "State compare "<< (state == state2)<< endl ;
	state[0].insert(33);
	cout << "State compare "<< (state == state2)<< endl ;
	state2[0].insert(32);
	cout << "State compare "<< (state == state2)<< endl ;
	state2[0].erase(32);
	state2[0].insert(33);
	cout << "State compare "<< (state == state2)<< endl ;





}


//
// Main function
//
int main(int argc, char* argv[]) {

	// TODO: add parameters
	if (argc < 4) {
		cout << "\nUsage: minband [graph file] [root BDD width] [BDD width]\n\n";
		cout << "\twhere: \n";
		cout << "\t\tgraph file: graph in DIMACS format\n";
		cout << "\t\troot BDD width: BDD width at the root node\n";
		cout << "\t\tBDD width: BDD width at branching nodes\n";
		cout << endl;
		cout << endl;
		cout << "** Set width to '-1' for either case to use size of state as max width" << endl;
		cout << endl;
		exit(1);
	}

	// read input
	int root_max_width = atoi(argv[2]);
	int max_width = atoi(argv[3]);

	// initialize time
	maxTime = TIME_LIMIT * CLOCKS_PER_SEC;
	init = clock();
	clock_t current;

	// ------------------------------------------------------
	// Initialize and process initial root node
	// ------------------------------------------------------

	//testStates();
	//testFilter();
	cout<< "Hi";


	MinBandBDD minband_bdd(root_max_width, max_width, argv[1]);
	current = clock();
	double totalTime = ((double)(current - init))/CLOCKS_PER_SEC;
	cout << endl << totalTime << endl;



	// Initialize global bounds
	if (minband_bdd.isBDDExact()) {
	  cout << "\troot relaxation is exact" << endl << endl;
	  global_ub = minband_bdd.getUB();
	  global_lb = global_ub;

	} else {
	  global_lb = minband_bdd.getBestLB();
	  global_ub = minband_bdd.getUB();
	}

	// approximate size of the pool
	size_pool = 0;

	// number of nodes explored
	nodes_explored = 1;

	// Create branching pool
	BranchNodeQueue branch_node_queue;
	minband_bdd.addBranchNodesQueue(branch_node_queue, size_pool);
	//cout<< " Brachnodespool size "<< branch_node_queue.size();

	// ------------------------------------------------------
	// Search
	// TODO: apply primal heuristic
	// ------------------------------------------------------

	// repeat until all search space is explored
	/*while (!branch_node_queue.empty()) {

		cout << "brnach_node_queue size after relaxation "<< branch_node_queue.size() << endl;
		for (BranchNodeQueue::iterator it = branch_node_queue. ; it != branch_node_queue.end(); ++it) {
			(*it)->printState();
		}
		cout<<endl;

		// switch temporally to DFS if necessary
		if (size_pool > MAX_NODES_MEMORY) {
			performDFS(minband_bdd, branch_node_queue);
		}

		// check optimality conditions
		if (branch_node_queue.empty()) {
			break;
		}

		// take branch node from the top
		BranchNode* branch_node = branch_node_queue.top();
		branch_node_queue.pop();

		// update statistics
		nodes_explored++;
		size_pool -= branch_node->state.size();

		if (nodes_explored % 1 == 0 ) {
		  //cout << "Explored = " << nodes_explored;
		  cout << endl;
		  cout << "to_explore = " << branch_node_queue.size();
		  cout << "  ::  lb = " << global_lb;
		  cout << "  ::  ub = " << global_ub;
		  cout << "  ::  gap = " << ((double)(global_ub - global_lb) / global_ub)*100.0 << "%";
		  //cout << "  ::  pool=" << size_pool;
		  cout << endl;
		}

		//cout<< branch_node->relax_lb << endl;

		// explores node if not pruned due to global upper bound
		if (branch_node->relax_lb < global_ub) {

			// explore node
			int lb = minband_bdd.generateRelaxation(branch_node);

			//cout << endl << "new lb " << lb << " " << minband_bdd.isBDDExact()	;

			// update global upper bound if BDD is exact
			if (minband_bdd.isBDDExact()) {
				global_ub = std::min(global_ub, lb);

			}
			if (minband_bdd.wasUBUpdated()) {
				//note that getUB returns upper_bound, not get_best_ub
				global_ub = std::min(global_ub, minband_bdd.getUB());

			}else {
				// Primal heuristic **************************************
				int ub = minband_bdd.generateRestriction(branch_node);
				//cout << "Found sol value " << ub;
				minband_bdd.updateLocalUB(ub, true);
				//cout << "updateing global_ub: " << global_ub << ub<< endl;
				global_ub = std::min(global_ub, ub);
				// *******************************************************
			}

			// add open nodes to pool
			minband_bdd.addBranchNodesQueue(branch_node_queue, size_pool);
		}

		// erase branch node
		delete branch_node;

		// update global upper bound
		if (!branch_node_queue.empty()) {
			//cout << "updateing global_lb: " << global_lb << branch_node_queue.top()->relax_lb<< endl;
			global_lb = MAX(global_lb, branch_node_queue.top()->relax_lb);
		} else {
			cout << "Branch nodes queue empty, setting global_lb = global_ub"<< endl;
			global_lb = global_ub;
		}

		// if bounds match, optimum has been found
		if (global_lb >= global_ub) {
		  global_ub = global_lb;
		  break;
		}

		// check time
		current = clock();
		if (current - init >= maxTime) {
			break;
		}
	}*/

	cout << endl;
	cout << "Lower bound = " << global_lb << endl;
	cout << "Upper bound = " << global_ub << endl;

	current = clock();
	totalTime = ((double)(current - init))/CLOCKS_PER_SEC;

	char statfilename[256];
	sprintf(statfilename, "stats.txt");

	ofstream statfile(statfilename, ios::app);
	statfile << argv[1];
	statfile << "\t" << nodes_explored;
	statfile << "\t" << root_max_width;
	statfile << "\t" << minband_bdd.inst->calculate_HalfDensity_Bound();
	statfile << "\t" << minband_bdd.inst->calculate_Caprara_Bound();
	statfile << "\t" << global_lb;
	statfile << "\t" << global_ub;
	statfile << "\t" << ((double)(global_ub - global_lb) / (double)global_ub)*100.0;
	statfile << "\t" << totalTime;
	statfile << endl;
	statfile.close();
	return 0;
}






