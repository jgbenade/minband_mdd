/*
 * ---------------------------------------------------------
 * BDD Structure definitions for independet set
 * ---------------------------------------------------------
 */

#ifndef BDD_HPP_
#define BDD_HPP_

#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <set>
#include <algorithm>
//#include <boost/dynamic_bitset.hpp>
#include "util.hpp"

using namespace std;


//
// State
//
//typedef boost::dynamic_bitset<> State;
typedef set<int> Domain;
typedef vector<Domain> State;

//
// Lexicographic state comparator
//
struct StateLessThan {
	bool operator()(const State* stateA, const State* stateB) const {
		return (*stateA) < (*stateB);
	} //TODO does this make sense here?
};


//
// BDD Node
//
struct Node {
	State			state;
	int				cost;
	bool			exact;

	//todo
	//vector<Node>	children;
	int 			layer;
	vector<std::vector<set<int> >::iterator> marked; //used for filtering


	Node(State &_state, double _cost)
	: state(_state), cost(_cost), exact(false), layer(-1)
	{ }

	/*Node(State &_state, double _cost, bool _exact)
	: state(_state), cost(_cost), exact(_exact),  layer(-1)
	{ state.reserve(24);}*/
	Node(State &_state, double _cost, bool _exact)
		: state(_state), cost(_cost), exact(_exact),  layer(-1)
		{ }

	int filterDomains();
	int filterDomains2();
	int filterDomains3();

	void printState();

};

inline int Node::filterDomains2(){
	//cout << "Before "; printState();
	int feasible = 0;
	//vector<std::vector<set<int> >::iterator> marked;

	for (std::vector<set<int> >::iterator domain=state.begin(); domain !=state.end(); ++domain	){
		//marked.clear();
		//marked.push_back(domain);

		//only for small sets, too expensive
		if ((*domain).size()==1){

			std::vector<set<int> >::iterator domain2 = state.begin();

			//count number of other domains that is a subset of this one, see if you can do an obv Hall set
			while (domain2!= state.end()){
				// if this domain is included in the previous one
				if (domain != domain2){
					//std::set<int>::iterator itlow,itup;
					//itlow = (*domain2).lower_bound(*( (*domain).begin() ) - ub* inst->graph->dist()

					if (std::includes( (*domain2).begin(), (*domain2).end(),
							(*domain).begin(), (*domain).end() )	){

						(*domain2).erase( *( (*domain).begin() ) );

						if ((*domain2).size() == 0 and domain2 < state.end()-1)
							return -1;
					}

				}
				++domain2;
			}
		}

		//printState();

	}
	//cout << "After "; printState();


	/*for (std::vector<set<int> >::iterator domain=state.begin(); domain !=state.end(); ++domain	){
		if ((*domain).size() == 0)
			return -1;
	}*/
	return feasible;
}

inline int Node::filterDomains3(){
	//last domain is full before filtering, use this to get n_vertices
	vector<bool> flag(state[state.size()-1].size(), false);
	set<int> fixed;

	for (int i=0; i<state.size(); i++){
		if (state[i].size()==1)
			fixed.insert(*(state[i].begin()) );
	}
	int prev_fixed_size;
	if (fixed.size() > 0){
		do{
			prev_fixed_size =fixed.size();
			for (int i=0; i<state.size(); i++){

				if (state[i].size()==1){
					if (flag[ *(state[i].begin() ) ]){
						return -1;
					}else{
						flag[ *(state[i].begin() ) ] = true;
					}
				}else{

					/*std::set_difference(state[i].begin(), state[i].end(),
							fixed.begin(), fixed.end(),
							inserter(state[i],state[i].begin() ) );*/
					for (set<int>::iterator val   = fixed.begin(); val != fixed.end(); ++ val)
						state[i].erase(*val);

					if (state[i].size() == 0)
						return -1;
					else
						if (state[i].size() ==1){
							fixed.insert(*(state[i].begin()) );
							flag[*(state[i]).begin()] = true;
						}
				}
			}

		}
		while(prev_fixed_size != fixed.size() /*and counter< 2*/);
	}
	return 1;

}

inline int Node::filterDomains(){
	//cout << "Before "; printState();
	//cout << "F " ; printState();
	int feasible = 0;
	//vector<std::vector<set<int> >::iterator> marked;

	for (std::vector<set<int> >::iterator domain=state.begin(); domain !=state.end(); ++domain	){
		marked.clear();
		marked.push_back(domain);

		//only for small sets, too expensive
		if ((*domain).size()<3){
			std::vector<set<int> >::iterator domain2 = state.begin();

			//count number of other domains that is a subset of this one, see if you can do an obv Hall set
			while (marked.size() < (*domain).size() and domain2!= state.end()){
				// if this domain is included in the previous one
				if (domain!= domain2){
					if (std::includes( (*domain).begin(), (*domain).end(),
							(*domain2).begin(), (*domain2).end()) ){
						std::vector<set<int> >::iterator it(domain2);
						marked.push_back(it);
					}
				}
				++domain2;
			}

			//correct #, filter obvious Hall set.
			if (marked.size() == (*domain).size()){
				domain2 = state.begin();
				while ( domain2!= state.end()){
					// if this domain is included in the previous one

					if (std::find(marked.begin(), marked.end(), domain2) == marked.end()){
						//Domain 	pruned_domain;

						//std::sort((*domain2).begin(), (*domain2).end());
						//std::sort((*domain).begin(), (*domain).end());
						//std::vector<int>::iterator it = std::set_difference( (*domain2).begin(), (*domain2).end(),
						//						(*domain).begin(), (*domain).end(),
						//						std::inserter(pruned_domain, pruned_domain.end() ) );

						for (std::set<int>::iterator value = (*domain).begin(); value!= (*domain).end(); ++value){
							if ( (*domain2).count(*value) )
								(*domain2).erase(*value);
						}

						//pruned_domain.resize(it - pruned_domain.begin());
						//(*domain2) = pruned_domain;
					}
					++domain2;
				}
			}
		}

		//printState();

	}
	// if the very last state is empty it means that we are ons the last layer
	for (std::vector<set<int> >::iterator domain=state.begin(); domain !=state.end()-1; ++domain	){
		if ((*domain).size() == 0 )
			return -1;
	}
	return feasible;
}

inline void Node::printState(){
	for( std::vector<set<int> >::const_iterator i = state.begin(); i != state.end(); ++i){
		cout << "/";
		for( std::set<int>::const_iterator j = (*i).begin(); j != (*i).end(); ++j){
			std::cout << *j ;
		}

	}
	cout << ":" << cost << ":"<< ( (exact) ? "Ex" : "Rel") << endl;
}


//
// BDD data structure
//
struct BDD {

};

//
// BDD Node pool type
//
//typedef map<State*, Node*, StateLessThan> BDDNodePool;
typedef unordered_map<State*, Node*> BDDNodePool;


//
// Marker of the end of a state (for iteration purposes)
//
//const int state_end = static_cast<int>(vector<>::npos);


/*
 * ---------------------------------------------------------
 * BDD Node comparators
 * ---------------------------------------------------------
 */

//
// Node comparator by longest path
//
struct CompareNodesCost {
	bool operator()(const Node* nodeA, const Node* nodeB) const {
		return nodeA->cost < nodeB->cost;
	}
};


#endif /* BDD_HPP_ */
