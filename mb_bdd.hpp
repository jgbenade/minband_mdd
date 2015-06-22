/*
 * ---------------------------------------------------------
 * BDD Structure definitions for independet set
 * ---------------------------------------------------------
 */

#ifndef BDD_HPP_
#define BDD_HPP_

#include <iostream>
#include <map>
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


	Node(State &_state, double _cost)
	: state(_state), cost(_cost), exact(false), layer(-1)
	{ }

	Node(State &_state, double _cost, bool _exact)
	: state(_state), cost(_cost), exact(_exact),  layer(-1)
	{ }

	int filterDomains();

	void printState();
};


inline int Node::filterDomains(){
	int feasible = 0;
	vector<std::vector<set<int> >::iterator> marked;

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

	for (std::vector<set<int> >::iterator domain=state.begin(); domain !=state.end(); ++domain	){
		if ((*domain).size() == 0)
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
typedef map<State*, Node*, StateLessThan> BDDNodePool;

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
