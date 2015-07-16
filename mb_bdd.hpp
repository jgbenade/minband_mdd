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
//unsigned char has rang 0..255
//typedef unsigned char myint;
typedef  unsigned char myint;
typedef set<myint> Domain;
typedef vector<Domain> State;

//
// Lexicographic state comparator
//
struct StateLessThan {
	bool operator()(const State* stateA, const State* stateB) const {
		return (*stateA) < (*stateB);
	}
};


//
// BDD Node
//
struct Node {
	State			state;
	int				cost;
	double			cost_delta;
	bool			exact;

	//vector<Node>	children;
	int 			layer;

	vector<std::vector<set<myint> >::iterator> marked; //used for filtering
	//vector<set<myint> >::iterator domain, domain2;


	Node(State &_state, double _cost)
	: state(_state), cost(_cost), exact(false), layer(-1), cost_delta(0)
	{ }

	/*Node(State &_state, double _cost, bool _exact)
	: state(_state), cost(_cost), exact(_exact),  layer(-1)
	{ state.reserve(24);}*/
	Node(State &_state, double _cost, bool _exact)
	: state(_state), cost(_cost), exact(_exact),  layer(-1), cost_delta(0)
		{ }

	int filterDomains();
	int filterDomains2();
	int filterDomains3();
	int filterDomains4();
	int filterDomains5();

	void printState();

};

inline int Node::filterDomains5(){
	if (exact){
		int n = (*state.end()).size();
		std::vector<set<myint> >::iterator domain2 = --state.end();

		//all we have to do is delete all the other things, cant be infeasible

		for (std::vector<set<myint> >::iterator domain=state.begin(); domain != --state.end(); ++domain	){

			(*domain2).erase( *( (*domain).begin() ) );
		}
		return 1;
	}
	else
		return filterDomains2();
}
inline int Node::filterDomains2(){
	//filters all states, including the last one (but doesnt return false if the last one is empty)
	int feasible = 0;
	//vector<std::vector<set<myint> >::iterator> marked;
	std::vector<set<myint> >::iterator domain2 ;

	for (std::vector<set<myint> >::iterator domain=state.begin(); domain !=state.end(); ++domain	){
		//marked.clear();
		//marked.push_back(domain);

		//only for small sets, too expensive
		if ((*domain).size()==1){

			domain2 = state.begin();

			//count number of other domains that is a subset of this one, see if you can do an obv Hall set
			while (domain2 != state.end()){
				// if this domain is included in the previous one
				if (domain != domain2){
					//std::set<myint>::iterator itlow,itup;
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
	}
	return feasible;
}

inline int Node::filterDomains4(){
	set<myint> fixed;
	std::pair<std::set<myint>::iterator,bool> ret;

	for (int i=0; i<state.size(); i++){
		if (state[i].size()==1){
			ret = fixed.insert(*(state[i].begin()) );
			//test if we inserted a new value (if not duplicat singleton domains)
			if (!ret.second){
				return -1;
			}
		}
		else{
			for (set<myint>::iterator val   = fixed.begin(); val != fixed.end(); ++ val)
				state[i].erase(*val);

			if (state[i].size() == 0)
				return -1;
			else
				if (state[i].size() ==1){
					ret = fixed.insert(*(state[i].begin()) );
					if (!ret.second)
						return -1;
				}
		}
	}

	//and in reverse
	fixed.clear();
	for (int i=state.size()-1; i>=0; i--){
		if (state[i].size()==1){
			ret = fixed.insert(*(state[i].begin()) );
			//test if we inserted a new value (if not duplicat singleton domains)
			if (!ret.second){
				return -1;
			}
		}
		else{
			for (set<myint>::iterator val   = fixed.begin(); val != fixed.end(); ++ val)
				state[i].erase(*val);

			if (state[i].size() == 0)
				return -1;
			else
				if (state[i].size() ==1){
					ret = fixed.insert(*(state[i].begin()) );
					if (!ret.second)
						return -1;
				}
		}
	}

	return 0;
}
inline int Node::filterDomains3(){
	//last domain is full before filtering, use this to get n_vertices
	vector<bool> flag(state[state.size()-1].size(), false);
	set<myint> fixed;

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
					for (set<myint>::iterator val   = fixed.begin(); val != fixed.end(); ++ val)
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
	//vector<std::vector<set<myint> >::iterator> marked;

	for (std::vector<set<myint> >::iterator domain=state.begin(); domain !=state.end(); ++domain	){
		marked.clear();
		marked.push_back(domain);

		//only for small sets, too expensive
		if ((*domain).size()<3){
			std::vector<set<myint> >::iterator domain2 = state.begin();

			//count number of other domains that is a subset of this one, see if you can do an obv Hall set
			while (marked.size() < (*domain).size() and domain2!= state.end()){
				// if this domain is included in the previous one
				if (domain!= domain2){
					if (std::includes( (*domain).begin(), (*domain).end(),
							(*domain2).begin(), (*domain2).end()) ){
						marked.push_back(domain2);
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

						for (std::set<myint>::iterator value = (*domain).begin(); value!= (*domain).end(); ++value){
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
	for (std::vector<set<myint> >::iterator domain=state.begin(); domain !=state.end()-1; ++domain	){
		if ((*domain).size() == 0 )
			return -1;
	}
	return feasible;
}

inline void Node::printState(){
	for( std::vector<set<myint> >::const_iterator i = state.begin(); i != state.end(); ++i){
		cout << "/";
		for( std::set<myint>::const_iterator j = (*i).begin(); j != (*i).end(); ++j){
			std::cout << (int)(*j) ;
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
typedef unordered_map<State*, Node*>  BDDNodePool;


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
struct CompareNodesCostDelta {
	bool operator()(const Node* nodeA, const Node* nodeB) const {
		return nodeA->cost+nodeA->cost_delta < nodeB->cost+nodeB->cost_delta;
	}
};

#endif /* BDD_HPP_ */
