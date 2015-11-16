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
#include <bitset>
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
typedef std::bitset<15> Domain;
typedef vector<Domain> State;

//
// Lexicographic state comparator
//
struct StateLessThan {
	bool operator()(const State* stateA, const State* stateB) const {
		for (int pos=0; pos < (*stateA).size(); pos++){
			for (int v =0; v< (*stateA).size(); v++){
				if ( (*stateA)[pos][v]!= (*stateB)[pos][v]){
					if ((*stateA)[pos][v]) return true;
						else return false;
				}
			}
		}
		return false;
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

	//int filterDomains();
	int filterDomains2();
	//int filterDomains3();
	//int filterDomains4();
	int filterDomains5(int pos);

	void printState();

};

inline int Node::filterDomains5(int pos){
	//cout << "fd";
	if (exact){
		int n = (*state.end()).size();
		int setidx = -1;

		for (int i=0; i<state[pos].size(); i++)
			if (state[pos][i])
				setidx = i;

		int counter = 0;
		//all we have to do is delete all the other things (only later domains, exact), cant be infeasible
		for (std::vector<Domain >::iterator domain=state.begin(); domain != state.end(); ++domain	){
			if (counter != pos)
				(*domain).reset(setidx);
			counter++;
		}
		return 1;
	}
	else
		return filterDomains2();
}
inline int Node::filterDomains2(){
	//cout << "fd2";
	//filters all states, including the last one (but doesnt return false if the last one is empty)
	int feasible = 0;
	//vector<std::vector<set<myint> >::iterator> marked;
	std::vector<Domain >::iterator domain2 ;

	for (std::vector<Domain >::iterator domain=state.begin(); domain !=state.end(); ++domain	){
		//marked.clear();
		//marked.push_back(domain);

		//only for small sets, too expensive
		if ((*domain).count()==1){

			domain2 = state.begin();
			while (domain2 != state.end()){

				if (domain != domain2){
					//flip the single set bit in domain, do and to prune
					(*domain2) &= ~(*domain);
					if (!(*domain2).any()){
						return -1;
					}


				}
				++domain2;
			}
		}
	}

	for (int vert = 0; vert < state.size(); vert++){
		//for every vertex, see how many feasible positions it has
		int countPos = 0;
		int setpos = 0;
		for (int pos = 0; pos < state.size(); pos++){
			if (state[pos][vert]){
				countPos++;
				setpos = pos;
			}
		}

		if (countPos==1){
			for (int vert2 = 0; vert2 < state.size(); vert2++){
				if (vert2!= vert){
					state[setpos].reset(vert2);
				}
			}
		}
		else if (countPos == 0){
			return -1;
		}
	}
	return feasible;
}

inline void Node::printState(){

	for( int i =0; i < state.size(); i++){
		cout << "/";
		for( int j=0; j< state[i].size(); j++){
			if  (state[i][j])
				std::cout <<j;
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
struct CompareNodesDomain {
	bool operator()(const Node* nodeA, const Node* nodeB) const {
		cout << "x"<<endl;
		int n = nodeA->state.size();
		for (int i=0; i< nodeA->state.size()/2; i++){
			cout << i << endl;
			if (nodeA->state[i].count() != nodeB->state[i].count())
				return nodeA->state[i].count() < nodeB->state[i].count();
			cout << "y"<< endl;

			if (nodeA->state[n-1-i].count() != nodeB->state[n-1-i].count() )
				return nodeA->state[n-1-i].count() < nodeB->state[n-1-i].count();
		}

		return true;
	}
};

#endif /* BDD_HPP_ */
