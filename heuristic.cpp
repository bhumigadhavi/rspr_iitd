#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <climits>
#include <vector>
#include <map>
#include <time.h>
#include <list>
#include "rspr.h"
#include <csignal>

#include "timeout.h"
#include "Forest.h"
#include "ClusterForest.h"
#include "LCA.h"
#include "ClusterInstance.h"
#include "UndoMachine.h"
#include "lgt.h"
#include "read.h"

using namespace std;

// This flag must be volatile sig_atomic_t to be safe for asynchronous modification
volatile sig_atomic_t timeout_signaled = 0;

// options to pick default
bool DEFAULT_ALGORITHM=true;
bool DEFAULT_OPTIMIZATIONS=true;


bool MULTI_APPROX = false;
bool FPT = false;
bool RF = false;
bool QUIET = true;
bool UNROOTED = false;
bool ALL_UNROOTED = false;
bool SIMPLE_UNROOTED = false;
bool SIMPLE_UNROOTED_RSPR = false;
bool LCA_TEST = false;
bool CLUSTER_TEST = true;
bool TOTAL = false;
bool PAIRWISE = false;
bool PAIRWISE_SYMMETRIC = true;
int PAIRWISE_START = 0;
int PAIRWISE_END = INT_MAX;
int PAIRWISE_COL_START = 0;
int PAIRWISE_COL_END = INT_MAX;
bool PAIRWISE_MAX = false;
int PAIRWISE_MAX_SPR = INT_MAX;
bool APPROX = false;
bool LOWER_BOUND = false;
bool REDUCE_ONLY = false;
bool PRINT_ROOTED_TREES = false;
bool SHOW_MOVES = false;
bool SEQUENCE = false;
bool DEBUG_REVERSE = false;
bool RANDOM_SPR = false;
int RANDOM_SPR_COUNT = 0;
int MULTI_TEST = 0;

// The function that the OS will call when time is up
void signal_handler(int signal) {
    timeout_signaled = 1;
}

int main() {
    cout << "# Hello, World!" << endl;

	// Register the signals
    signal(SIGTERM, signal_handler);
    signal(SIGINT, signal_handler);

	// # Core B&B pruning rules: O(1) checks per node, proven safe
	CUT_ONE_B = true;
	REVERSE_CUT_ONE_B = true;
	REVERSE_CUT_ONE_B_3 = true;
	CUT_TWO_B = true;

	// # Approx algorithm pruning: tightens lower bound, fewer k iterations needed
	APPROX_CUT_ONE_B = true;
	APPROX_REVERSE_CUT_ONE_B = true;

	// # Skip clustering overhead for instances where approx_spr < 30
	CLUSTER_TUNE = 30;

    // Label maps to allow string labels
	map<string, int> label_map= map<string, int>();
	map<int, string> reverse_label_map = map<int, string>();

	// set random seed
	srand(unsigned(time(0)));

	Input2 inputData = readInput();
	
	string T1_line = inputData.trees[0];
	string T2_line = inputData.trees[1];
	int n = inputData.numLeaves;
	Node *T1;
	Node *T2;
	
	T1 = build_tree(T1_line);
	T2 = build_tree(T2_line);

	
	// checking reduction rules
	if (REDUCE_ONLY) {

		T1->preorder_number();
		T1->edge_preorder_interval();
		T2->preorder_number();
		T2->edge_preorder_interval();

		
		Forest F1 = Forest(T1);
		Forest F2 = Forest(T2);
		
		F1.print_components();
		F2.print_components();

		sync_twins(&F1, &F2);
		
		if (MULTIFURCATING) {
			reduction_leaf_mult(&F1, &F2);
		}
		else {
			reduction_leaf(&F1, &F2);
		}
		
		F1.print_components();
		F2.print_components();
	}

	//print trees
	if (!QUIET) {
		cout << "T1: ";
		T1->print_subtree();
		cout << "T2: ";
		T2->print_subtree();
		cout << endl;
	}

	//check this out
	if (LCA_TEST) {
		LCA lca_query = LCA(T1);
		cout << endl;
		lca_query.debug();
		cout << endl;
		vector<Node *> leaves = T1->find_leaves();
		for(vector<Node *>::iterator i = leaves.begin(); i != leaves.end(); i++) {
			for(vector<Node *>::iterator j = i; j != leaves.end(); j++) {
				if (j==i)
					continue;
				Node *lca = lca_query.get_lca(*i, *j);
				(*i)->print_subtree_hlpr();
				cout << "\t";
				(*j)->print_subtree_hlpr();
				cout << "\t";
				lca->print_subtree();
			}
		}
		T1->delete_tree();
		T2->delete_tree();
		return 0;
	}
	
	T1->labels_to_numbers(&label_map, &reverse_label_map);
	T2->labels_to_numbers(&label_map, &reverse_label_map);

	// for debugging
	if (SHOW_MOVES) {
		show_moves(T1, T2, &label_map, &reverse_label_map);
		T1->delete_tree();
		T2->delete_tree();
		return 0;
	}

	//heuristic algorithm
	if (CLUSTER_TEST) {
		T1->preorder_number();
		T1->edge_preorder_interval();
		T2->preorder_number();
		T2->edge_preorder_interval();

		int exact_k = rSPR_branch_and_bound_simple_clustering(T1,T2,false, &label_map, &reverse_label_map, -1, n - 1, NULL, NULL);
		cout << "# exact_k=" << exact_k << endl;
		T1->delete_tree();
		T2->delete_tree();
		return 0;
	}

	// cleanup
	T1->delete_tree();
	T2->delete_tree();
}