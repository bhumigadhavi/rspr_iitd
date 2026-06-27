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
	LEAF_REDUCTION = true;
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

	int m = inputData.numTrees;
	int n = inputData.numLeaves;
	if ((int)inputData.trees.size() < m)
		m = (int)inputData.trees.size();
	if (m < 2) {
		cerr << "# Error: need at least 2 trees." << endl;
		return 1;
	}

	// Build and label all trees
	vector<Node *> trees;
	for (int i = 0; i < m; i++)
		trees.push_back(build_tree(inputData.trees[i]));
	for (Node *t : trees)
		t->labels_to_numbers(&label_map, &reverse_label_map);
	if (n == 0)
		n = (int)label_map.size();

	// Among all pairs find the one with the highest 3-approx distance.
	// That pair gives the tightest lower-bound / best heuristic upper-bound.
	int best_approx = -1;
	int best_i = 0, best_j = 1;
	if (m == 2) {
		Forest fa = Forest(trees[0]);
		Forest fb = Forest(trees[1]);
		best_approx = rSPR_worse_3_approx(&fa, &fb);
	} else {
		for (int i = 0; i < m; i++) {
			for (int j = i + 1; j < m; j++) {
				Forest fa = Forest(trees[i]);
				Forest fb = Forest(trees[j]);
				int approx = rSPR_worse_3_approx(&fa, &fb);
				if (approx > best_approx) {
					best_approx = approx;
					best_i = i;
					best_j = j;
				}
			}
		}
	}

	Node *T1 = trees[best_i];
	Node *T2 = trees[best_j];

	// heuristic algorithm: exact B&B on the hardest pair
	if (CLUSTER_TEST) {
		T1->preorder_number();
		T1->edge_preorder_interval();
		T2->preorder_number();
		T2->edge_preorder_interval();

		Forest *out_F1 = NULL;
		Forest *out_F2 = NULL;
		int exact_k = rSPR_branch_and_bound_simple_clustering(T1, T2, false,
			&label_map, &reverse_label_map, best_approx / 3, min(best_approx, n - 1), &out_F1, &out_F2);
		if (out_F1 != NULL) {
			out_F1->numbers_to_labels(&reverse_label_map);
			out_F1->print_components();
			delete out_F1;
		}
		if (out_F2 != NULL)
			delete out_F2;
		cout << "# exact_k=" << exact_k << endl;
		for (Node *t : trees) t->delete_tree();
		return 0;
	}

	for (Node *t : trees) t->delete_tree();
}