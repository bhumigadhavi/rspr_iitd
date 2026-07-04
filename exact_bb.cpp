// exact_bb.cpp — Simultaneous multi-forest exact solver for drSPR
//
// Differs from exact.cpp (Rt-MAF pairwise approach) in that it operates on
// ALL forests at once rather than eliminating them one by one against F1.
//
// Algorithm (applied at every recursive step):
//   1. Singleton sync (free): if any forest has a singleton leaf, isolate
//      that leaf in every other forest.
//   2. Pick sibling pair (a,b) from F1 = forests[0].
//      * If (a,b) are siblings in EVERY forest -> contract/shrink them in all
//        forests for free (they will definitely be together in the MAF).
//      * Otherwise -> 3-way branch:
//          A. Isolate a in all forests  (cost 1)
//          B. Isolate b in all forests  (cost 1)
//          C. Force a and b to be siblings in every forest by removing all
//             non-path edges on the a<->b path in each forest  (cost = total
//             edges removed).  Infeasible if a and b are in different
//             components in any forest.
//   3. Base case: F1 has no sibling pairs -> return F1->ord().

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <climits>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <functional>
#include <numeric>
#include <unordered_map>
#include <time.h>
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
#include "rootless_maf.h"

using namespace std;

volatile sig_atomic_t timeout_signaled = 0;
void signal_handler(int signal) { timeout_signaled = 1; }

// ---------- required globals ----------
bool DEFAULT_ALGORITHM    = true;
bool DEFAULT_OPTIMIZATIONS= true;
bool MULTI_APPROX   = false;
bool FPT            = false;
bool RF             = false;
bool QUIET          = true;
bool UNROOTED       = false;
bool ALL_UNROOTED   = false;
bool SIMPLE_UNROOTED      = false;
bool SIMPLE_UNROOTED_RSPR = false;
bool LCA_TEST       = false;
bool CLUSTER_TEST   = false;
bool TOTAL          = false;
bool PAIRWISE       = false;
bool PAIRWISE_SYMMETRIC   = true;
int  PAIRWISE_START       = 0;
int  PAIRWISE_END         = INT_MAX;
int  PAIRWISE_COL_START   = 0;
int  PAIRWISE_COL_END     = INT_MAX;
bool PAIRWISE_MAX   = false;
int  PAIRWISE_MAX_SPR     = INT_MAX;
bool APPROX         = false;
bool LOWER_BOUND    = false;
bool REDUCE_ONLY    = false;
bool PRINT_ROOTED_TREES   = false;
bool SHOW_MOVES     = false;
bool SEQUENCE       = false;
bool DEBUG_REVERSE  = false;
bool RANDOM_SPR     = false;
int  RANDOM_SPR_COUNT     = 0;
int  MULTI_TEST     = 0;

// ===========================================================================
//  main
// ===========================================================================
int main(int argc, char* argv[]) {
    for (int i = 1; i < argc; i++) {
        string arg(argv[i]);
        if (arg == "-v" || arg == "--verbose") VERBOSE = true;
    }

    cout << "# Hello, World!" << endl;

    signal(SIGTERM, signal_handler);
    signal(SIGINT,  signal_handler);

    CUT_ONE_B = REVERSE_CUT_ONE_B = REVERSE_CUT_ONE_B_3 = CUT_TWO_B = true;
    APPROX_CUT_ONE_B = APPROX_REVERSE_CUT_ONE_B = true;
    CLUSTER_TUNE = 30;

    srand(unsigned(time(0)));

    map<string, int> label_map;
    map<int, string> reverse_label_map;

    Input2 inputData = readInput();
    int m = inputData.numTrees;
    int n = inputData.numLeaves;

    if (m < 2) {
        cerr << "# Error: need at least 2 trees." << endl;
        return 1;
    }

    vector<Node*> trees;
    for (int i = 0; i < m && i < (int)inputData.trees.size(); i++)
        trees.push_back(build_tree(inputData.trees[i]));

    if (trees.empty()) { cerr << "# Error: no trees parsed." << endl; return 1; }

    for (Node* t : trees)
        t->labels_to_numbers(&label_map, &reverse_label_map);

    if (n == 0) n = (int)label_map.size();
    cout << "# num_trees=" << trees.size() << " num_leaves=" << n << endl;

    // Build template forests directly from the parsed trees.  No root
    // marker is added: an agreement-forest component is not required to
    // stay anchored to any tree's actual root, in any tree, so there is no
    // extra bookkeeping leaf needed to track "the" root-containing piece.
    vector<Forest*> templates;
    for (Node* t : trees) templates.push_back(new Forest(t));

    vector<Forest*> result_forests;
    int dist = rm_solve(templates, n, &result_forests);

    if (dist <= n) {
        if (!result_forests.empty()) {
            Forest* rf = result_forests[0];
            rf->numbers_to_labels(&reverse_label_map);
            rf->print_components();
        }
        cout << "# drSPR=" << dist << endl;
    } else {
        cout << "# drSPR=?" << endl;
    }

    for (Forest* f : result_forests) delete f;
    for (Forest* f : templates) delete f;
    for (Node* t : trees) t->delete_tree();
    return 0;
}
