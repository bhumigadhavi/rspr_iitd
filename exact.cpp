// exact.cpp — Rt-MAF algorithm for multiple rooted binary phylogenetic trees
// Based on: Shi, Wang, Chen, Feng, Guo
//   "Algorithms for parameterized maximum agreement forest problem on multiple trees"
//   Theoretical Computer Science 554, 2014, pp. 207-216
//
// Optimization: UndoMachine replaces copy_forests.
//   Previously every branch in Cases 1 and 3 deep-copied all m forests — O(n·m).
//   Now we mutate in-place and push undo events, then call um.undo_to(checkpoint)
//   to backtrack.  O(1) amortized bookkeeping per cut instead of O(n·m) copy.

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

typedef map<string, set<string>>         CombMap;
typedef map<string, pair<string,string>> ShrinkMap;  // combined_name → (lchild, rchild)

static set<string> expand_label(const string& lbl, const CombMap& cmap) {
    auto it = cmap.find(lbl);
    return (it == cmap.end()) ? set<string>{lbl} : it->second;
}

// Build a name→Node map from all leaves of a forest in one pass.
static unordered_map<string, Node*> build_leaf_map(Forest* f) {
    unordered_map<string, Node*> m;
    for (Node* lf : f->forest_leaves()) m[lf->str()] = lf;
    return m;
}

// Wrap a forest's component-0 in a new root with "p" as a leaf sibling.
// This ρ-extension enforces the rooted constraint when rt_maf_hlpr processes
// the resulting forest: any component that must sit "above the root" of some
// tree now requires an explicit root-cut, adding 1 to the drSPR count.
static Forest* make_rho_forest(Forest* f) {
    Node* comp = f->get_component(0);
    Node* rho_root = new Node("");
    rho_root->add_child(new Node("p"));
    rho_root->add_child(new Node(*comp));
    Forest* result = new Forest(rho_root);
    rho_root->delete_tree();
    return result;
}

// Returns true for the synthetic rho placeholder labels (>= 100000) used
// during cluster decomposition.  These must never appear in user output.
static bool is_rho_name(const string& nm) {
    if (nm.empty() || nm == "p") return false;
    if (nm.find_first_not_of("0123456789") != string::npos) return false;
    try { return stoi(nm) >= 100000; } catch (...) { return false; }
}

// Parse a combined-label string like "(a,(b,c))" and build a real Node
// subtree keeping only non-rho leaves.  Returns nullptr if nothing remains.
// Caller owns the returned Node*.
static Node* build_real_subtree(const string& name) {
    if (name.empty()) return nullptr;
    if (name[0] != '(') {
        if (is_rho_name(name) || name == "p") return nullptr;
        return new Node(name);
    }
    int depth = 0;
    size_t split = string::npos;
    for (size_t i = 1; i + 1 < name.size(); i++) {
        if (name[i] == '(')       depth++;
        else if (name[i] == ')')  depth--;
        else if (name[i] == ',' && depth == 0) { split = i; break; }
    }
    if (split == string::npos) return nullptr;
    Node* lchild = build_real_subtree(name.substr(1, split - 1));
    Node* rchild = build_real_subtree(name.substr(split + 1, name.size() - split - 2));
    if (!lchild && !rchild) return nullptr;
    if (!lchild) return rchild;
    if (!rchild) return lchild;
    Node* par = new Node("");
    par->add_child(lchild);
    par->add_child(rchild);
    return par;
}

// Prune rho placeholder leaves from a deep-copied component subtree in-place.
// parent is n's actual parent node (nullptr when n is a component root).
// Returns true if n itself was removed (was a rho leaf, or all children pruned).
static bool prune_rho_subtree(Node* n, Node* parent) {
    if (!n) return false;
    if (n->is_leaf()) {
        if (n->str() == "p" || is_rho_name(n->str())) {
            if (parent) n->cut_parent();
            delete n;
            return true;
        }
        return false;
    }
    vector<Node*> ch(n->get_children().begin(), n->get_children().end());
    for (Node* c : ch) prune_rho_subtree(c, n);
    size_t rem = n->get_children().size();
    if (rem == 0) {
        if (parent) n->cut_parent();
        delete n;
        return true;
    }
    if (rem == 1 && parent) {
        // Contract degree-1 node: bypass n, connect its only child to n's parent.
        Node* only = *n->get_children().begin();
        only->cut_parent();
        n->cut_parent();
        parent->add_child(only);
        delete n;
    }
    return false;
}

// Deep-copy comp, strip all rho placeholder leaves/labels, and add to dest.
// Handles: (1) rho singleton leaves, (2) combined-label leaves embedding rho
// names (from Rt-MAF shrinking), (3) internal nodes with rho leaf children.
static void add_component_filtered(Node* comp, Forest* dest) {
    if (!comp) return;
    Node* copy = new Node(*comp);
    if (copy->is_leaf()) {
        const string& nm = copy->str();
        if (nm == "p" || is_rho_name(nm)) { delete copy; return; }
        if (!nm.empty() && nm[0] == '(') {
            Node* real_tree = build_real_subtree(nm);
            delete copy;
            if (real_tree) dest->add_component(real_tree);
            return;
        }
        dest->add_component(copy);
        return;
    }
    {
        vector<Node*> ch(copy->get_children().begin(), copy->get_children().end());
        for (Node* c : ch) prune_rho_subtree(c, copy);
    }
    if (copy->get_children().empty()) { delete copy; return; }
    while (!copy->is_leaf() && copy->get_children().size() == 1) {
        Node* only = *copy->get_children().begin();
        only->cut_parent();
        delete copy;
        copy = only;
    }
    if (copy->is_leaf() && is_rho_name(copy->str())) { delete copy; return; }
    dest->add_component(copy);
}

// ===========================================================================
//  Undo event: re-insert a forest at a specific index in the forests vector.
//  Undo of "forests.erase(begin+idx)": re-insert f at idx.
// ===========================================================================
class EraseForestFromVector : public Undoable {
public:
    vector<Forest*>* vec;
    Forest* f;
    size_t idx;
    EraseForestFromVector(vector<Forest*>* v, size_t i)
        : vec(v), idx(i), f((*v)[i]) {}
    void undo() override {
        vec->insert(vec->begin() + idx, f);
    }
};

// ===========================================================================
//  Undo event: reverse the in-place expansion of a combined-label leaf.
//  undo() cuts each added child from root, deletes its subtree, restores name.
//  LIFO ordering guarantees inner expansions are undone before outer ones.
// ===========================================================================
class ExpandLeafUndo : public Undoable {
public:
    Node* root;
    string original_name;
    vector<Node*> added_children;
    ExpandLeafUndo(Node* r, const string& n, const vector<Node*>& ch)
        : root(r), original_name(n), added_children(ch) {}
    void undo() override {
        for (Node* ch : added_children) {
            ch->cut_parent();
            ch->delete_tree();
        }
        root->set_name(original_name);
    }
};

// ===========================================================================
//  Undoable make_singleton: cut leaf from parent, register as new component,
//  contract the now-degree-1 parent.  All changes logged in um.
// ===========================================================================
static void make_singleton_undoable(Node* leaf, Forest* F, UndoMachine& um) {
    Node* par = leaf->parent();
    if (!par) return;
    um.add_event(new CutParent(leaf));
    leaf->cut_parent();
    um.add_event(new AddComponent(F));
    F->add_component(leaf);
    ContractEvent(&um, par);
    par->contract(false);
}

// ===========================================================================
//  Undoable shrink: rename parent of sibling pair (a,b) to combined label,
//  cut both children.  Events: ChangeName(par), CutParent(a), CutParent(b).
//  Undo re-attaches children and restores par's name — no allocations needed.
// ===========================================================================
static string shrink_in_forest_undoable(Node* a, Node* b, CombMap& cmap,
                                         ShrinkMap& smap, UndoMachine& um) {
    string an = a->str(), bn = b->str();
    Node* par = a->parent();
    // Canonical order: smaller name left so combined labels are identical across forests.
    string lname = an < bn ? an : bn;
    string rname = an < bn ? bn : an;
    string cname = "(" + lname + "," + rname + ")";
    um.add_event(new ChangeName(par));
    par->set_name(cname);
    um.add_event(new CutParent(a));
    a->cut_parent();
    um.add_event(new CutParent(b));
    b->cut_parent();
    set<string> s = expand_label(an, cmap);
    for (const string& x : expand_label(bn, cmap)) s.insert(x);
    cmap[cname]  = s;
    smap[cname]  = {lname, rname};
    return cname;
}

// ===========================================================================
//  Undoable twin sync between F1 and F2.  Mirrors name_sync() but:
//   • twin-pointer changes are recorded via SetTwin events,
//   • unmatched F2 leaves are singletonised (undoable) instead of deleted.
// ===========================================================================
static void name_sync_undoable(Forest* F1, Forest* F2, const CombMap& cmap,
                                UndoMachine& um) {
    map<string, Node*> f2map;
    for (Node* lf : F2->forest_leaves()) f2map[lf->str()] = lf;

    // Record and clear all current twin pointers.
    for (Node* lf : F1->forest_leaves()) {
        um.add_event(new SetTwin(lf));
        lf->set_twin(nullptr);
    }
    for (auto& kv : f2map) {
        um.add_event(new SetTwin(kv.second));
        kv.second->set_twin(nullptr);
    }

    set<Node*> matched;

    for (Node* lf1 : F1->forest_leaves()) {
        const string& name = lf1->str();
        set<string> originals = expand_label(name, cmap);
        Node* best = nullptr;

        // Exact name match first (handles shrunk combined labels).
        auto it = f2map.find(name);
        if (it != f2map.end()) best = it->second;

        if (!best && originals.size() > 1) {
            // Combined label: match constituent individual leaves in F2.
            vector<Node*> found;
            for (const string& o : originals) {
                auto jt = f2map.find(o);
                if (jt != f2map.end()) found.push_back(jt->second);
            }
            if (!found.empty()) {
                bool same = true;
                for (size_t i = 1; i < found.size(); i++)
                    if (!found[0]->same_component(found[i])) { same = false; break; }
                if (same && found.size() > 1) {
                    Node* lca = found[0];
                    for (size_t i = 1; i < found.size(); i++)
                        lca = find_lca(lca, found[i]);
                    best = lca ? lca : found[0];
                } else {
                    best = found[0];
                }
                for (Node* nd : found) {
                    um.add_event(new SetTwin(nd));
                    nd->set_twin(lf1);
                    matched.insert(nd);
                }
            }
        }

        if (best) {
            um.add_event(new SetTwin(lf1));
            lf1->set_twin(best);
            um.add_event(new SetTwin(best));
            best->set_twin(lf1);
            matched.insert(best);
        }
    }

    // Singletonise unmatched F2 leaves instead of deleting them.
    // Snapshot the leaf list first (forest_leaves() is a fresh allocation).
    vector<Node*> f2leaves = F2->forest_leaves();
    for (Node* lf2 : f2leaves) {
        if (!matched.count(lf2) && lf2->str() != "p") {
            if (lf2->parent())
                make_singleton_undoable(lf2, F2, um);
            // If already a singleton root, twin was cleared above — nothing more.
        }
    }
}

// ===========================================================================
//  Undoable compute_maximal_af: removes conflicting sibling pairs from F1
//  until no more conflicts exist.  Each removal recorded via
//  make_singleton_undoable so the parent branch can undo them.
// ===========================================================================
static void compute_maximal_af_undoable(Forest* F1, Forest* F2,
                                         int* cuts, UndoMachine& um) {
    bool changed = true;
    while (changed) {
        changed = false;
        list<Node*>* sp = F1->find_sibling_pairs();
        for (auto it = sp->begin(); it != sp->end(); ) {
            Node* a = *it; ++it;
            Node* b = *it; ++it;
            Node* ta = a->get_twin();
            Node* tb = b->get_twin();
            if (!ta || !tb || !ta->same_component(tb)) {
                make_singleton_undoable(a, F1, um);
                if (cuts) (*cuts)++;
                changed = true;
                break;
            }
        }
        delete sp;
    }
}

// ===========================================================================
//  Undoable remove_non_path_edges: cuts off-path children of path nodes in
//  F1, each becoming a new component; contracts degree-1 path nodes.
//  *k decremented per cut.
// ===========================================================================
static void remove_non_path_edges_undoable(Forest* F1, Node* a1, Node* b1,
                                            int* k, UndoMachine& um) {
    Node* lca = find_lca(a1, b1);
    if (!lca) return;

    set<Node*> on_path;
    for (Node* c = a1; c && c != lca; c = c->parent()) on_path.insert(c);
    for (Node* c = b1; c && c != lca; c = c->parent()) on_path.insert(c);
    on_path.insert(lca);

    // Precompute depths once; computing depth inside the comparator would be
    // O(depth) per comparison → O(n·depth·log n) total.
    map<Node*, int> depth_map;
    for (Node* nd : on_path) {
        int d = 0;
        for (Node* c = nd; c; c = c->parent()) d++;
        depth_map[nd] = d;
    }

    // Sort deepest-first so we cut leaves before their parents.
    vector<Node*> path_vec(on_path.begin(), on_path.end());
    sort(path_vec.begin(), path_vec.end(), [&](Node* x, Node* y) {
        return depth_map[x] > depth_map[y];
    });

    for (Node* node : path_vec) {
        list<Node*> ch(node->get_children().begin(), node->get_children().end());
        for (Node* child : ch) {
            if (!on_path.count(child)) {
                um.add_event(new CutParent(child));
                child->cut_parent();
                um.add_event(new AddComponent(F1));
                F1->add_component(child);
                if (k) (*k)--;
            }
        }
    }
    for (Node* node : path_vec) {
        if (node->get_children().size() == 1 && node->parent()) {
            ContractEvent(&um, node);
            node->contract(false);
        }
    }
}

// ===========================================================================
//  Undoable expand_combined_leaf: expands a "(a,b)"-named leaf back into a
//  proper internal node with children reconstructed via build_tree.
//  ExpandLeafUndo is pushed AFTER the mutation; LIFO ensures inner expansions
//  undo before outer ones.
// ===========================================================================
static void expand_node_undoable(Node* root, const ShrinkMap& smap, UndoMachine& um) {
    if (!root) return;

    if (root->is_leaf() && !root->str().empty() && root->str()[0] == '(') {
        string orig_name = root->str();
        vector<Node*> added;

        auto it = smap.find(orig_name);
        if (it != smap.end()) {
            // Fast path: directly construct children from recorded pair names.
            // add_child() sets depth = parent->depth+1, so no fix_depths() needed.
            Node* lchild = new Node(it->second.first);
            Node* rchild = new Node(it->second.second);
            root->add_child(lchild);
            root->add_child(rchild);
            root->set_name("");
            added = {lchild, rchild};
        } else {
            // Fallback: combined label not in smap (e.g. from external shrink).
            Node* subtree = build_tree(orig_name);
            if (!subtree) return;
            list<Node*> sub_ch(subtree->get_children().begin(),
                               subtree->get_children().end());
            for (Node* child : sub_ch) {
                child->cut_parent();
                root->add_child(child);
                added.push_back(child);
            }
            root->set_name("");
            delete subtree;
            root->fix_depths();
        }

        // Push undo AFTER mutation — LIFO ensures inner expansions undo first.
        um.add_event(new ExpandLeafUndo(root, orig_name, added));

        // Recurse: children may carry nested combined labels.
        for (Node* child : added)
            expand_node_undoable(child, smap, um);
    } else {
        list<Node*> ch(root->get_children().begin(), root->get_children().end());
        for (Node* child : ch)
            expand_node_undoable(child, smap, um);
    }
}

static void expand_forest_undoable(Forest* F, const ShrinkMap& smap, UndoMachine& um) {
    for (Node* comp : F->components)
        if (comp) expand_node_undoable(comp, smap, um);
}

// ===========================================================================
//  Core Rt-MAF recursive helper — Figure 3 of Shi et al. 2014.
//
//  forests: shared by reference; ALL mutations are recorded in um so that
//           the caller can call um.undo_to(checkpoint) to backtrack.
//  cmap:    passed by value — each call owns its own copy; shrink entries
//           are naturally discarded when the call returns.
//  out:     if non-null, filled with a deep-copy of the winning forest state
//           at the base case, before any parent undo rolls back.
//
//  Returns the drSPR distance if achievable within budget k, else total_k+1.
// ===========================================================================
static int rt_maf_hlpr(vector<Forest*>& forests, int k, int total_k,
                        CombMap cmap, ShrinkMap smap, UndoMachine& um,
                        vector<Forest*>* out = nullptr) {
    if (timeout_signaled) return total_k + 1;

    // Step 1: base case — only one forest left.
    if (forests.size() == 1) {
        int ord = forests[0]->ord();
        if (ord <= total_k) {
            if (out) {
                out->clear();
                for (Forest* f : forests) out->push_back(new Forest(f));
            }
            return ord;
        }
        return total_k + 1;
    }

    Forest* F1 = forests[0];
    Forest* F2 = forests[1];

    // Step 2: prune — F1 already fragmented beyond total budget.
    if (F1->ord() > total_k) return total_k + 1;

    // Step 3: sync singletons — FREE, no k decrement (Theorem 3.2).
    {
        list<Node*> s2 = F2->find_singletons();
        for (Node* s : s2) {
            Node* t = s->get_twin();
            if (t && t->parent()) make_singleton_undoable(t, F1, um);
        }
        list<Node*> s1 = F1->find_singletons();
        for (Node* s : s1) {
            Node* t = s->get_twin();
            if (t && t->parent()) make_singleton_undoable(t, F2, um);
        }
    }
    if (F1->ord() > total_k) return total_k + 1;

    // Steps 4+5: find a sibling pair in F2; also compute a lower bound.
    Node* a2 = nullptr, *b2 = nullptr;
    {
        list<Node*>* sp2 = F2->find_sibling_pairs();

        if (sp2->empty()) {
            delete sp2;

            // No sibling pairs in F2.  Find remaining F1 conflicts w.r.t. F2.
            // For m=2 we can cut deterministically; for m>2 the choice of which
            // pair member to cut affects F3,...,Fm so we must branch.
            {
                list<Node*>* sp1 = F1->find_sibling_pairs();
                Node* conf_a = nullptr, *conf_b = nullptr;
                // First pass: singletonise free (null-twin) pairs w/o branching.
                bool changed = true;
                while (changed) {
                    changed = false;
                    delete sp1;
                    sp1 = F1->find_sibling_pairs();
                    for (auto it = sp1->begin(); it != sp1->end(); ) {
                        Node* a = *it; ++it;
                        Node* b = *it; ++it;
                        Node* ta = a->get_twin(), *tb = b->get_twin();
                        if (!ta || !tb) {
                            // one twin missing — free cut (leaf absent in F2)
                            make_singleton_undoable(a, F1, um);
                            k--;
                            if (k < 0) { delete sp1; return total_k + 1; }
                            changed = true; break;
                        }
                    }
                }
                // Find first costly conflict (both twins exist, diff components).
                for (auto it = sp1->begin(); it != sp1->end(); ) {
                    Node* a = *it; ++it;
                    Node* b = *it; ++it;
                    Node* ta = a->get_twin(), *tb = b->get_twin();
                    if (ta && tb && !ta->same_component(tb)) {
                        conf_a = a; conf_b = b; break;
                    }
                }
                delete sp1;

                if (conf_a) {
                    // Branch: cut conf_a OR conf_b from F1 (costs 1 each).
                    if (k <= 0) return total_k + 1;
                    string an = conf_a->str(), bn = conf_b->str();
                    int best = total_k + 1;
                    vector<Forest*> best_out;

                    // Branch A: cut a
                    {
                        int cp = um.num_events();
                        CombMap cb = cmap;
                        Node* xa = F1->find_leaf(an);
                        if (xa && xa->parent()) make_singleton_undoable(xa, F1, um);
                        vector<Forest*> br_out;
                        int r = rt_maf_hlpr(forests, k - 1, total_k, cb, smap, um,
                                            out ? &br_out : nullptr);
                        if (r < best) { best = r; best_out = move(br_out); }
                        else { for (Forest* f : br_out) delete f; }
                        um.undo_to(cp);
                    }
                    // Branch B: cut b
                    {
                        int cp = um.num_events();
                        CombMap cb = cmap;
                        Node* xb = F1->find_leaf(bn);
                        if (xb && xb->parent()) make_singleton_undoable(xb, F1, um);
                        vector<Forest*> br_out;
                        int r = rt_maf_hlpr(forests, k - 1, total_k, cb, smap, um,
                                            out ? &br_out : nullptr);
                        if (r < best) { best = r; best_out = move(br_out); }
                        else { for (Forest* f : br_out) delete f; }
                        um.undo_to(cp);
                    }
                    if (out) *out = move(best_out);
                    return best;
                }
            }

            // No conflicts remain — remove F2.
            um.add_event(new EraseForestFromVector(&forests, 1));
            forests.erase(forests.begin() + 1);

            if (forests.size() == 1) {
                int ord = forests[0]->ord();
                if (ord <= total_k) {
                    if (out) {
                        out->clear();
                        for (Forest* f : forests) out->push_back(new Forest(f));
                    }
                    return ord;
                }
                return total_k + 1;
            }

            // Expand combined-label leaves and re-sync twins for the next pair.
            for (Forest* f : forests)
                expand_forest_undoable(f, smap, um);
            name_sync_undoable(forests[0], forests[1], cmap, um);
            if (VERBOSE)
                fprintf(stderr,"# dbg F2removed: forests=%zu F1.ord=%d k=%d total_k=%d\n",
                        forests.size(), forests[0]->ord(), k, total_k);

            return rt_maf_hlpr(forests, k, total_k, cmap, smap, um, out);
        }

        // Scan all pairs: prefer a Case-2 (free) pair over a branching pair.
        // Picking Case-2 first eliminates free reductions before any branch is
        // taken, shrinking the search space at zero cost.
        int lb_conflicts = 0;
        Node* free_a = nullptr, *free_b = nullptr;
        for (auto it = sp2->begin(); it != sp2->end(); ) {
            Node* a = *it; ++it;
            Node* b = *it; ++it;
            if (!a2) { a2 = a; b2 = b; }   // default: first pair
            Node* ta = a->get_twin();
            Node* tb = b->get_twin();
            bool case2    = ta && tb && ta != tb
                            && ta->parent() && ta->parent() == tb->parent();
            bool combined = (ta && ta == tb);
            if ((case2 || combined) && !free_a) { free_a = a; free_b = b; }
            // Only count pairs where both twins exist: null-twin pairs are
            // handled for free (singletonised without decrementing k) and
            // must not inflate the lower bound.
            if (!case2 && !combined && ta && tb) lb_conflicts++;
        }
        delete sp2;
        // Prefer the free pair so the recursion shrinks without branching.
        if (free_a) { a2 = free_a; b2 = free_b; }

        if (k < lb_conflicts) {
            if (VERBOSE)
                fprintf(stderr,"# dbg lb_prune: k=%d lb=%d F1.ord=%d forests=%zu\n",
                        k, lb_conflicts, forests[0]->ord(), forests.size());
            return total_k + 1;
        }
    }
    if (!a2 || !b2) return total_k + 1;

    Node* a1 = a2->get_twin();
    Node* b1 = b2->get_twin();

    // Guard: twin missing — make the orphaned F2 node a singleton.
    if (!a1) {
        if (a2->parent()) make_singleton_undoable(a2, F2, um);
        return rt_maf_hlpr(forests, k, total_k, cmap, smap, um, out);
    }
    if (!b1) {
        if (b2->parent()) make_singleton_undoable(b2, F2, um);
        return rt_maf_hlpr(forests, k, total_k, cmap, smap, um, out);
    }

    bool diff_comp  = !a1->same_component(b1);
    bool sibs_in_F1 = (a1 != b1) && a1->parent() && (a1->parent() == b1->parent());

    // ------ Case 1: different components in F1 — 2-way branch ------
    if (diff_comp) {
        string a1n = a1->str(), b1n = b1->str();
        string a2n = a2->str(), b2n = b2->str();
        int best = total_k + 1;
        vector<Forest*> best_out;

        // Branch A: cut a in F1, F2, and all extra forests.
        {
            int cp = um.num_events();
            CombMap cb = cmap;
            Node* xa1 = F1->find_leaf(a1n);
            Node* xa2 = F2->find_leaf(a2n);
            if (xa1 && xa1->parent()) make_singleton_undoable(xa1, F1, um);
            if (xa2 && xa2->parent()) make_singleton_undoable(xa2, F2, um);
            for (size_t i = 2; i < forests.size(); i++) {
                Node* xi = forests[i]->find_leaf(a1n);
                if (!xi) xi = forests[i]->find_leaf(a2n);
                if (xi && xi->parent()) make_singleton_undoable(xi, forests[i], um);
            }
            vector<Forest*> branch_out;
            int r = rt_maf_hlpr(forests, k - 1, total_k, cb, smap, um,
                                 out ? &branch_out : nullptr);
            if (r < best) { best = r; best_out = move(branch_out); }
            else { for (Forest* f : branch_out) delete f; }
            um.undo_to(cp);
        }
        // Branch B: cut b in F1, F2, and all extra forests.
        {
            int cp = um.num_events();
            CombMap cb = cmap;
            Node* xb1 = F1->find_leaf(b1n);
            Node* xb2 = F2->find_leaf(b2n);
            if (xb1 && xb1->parent()) make_singleton_undoable(xb1, F1, um);
            if (xb2 && xb2->parent()) make_singleton_undoable(xb2, F2, um);
            for (size_t i = 2; i < forests.size(); i++) {
                Node* xi = forests[i]->find_leaf(b1n);
                if (!xi) xi = forests[i]->find_leaf(b2n);
                if (xi && xi->parent()) make_singleton_undoable(xi, forests[i], um);
            }
            vector<Forest*> branch_out;
            int r = rt_maf_hlpr(forests, k - 1, total_k, cb, smap, um,
                                 out ? &branch_out : nullptr);
            if (r < best) { best = r; best_out = move(branch_out); }
            else { for (Forest* f : branch_out) delete f; }
            um.undo_to(cp);
        }

        if (out) *out = move(best_out);
        return best;
    }

    // ------ Case 2: siblings in F1 (or already a combined leaf) — shrink ------
    if (sibs_in_F1 || a1 == b1) {
        if (a1 == b1) {
            // F1 already contracted these two; only shrink in F2.
            shrink_in_forest_undoable(a2, b2, cmap, smap, um);
        } else {
            string a1n = a1->str(), b1n = b1->str();
            shrink_in_forest_undoable(a1, b1, cmap, smap, um);
            shrink_in_forest_undoable(a2, b2, cmap, smap, um);
            for (size_t i = 2; i < forests.size(); i++) {
                // Build leaf map once per extra forest (2 lookups per forest
                // vs 2 full traversals with find_leaf).
                auto lmap = build_leaf_map(forests[i]);
                auto ita = lmap.find(a1n), itb = lmap.find(b1n);
                Node* ai = (ita != lmap.end()) ? ita->second : nullptr;
                Node* bi = (itb != lmap.end()) ? itb->second : nullptr;
                if (ai && bi && ai->parent() && ai->parent() == bi->parent())
                    shrink_in_forest_undoable(ai, bi, cmap, smap, um);
            }
        }
        name_sync_undoable(F1, F2, cmap, um);
        return rt_maf_hlpr(forests, k, total_k, cmap, smap, um, out);
    }

    // ------ Case 3: same component, not siblings — 3-way branch ------
    {
        string a1n = a1->str(), b1n = b1->str();
        string a2n = a2->str(), b2n = b2->str();
        int best = total_k + 1;
        vector<Forest*> best_out;

        // Branch A: cut a in F1, F2, and all extra forests.
        {
            int cp = um.num_events();
            CombMap cb = cmap;
            Node* xa1 = F1->find_leaf(a1n);
            Node* xa2 = F2->find_leaf(a2n);
            if (xa1 && xa1->parent()) make_singleton_undoable(xa1, F1, um);
            if (xa2 && xa2->parent()) make_singleton_undoable(xa2, F2, um);
            for (size_t i = 2; i < forests.size(); i++) {
                Node* xi = forests[i]->find_leaf(a1n);
                if (!xi) xi = forests[i]->find_leaf(a2n);
                if (xi && xi->parent()) make_singleton_undoable(xi, forests[i], um);
            }
            vector<Forest*> branch_out;
            int r = rt_maf_hlpr(forests, k - 1, total_k, cb, smap, um,
                                 out ? &branch_out : nullptr);
            if (r < best) { best = r; best_out = move(branch_out); }
            else { for (Forest* f : branch_out) delete f; }
            um.undo_to(cp);
        }
        // Branch B: cut b in F1, F2, and all extra forests.
        {
            int cp = um.num_events();
            CombMap cb = cmap;
            Node* xb1 = F1->find_leaf(b1n);
            Node* xb2 = F2->find_leaf(b2n);
            if (xb1 && xb1->parent()) make_singleton_undoable(xb1, F1, um);
            if (xb2 && xb2->parent()) make_singleton_undoable(xb2, F2, um);
            for (size_t i = 2; i < forests.size(); i++) {
                Node* xi = forests[i]->find_leaf(b1n);
                if (!xi) xi = forests[i]->find_leaf(b2n);
                if (xi && xi->parent()) make_singleton_undoable(xi, forests[i], um);
            }
            vector<Forest*> branch_out;
            int r = rt_maf_hlpr(forests, k - 1, total_k, cb, smap, um,
                                 out ? &branch_out : nullptr);
            if (r < best) { best = r; best_out = move(branch_out); }
            else { for (Forest* f : branch_out) delete f; }
            um.undo_to(cp);
        }
        // Branch C: remove all edges not on path a1→b1 in F1.
        {
            int cp = um.num_events();
            CombMap cb = cmap;
            // Two lookups on the same forest: build the map once.
            auto f1map = build_leaf_map(F1);
            auto _it_a1 = f1map.find(a1n), _it_b1 = f1map.find(b1n);
            Node* xa1 = (_it_a1 != f1map.end()) ? _it_a1->second : nullptr;
            Node* xb1 = (_it_b1 != f1map.end()) ? _it_b1->second : nullptr;
            if (xa1 && xb1) {
                int rem = k;
                remove_non_path_edges_undoable(F1, xa1, xb1, &rem, um);
                if (rem >= 0) {
                    name_sync_undoable(F1, F2, cb, um);
                    vector<Forest*> branch_out;
                    int r = rt_maf_hlpr(forests, rem, total_k, cb, smap, um,
                                        out ? &branch_out : nullptr);
                    if (r < best) { best = r; best_out = move(branch_out); }
                    else { for (Forest* f : branch_out) delete f; }
                }
            }
            um.undo_to(cp);
        }

        if (out) *out = move(best_out);
        return best;
    }
}

// ===========================================================================
//  Multi-tree cluster decomposition  (m > 2)
//
//  A leaf set L is a "common cluster" of T1,...,Tm when, for every Ti,
//  the leaves of Ti that fall in L form an exact clade of Ti.  The
//  generalised Bordewich–Semple additivity theorem then gives:
//
//      drSPR(T1,...,Tm) = drSPR(T1[L],...,Tm[L]) + drSPR(T1\L,...,Tm\L)
//
//  so each common cluster becomes an independent sub-problem.  On real
//  biological datasets that share many clades across all trees, this can
//  give an exponential speedup over the undecomposed branch-and-bound.
// ===========================================================================

// Collect non-rho leaf names in the subtree rooted at n.
static void collect_leaf_names(Node* n, set<string>& out) {
    if (!n) return;
    if (n->is_leaf()) {
        if (n->str() != "p") out.insert(n->str());
        return;
    }
    for (Node* c : n->get_children()) collect_leaf_names(c, out);
}

// In forest F find the root of the clade whose leaf set equals exactly the
// intersection of F's leaves with lset.  Returns nullptr when:
//   • no F-leaf is in lset;
//   • the matching leaves span multiple components;
//   • their LCA carries extra leaves (not all in lset); or
//   • the LCA has no parent (whole-tree / trivial cluster).
static Node* find_exact_clade_root(Forest* F, const set<string>& lset) {
    vector<Node*> matching;
    for (Node* lf : F->forest_leaves())
        if (lf->str() != "p" && lset.count(lf->str()))
            matching.push_back(lf);

    if (matching.empty()) return nullptr;

    for (size_t i = 1; i < matching.size(); i++)
        if (!matching[0]->same_component(matching[i])) return nullptr;

    Node* lca = matching[0];
    for (size_t i = 1; i < matching.size(); i++) {
        lca = find_lca(lca, matching[i]);
        if (!lca) return nullptr;
    }

    if (!lca->parent()) return nullptr;   // trivial – covers the whole tree

    // Exact: the LCA must not carry leaves outside lset.
    vector<Node*> lca_lv = lca->find_leaves();
    if (lca_lv.size() != matching.size()) return nullptr;

    return lca;
}

// Post-order scan of forests[0] component 0.  Returns every internal node n
// (with a parent and ≥2 children) whose leaf set forms an exact common clade
// in EVERY forest in the vector.
static vector<pair<Node*, vector<Node*>>>
find_common_cluster_points(vector<Forest*>& forests) {
    vector<pair<Node*, vector<Node*>>> result;
    if (forests.size() < 2 || forests[0]->num_components() == 0) return result;

    Node* root0 = forests[0]->get_component(0);
    if (!root0) return result;

    vector<Node*> postorder;
    function<void(Node*)> po = [&](Node* n) {
        for (Node* c : n->get_children()) po(c);
        postorder.push_back(n);
    };
    po(root0);

    for (Node* n : postorder) {
        if (n->is_leaf() || !n->parent()) continue;
        if ((int)n->get_children().size() < 2) continue;

        set<string> lset;
        collect_leaf_names(n, lset);
        if (lset.size() < 2) continue;   // singleton clusters give no speedup

        vector<Node*> clade_roots;
        bool ok = true;
        for (size_t i = 1; i < forests.size(); i++) {
            Node* cr = find_exact_clade_root(forests[i], lset);
            if (!cr) { ok = false; break; }
            clade_roots.push_back(cr);
        }
        if (ok) {
            if (VERBOSE) {
                set<string> ls2; collect_leaf_names(n, ls2);
                fprintf(stderr,"# dbg cpoint found: {");
                bool ff=true; for(auto& s:ls2){if(!ff)fprintf(stderr,",");fprintf(stderr,"%s",s.c_str());ff=false;}
                fprintf(stderr,"}\n");
            }
            result.push_back({n, clade_roots});
        }
    }
    return result;
}

// Cut common clusters out of forests in-place.
// Returns a list of independent sub-problem forest-vectors.
// forests itself becomes the remainder sub-problem (re-sync'd before return).
// Only top-level (non-nested) clusters are extracted; nested ones are handled
// recursively when each returned sub-problem is subsequently solved.
static list<vector<Forest*>> cluster_reduction_multi(vector<Forest*>& forests) {
    list<vector<Forest*>> sub_problems;
    if (forests.size() < 2) return sub_problems;

    auto cpoints = find_common_cluster_points(forests);
    if (VERBOSE) fprintf(stderr,"# dbg cpoints.size=%zu\n", cpoints.size());
    if (cpoints.empty()) return sub_problems;

    // Deepest first: inner nodes are cut before their ancestors are processed,
    // keeping ancestor pointers valid throughout.
    sort(cpoints.begin(), cpoints.end(),
         [](const pair<Node*, vector<Node*>>& a,
            const pair<Node*, vector<Node*>>& b) {
             return a.first->get_depth() > b.first->get_depth();
         });

    // Pre-collect every component root so we never free a node that is still
    // serving as a forest component — it may have contracted to a leaf singleton.
    set<Node*> comp_roots;
    for (Forest* f : forests)
        for (Node* comp : f->components)
            if (comp) comp_roots.insert(comp);

    set<Node*> used_f1_nodes;
    // Each extracted cluster is replaced by a unique placeholder leaf.
    // The counter is static so nested solve_sub_problem calls never reuse a
    // label, preventing name collisions between different recursion levels.
    static int rho_label = 100000;

    for (size_t _ci = 0; _ci < cpoints.size(); _ci++) {
        Node* n1 = cpoints[_ci].first;
        vector<Node*>& clade_roots = cpoints[_ci].second;
        // Skip nodes nested inside an already-extracted cluster.
        bool nested = false;
        for (Node* p = n1->parent(); p; p = p->parent())
            if (used_f1_nodes.count(p)) { nested = true; break; }
        if (nested) continue;

        // Verify n1 and all clade roots are still attached.
        // A prior cut can cause a root contraction that absorbs n1 into its
        // parent, leaving n1 with no parent and no children — skip if so.
        if (!n1->parent()) continue;
        bool valid = true;
        for (Node* cr : clade_roots)
            if (!cr->parent()) { valid = false; break; }
        if (!valid) continue;

        used_f1_nodes.insert(n1);

        if (VERBOSE) {
            set<string> ls;
            collect_leaf_names(n1, ls);
            fprintf(stderr,"# dbg extract cluster: {");
            bool first = true;
            for (const string& s : ls) { if (!first) fprintf(stderr,","); fprintf(stderr,"%s",s.c_str()); first=false; }
            fprintf(stderr,"}\n");
        }

        // Cut the cluster clade and replace it with a placeholder rho leaf.
        // All forests for the same cluster use the same rho label so that
        // sync_twins pairs them as twins in the outer problem.
        string rho_name = to_string(rho_label);
        auto cut_clade = [&](Node* clade_root) {
            Node* par = clade_root->parent();
            clade_root->cut_parent();
            Node* rho = new Node(rho_name);
            par->add_child(rho);
        };

        vector<Forest*> sub;
        cut_clade(n1);
        sub.push_back(new Forest(vector<Node*>{n1}));

        for (size_t i = 1; i < forests.size(); i++) {
            cut_clade(clade_roots[i - 1]);
            sub.push_back(new Forest(vector<Node*>{clade_roots[i - 1]}));
        }

        rho_label++;
        sync_twins(sub[0], sub[1]);
        sub_problems.push_back(move(sub));
    }

    if (!sub_problems.empty())
        sync_twins(forests[0], forests[1]);

    return sub_problems;
}

// Forward declaration — defined immediately below.
static int solve_sub_problem(vector<Forest*>& templates, int max_k,
                              vector<Forest*>* out_forests = nullptr);

// Solve one sub-problem given as pre-built template forests.
// Recursively applies cluster decomposition; falls back to the Rt-MAF
// k-bounded branch-and-bound when no common clusters exist.
// Returns the drSPR distance, or (max_k+1) if infeasible within budget.
// If out_forests is non-null, *out_forests[0] receives the F1 MAF forest.
static int solve_sub_problem(vector<Forest*>& templates, int max_k,
                              vector<Forest*>* out_forests) {
    if (templates.empty() || max_k < 0) return 0;
    if (templates[0]->forest_leaves().empty()) return 0;

    // ---- Try cluster decomposition on fresh copies of the templates ----
    {
        vector<Forest*> ct;
        for (Forest* t : templates) ct.push_back(new Forest(t));
        sync_twins(ct[0], ct[1]);

        list<vector<Forest*>> subs_list = cluster_reduction_multi(ct);

        if (VERBOSE)
            fprintf(stderr,"# dbg clust: subs=%zu orig_leaves=%d rem_leaves=%d\n",
                    subs_list.size(), (int)templates[0]->forest_leaves().size(),
                    (int)ct[0]->forest_leaves().size());
        if (!subs_list.empty()) {
            // Convert to vector for indexed access and safe cleanup.
            vector<vector<Forest*>> subs(
                make_move_iterator(subs_list.begin()),
                make_move_iterator(subs_list.end()));

            int total = 0;
            bool feasible = true;

            // Collect F1 components from each sub-problem and the remainder.
            Forest* merged_f1 = out_forests ? new Forest() : nullptr;

            for (auto& sub : subs) {
                if (timeout_signaled) { feasible = false; break; }
                int budget = max_k - total;
                if (budget < 0) { feasible = false; break; }
                int sub_n = (int)sub[0]->forest_leaves().size();
                vector<Forest*> sub_out;
                int d = solve_sub_problem(sub, min(sub_n, budget),
                                          merged_f1 ? &sub_out : nullptr);
                if (d > budget) { feasible = false; break; }
                total += d;
                if (merged_f1 && !sub_out.empty()) {
                    for (Node* comp : sub_out[0]->components)
                        add_component_filtered(comp, merged_f1);
                    for (Forest* f : sub_out) delete f;
                }
            }

            // Solve the remainder (ct modified by cluster_reduction_multi).
            if (feasible) {
                vector<Node*> rem_leaves = ct[0]->forest_leaves();
                int rem_n = (int)rem_leaves.size();
                if (rem_n > 1) {
                    int budget = max_k - total;
                    vector<Forest*> rem_out;
                    int d = solve_sub_problem(ct, budget,
                                              merged_f1 ? &rem_out : nullptr);
                    if (d > budget) feasible = false;
                    else {
                        total += d;
                        if (merged_f1 && !rem_out.empty()) {
                            for (Node* comp : rem_out[0]->components)
                                add_component_filtered(comp, merged_f1);
                            for (Forest* f : rem_out) delete f;
                        }
                    }
                } else if (rem_n == 1 && merged_f1) {
                    // Single leaf left after cluster extraction (e.g. an outermost
                    // cherry leaf that was never part of any extracted cluster).
                    // It agrees trivially across all trees (drSPR += 0) but must
                    // appear in the MAF output as a singleton component.
                    for (Node* lf : rem_leaves)
                        if (lf->str() != "p" && !is_rho_name(lf->str()))
                            merged_f1->add_component(new Node(lf->str()));
                }
            }

            for (auto& sub : subs) for (Forest* f : sub) delete f;
            for (Forest* f : ct) delete f;

            if (feasible && out_forests) out_forests->push_back(merged_f1);
            else delete merged_f1;

            return feasible ? total : max_k + 1;
        }

        // No common clusters found — fall through to the k-loop.
        for (Forest* f : ct) delete f;
    }

    // ---- k-bounded Rt-MAF branch-and-bound ----

    // Lower bound: 3-approximation / 3.
    // Suppress stdout during approx to prevent internal warnings from polluting output.
    int lb = 0;
    {
        streambuf* orig = cout.rdbuf();
        ostringstream sink;
        cout.rdbuf(sink.rdbuf());
        Forest f0(*templates[0]), f1(*templates[1]);
        lb = max(0, rSPR_worse_3_mult_approx(&f0, &f1) / 3);
        cout.rdbuf(orig);
    }

    for (int k = lb; k <= max_k; k++) {
        if (timeout_signaled) return max_k + 1;

        // Build working copies from templates for this k-iteration.
        vector<Forest*> wf;
        for (Forest* t : templates) wf.push_back(new Forest(t));
        sync_twins(wf[0], wf[1]);

        CombMap cmap; ShrinkMap smap; UndoMachine um;
        vector<Forest*> k_out;
        int result = rt_maf_hlpr(wf, k, k, cmap, smap, um,
                                 out_forests ? &k_out : nullptr);
        um.undo_all();
        for (Forest* f : wf) delete f;

        if (result <= k) {
            if (out_forests) *out_forests = move(k_out);
            return result;
        }
        for (Forest* f : k_out) delete f;
    }
    return max_k + 1;
}

// ===========================================================================
//  Solve the multi-tree MAF.  Returns drSPR distance, or -1 on error.
// ===========================================================================
static int solve_multi_maf(vector<Node*>& trees, int n,
                            map<string,int>& label_map,
                            map<int,string>& reverse_label_map,
                            vector<Forest*>* out_forests = nullptr)
{
    int m = (int)trees.size();
    if (m < 2) return 0;

    // m == 2: delegate to the exact 2-tree algorithm with ClusterForest.
    if (m == 2) {
        Forest* out_F1 = nullptr;
        Forest* out_F2 = nullptr;
        int dist = rSPR_branch_and_bound_simple_clustering(
            trees[0], trees[1], false,
            &label_map, &reverse_label_map,
            -1, n - 1,
            out_forests ? &out_F1 : nullptr,
            out_forests ? &out_F2 : nullptr);
        if (timeout_signaled) {
            delete out_F1; delete out_F2;
            return -1;
        }
        if (out_forests && out_F1 && out_F2) {
            out_forests->push_back(out_F1);
            out_forests->push_back(out_F2);
        }
        return dist;
    }

    // m > 2: Rt-MAF with multi-tree cluster decomposition.
    CUT_ONE_B = REVERSE_CUT_ONE_B = REVERSE_CUT_ONE_B_3 = CUT_TWO_B = true;
    APPROX_CUT_ONE_B = APPROX_REVERSE_CUT_ONE_B = true;
    CLUSTER_TUNE = 30;

    // Build template forests then hand off to the cluster-aware solver.
    // solve_sub_problem recursively decomposes common clades and runs the
    // Rt-MAF k-loop on each independent piece.
    vector<Forest*> templates;
    for (Node* t : trees) templates.push_back(new Forest(t));

    // Step 1: compute distance + fallback MAF via cluster decomposition.
    // The cluster decomp can underestimate the true drSPR for m>2 trees when
    // its "remove-cluster" reduction misses structural differences, so we keep
    // the MAF it builds as a valid fallback even when the distance is wrong.
    vector<Forest*> fallback_maf;
    int answer = solve_sub_problem(templates, n,
                                   out_forests ? &fallback_maf : nullptr);

    // Step 2: find the true drSPR via rt_maf_hlpr, starting from the cluster-
    // decomp lower bound lb=answer.  The cluster decomp can underestimate for
    // m>2 trees, so we increase k until rt_maf_hlpr succeeds.  All (m-1)!
    // orderings of the non-T1 trees are tried per k level to maximise the
    // chance of finding the optimum with the sequential Rt-MAF scan.
    if (answer >= 0 && answer <= n) {
        int true_dist = -1;
        vector<Forest*> maf_out;

        // Heuristic ordering: sort T2…Tm by pairwise approx drSPR from T1
        // (ascending).  The most similar tree placed as T2 reduces branching.
        // At k == answer (lower bound) only this one ordering is tried;
        // at k > answer all (m-1)! permutations are tried as a fallback.
        vector<int> heuristic_perm(m - 1);
        iota(heuristic_perm.begin(), heuristic_perm.end(), 1);
        {
            vector<pair<int,int>> pair_dists;
            pair_dists.reserve(m - 1);
            streambuf* orig = cout.rdbuf();
            ostringstream sink;
            cout.rdbuf(sink.rdbuf());
            for (int i = 1; i < m; i++) {
                Forest f0(*templates[0]), fi(*templates[i]);
                int d = rSPR_worse_3_mult_approx(&f0, &fi);
                pair_dists.push_back({d, i});
            }
            cout.rdbuf(orig);
            sort(pair_dists.begin(), pair_dists.end());
            for (int i = 0; i < m - 1; i++) heuristic_perm[i] = pair_dists[i].second;
        }

        for (int k = answer; k <= n && !timeout_signaled; k++) {
            int best_result = k + 1;
            vector<Forest*> best_out;

            // Run rt_maf_hlpr for cur_perm and update best_result / best_out.
            auto try_perm = [&](const vector<int>& cur_perm) {
                // Use ρ-extended forests so rt_maf_hlpr enforces the rooted
                // constraint: a component cannot sit "above the root" of any
                // tree without an explicit root-cut costing 1.
                vector<Forest*> wf;
                wf.push_back(make_rho_forest(templates[0]));
                for (int idx : cur_perm) wf.push_back(make_rho_forest(templates[idx]));
                sync_twins(wf[0], wf[1]);

                CombMap cmap; ShrinkMap smap; UndoMachine um;
                vector<Forest*> k_out;
                int result = rt_maf_hlpr(wf, k, k, cmap, smap, um,
                                         out_forests ? &k_out : nullptr);
                um.undo_all();
                for (Forest* f : wf) delete f;

                if (result < best_result) {
                    best_result = result;
                    for (Forest* f : best_out) delete f;
                    best_out = move(k_out);
                } else {
                    for (Forest* f : k_out) delete f;
                }
            };

            // Always try the heuristic ordering first (O(1) per k-iteration).
            try_perm(heuristic_perm);

            // Fall back to all (m-1)! permutations only when the budget has
            // grown past the cluster-decomp lower bound.  This slow path is
            // rare; on most instances the heuristic ordering suffices at k==answer.
            if (best_result > k && k > answer) {
                vector<int> perm(m - 1);
                iota(perm.begin(), perm.end(), 1);
                do {
                    if (perm == heuristic_perm) continue;
                    try_perm(perm);
                    if (best_result <= k) break;
                } while (next_permutation(perm.begin(), perm.end()));
            }

            if (best_result <= k) {
                true_dist = best_result;
                maf_out = move(best_out);
                break;
            }
            for (Forest* f : best_out) delete f;
        }

        if (true_dist >= 0) {
            answer = true_dist;
            if (out_forests) {
                // Filter "p" (virtual ρ leaf) from rt_maf_hlpr output.
                Forest* filtered = new Forest();
                if (!maf_out.empty())
                    for (Node* comp : maf_out[0]->components)
                        add_component_filtered(comp, filtered);
                for (Forest* f : maf_out) delete f;
                *out_forests = {filtered};
            }
            for (Forest* f : fallback_maf) delete f;
        } else if (!timeout_signaled && out_forests) {
            // Budget exhausted normally — cluster-decomp fallback is a valid result.
            *out_forests = move(fallback_maf);
        } else {
            // Timed out without an exact answer — discard partial fallback.
            for (Forest* f : fallback_maf) delete f;
            if (timeout_signaled) answer = n + 1;  // force -1 return below
        }
    } else {
        for (Forest* f : fallback_maf) delete f;
    }

    for (Forest* f : templates) delete f;
    return (answer <= n) ? answer : -1;
}

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

    map<string,int> label_map;
    map<int,string> reverse_label_map;

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

    vector<Forest*> result_forests;
    int dist = solve_multi_maf(trees, n, label_map, reverse_label_map, &result_forests);

    if (dist >= 0) {
        if (!result_forests.empty()) {
            result_forests[0]->numbers_to_labels(&reverse_label_map);
            result_forests[0]->print_components();
        }
        cout << "# drSPR=" << dist << endl;
        for (Forest* f : result_forests) delete f;
    } else {
        cout << "# drSPR=?" << endl;
    }

    for (Node* t : trees) t->delete_tree();
    return 0;
}
