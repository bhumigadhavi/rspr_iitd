#include <iostream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

// Struct to hold the parsed input
struct Input2 {
    int numTrees;   // 't' from the #p line
    int numLeaves;  // 'n' from the #p line
    vector<string> trees;
};

Input2 readInput() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string line;
    int numTrees = 0, numLeaves = 0;
    vector<string> trees;
    while (getline(cin, line)) {
        if (line.empty()) continue;  // ignore empty lines
        if (line[0] == '#') {
            if (line.size() > 2 && line[1] == 'p') {
                // parse "#p t n"
                stringstream ss(line.substr(2)); // skip "#p"
                ss >> numTrees >> numLeaves;
            }
            continue; // ignore all comment lines
        }
        cout << "#Parsed line: " << line << "\n"; // Debug: print non-comment lines

        // any other non-empty line is a tree in Newick format
        trees.push_back(line);
    }

    // // CRITICAL: Print summary to cerr, NOT cout. 
    // // The checker will fail if anything other than the forest is printed to stdout.
    // cerr << "# Number of trees declared: " << numTrees << "\n";
    // cerr << "# Number of leaves per tree: " << numLeaves << "\n";
    // cerr << "# Actual trees read: " << trees.size() << "\n";
    
    // Optional safeguard if you only want exactly 2 trees
    // if (trees.size()!= 2) {
    //     cerr << " Expected exactly 2 trees, but read " << trees.size() << "\n";
    // }
    return {numTrees, numLeaves, trees};
}