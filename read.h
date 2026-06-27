#include <iostream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

// Struct to hold the parsed input
struct Input2 {
    int numTrees = 0;       // 't' from the #p line
    int numLeaves = 0;      // 'n' from the #p line
    
    // Lower bound track parameters
    float a = 1.0f;         // 'a' from the #a line
    int b = 0;              // 'b' from the #a line
    
    // Tree decomposition parameters
    int treeWidth = -1;     // 'tw' from the #x treedecomp line (-1 means not provided)
    string treeDecompRaw = ""; // Stores the raw "[{tw},{bags},{edges}]" string
    
    string idigest = """";
    string name = "";
    vector<string> trees;
};

Input2 readInput() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string line;
    Input2 inputData;

    while (getline(cin, line)) {
        if (line.empty()) continue;  // ignore empty lines

        if (line[0] == '#') {
            if (line.size() > 2 && line[1] == 'p') {
                // parse "#p t n"
                stringstream ss(line.substr(2)); // skip "#p"
                ss >> inputData.numTrees >> inputData.numLeaves;
            }
            else if (line.size() > 2 && line[1] == 'a') {
                // parse "#a a b"
                stringstream ss(line.substr(2)); // skip "#a"
                ss >> inputData.a >> inputData.b;
            }
            else if (line.size() > 2 && line[1] == 'x') {
                // parse "#x key value"
                string prefix = "#x treedecomp ";
                if (line.substr(0, prefix.size()) == prefix) {
                    string value = line.substr(prefix.size());
                    inputData.treeDecompRaw = value;
                    
                    // Extract {tw} from the JSON subset "[tw,...]"
                    size_t bracketPos = value.find('[');
                    size_t commaPos = value.find(',');
                    if (bracketPos!= string::npos && commaPos!= string::npos && commaPos > bracketPos) {
                        string twStr = value.substr(bracketPos + 1, commaPos - bracketPos - 1);
                        try {
                            inputData.treeWidth = stoi(twStr);
                        } catch (...) {
                            inputData.treeWidth = -1; // Fallback if parsing fails
                        }
                    }
                }
            }
            else if (line.size() > 2 && line[1] == 's') {
                const string pfx_idigest = "#s idigest ";
                const string pfx_name    = "#s name ";
                if (line.substr(0, pfx_idigest.size()) == pfx_idigest)
                    inputData.idigest = line.substr(pfx_idigest.size());
                else if (line.substr(0, pfx_name.size()) == pfx_name)
                    inputData.name = line.substr(pfx_name.size());
            }
            continue; // ignore all other comment lines
        }

        // any other non-empty line is a tree in Newick format
        inputData.trees.push_back(line);
    }

    return inputData;
}