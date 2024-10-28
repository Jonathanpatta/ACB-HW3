#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <algorithm>

using namespace std;

std::unordered_map<char, std::unordered_map<char, int>> blosum62 = {
    {'A', {
        {'A', 4}, {'R', -1}, {'N', -2}, {'D', -2}, {'C', 0}, {'Q', -1}, {'E', -1}, {'G', 0},  
        {'H', -2}, {'I', -1}, {'L', -1}, {'K', -1}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 1},
        {'T', 0}, {'W', -3}, {'Y', -2}, {'V', 0}, {'B', -2}, {'J', -1}, {'Z', -1}, {'X', -1}, 
        {'*', -4}
    }},
    {'R', {
        {'A', -1}, {'R', 5}, {'N', 0}, {'D', -2}, {'C', -3}, {'Q', 1}, {'E', 0}, {'G', -2},   
        {'H', 0}, {'I', -3}, {'L', -2}, {'K', 2}, {'M', -1}, {'F', -3}, {'P', -2}, {'S', -1}, 
        {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -3}, {'B', -1}, {'J', -2}, {'Z', 0}, {'X', -1},
        {'*', -4}
    }},
    {'N', {
        {'A', -2}, {'R', 0}, {'N', 6}, {'D', 1}, {'C', -3}, {'Q', 0}, {'E', 0}, {'G', 0},     
        {'H', 1}, {'I', -3}, {'L', -3}, {'K', 0}, {'M', -2}, {'F', -3}, {'P', -2}, {'S', 1},  
        {'T', 0}, {'W', -4}, {'Y', -2}, {'V', -3}, {'B', 4}, {'J', -3}, {'Z', 0}, {'X', -1},  
        {'*', -4}
    }},
    {'D', {
        {'A', -2}, {'R', -2}, {'N', 1}, {'D', 6}, {'C', -3}, {'Q', 0}, {'E', 2}, {'G', -1},
        {'H', -1}, {'I', -3}, {'L', -4}, {'K', -1}, {'M', -3}, {'F', -3}, {'P', -1}, {'S', 0},
        {'T', -1}, {'W', -4}, {'Y', -3}, {'V', -3}, {'B', 4}, {'J', -3}, {'Z', 1}, {'X', -1},
        {'*', -4}
    }},
    {'C', {
        {'A', 0}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', 9}, {'Q', -3}, {'E', -4}, {'G', -3},
        {'H', -3}, {'I', -1}, {'L', -1}, {'K', -3}, {'M', -1}, {'F', -2}, {'P', -3}, {'S', -1},
        {'T', -1}, {'W', -2}, {'Y', -2}, {'V', -1}, {'B', -3}, {'J', -1}, {'Z', -3}, {'X', -1},
        {'*', -4}
    }},
    {'Q', {
        {'A', -1}, {'R', 1}, {'N', 0}, {'D', 0}, {'C', -3}, {'Q', 5}, {'E', 2}, {'G', -2},
        {'H', 0}, {'I', -3}, {'L', -2}, {'K', 1}, {'M', 0}, {'F', -3}, {'P', -1}, {'S', 0},
        {'T', -1}, {'W', -2}, {'Y', -1}, {'V', -2}, {'B', 0}, {'J', -2}, {'Z', 4}, {'X', -1},
        {'*', -4}
    }},
    {'E', {
        {'A', -1}, {'R', 0}, {'N', 0}, {'D', 2}, {'C', -4}, {'Q', 2}, {'E', 5}, {'G', -2},
        {'H', 0}, {'I', -3}, {'L', -3}, {'K', 1}, {'M', -2}, {'F', -3}, {'P', -1}, {'S', 0},
        {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -2}, {'B', 1}, {'J', -3}, {'Z', 4}, {'X', -1},
        {'*', -4}
    }},
    {'G', {
        {'A', 0}, {'R', -2}, {'N', 0}, {'D', -1}, {'C', -3}, {'Q', -2}, {'E', -2}, {'G', 6},
        {'H', -2}, {'I', -4}, {'L', -4}, {'K', -2}, {'M', -3}, {'F', -3}, {'P', -2}, {'S', 0},
        {'T', -2}, {'W', -2}, {'Y', -3}, {'V', -3}, {'B', -1}, {'J', -4}, {'Z', -2}, {'X', -1},
        {'*', -4}
    }},
    {'H', {
        {'A', -2}, {'R', 0}, {'N', 1}, {'D', -1}, {'C', -3}, {'Q', 0}, {'E', 0}, {'G', -2},
        {'H', 8}, {'I', -3}, {'L', -3}, {'K', -1}, {'M', -2}, {'F', -1}, {'P', -2}, {'S', -1},
        {'T', -2}, {'W', -2}, {'Y', 2}, {'V', -3}, {'B', 0}, {'J', -3}, {'Z', 0}, {'X', -1},
        {'*', -4}
    }},
    {'I', {
        {'A', -1}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', -1}, {'Q', -3}, {'E', -3}, {'G', -4},
        {'H', -3}, {'I', 4}, {'L', 2}, {'K', -3}, {'M', 1}, {'F', 0}, {'P', -3}, {'S', -2},
        {'T', -1}, {'W', -3}, {'Y', -1}, {'V', 3}, {'B', -3}, {'J', 3}, {'Z', -3}, {'X', -1},
        {'*', -4}
    }},
    {'L', {
        {'A', -1}, {'R', -2}, {'N', -3}, {'D', -4}, {'C', -1}, {'Q', -2}, {'E', -3}, {'G', -4},
        {'H', -3}, {'I', 2}, {'L', 4}, {'K', -2}, {'M', 2}, {'F', 0}, {'P', -3}, {'S', -2},
        {'T', -1}, {'W', -2}, {'Y', -1}, {'V', 1}, {'B', -4}, {'J', 3}, {'Z', -3}, {'X', -1},
        {'*', -4}
    }},
    {'K', {
        {'A', -1}, {'R', 2}, {'N', 0}, {'D', -1}, {'C', -3}, {'Q', 1}, {'E', 1}, {'G', -2},
        {'H', -1}, {'I', -3}, {'L', -2}, {'K', 5}, {'M', -1}, {'F', -3}, {'P', -1}, {'S', 0},
        {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -2}, {'B', 0}, {'J', -3}, {'Z', 1}, {'X', -1},
        {'*', -4}
    }},
    {'M', {
        {'A', -1}, {'R', -1}, {'N', -2}, {'D', -3}, {'C', -1}, {'Q', 0}, {'E', -2}, {'G', -3},
        {'H', -2}, {'I', 1}, {'L', 2}, {'K', -1}, {'M', 5}, {'F', 0}, {'P', -2}, {'S', -1},
        {'T', -1}, {'W', -1}, {'Y', -1}, {'V', 1}, {'B', -3}, {'J', 2}, {'Z', -1}, {'X', -1},
        {'*', -4}
    }},
    {'F', {
        {'A', -2}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', -2}, {'Q', -3}, {'E', -3}, {'G', -3},
        {'H', -1}, {'I', 0}, {'L', 0}, {'K', -3}, {'M', 0}, {'F', 6}, {'P', -4}, {'S', -2},
        {'T', -2}, {'W', 1}, {'Y', 3}, {'V', -1}, {'B', -3}, {'J', 0}, {'Z', -3}, {'X', -1},
        {'*', -4}
    }},
    {'P', {
        {'A', -1}, {'R', -2}, {'N', -2}, {'D', -1}, {'C', -3}, {'Q', -1}, {'E', -1}, {'G', -2},
        {'H', -2}, {'I', -3}, {'L', -3}, {'K', -1}, {'M', -2}, {'F', -4}, {'P', 7}, {'S', -1},
        {'T', -1}, {'W', -4}, {'Y', -3}, {'V', -2}, {'B', -2}, {'J', -3}, {'Z', -1}, {'X', -1},
        {'*', -4}
    }},
    {'S', {
        {'A', 1}, {'R', -1}, {'N', 1}, {'D', 0}, {'C', -1}, {'Q', 0}, {'E', 0}, {'G', 0},
        {'H', -1}, {'I', -2}, {'L', -2}, {'K', 0}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 4},
        {'T', 1}, {'W', -3}, {'Y', -2}, {'V', -2}, {'B', 0}, {'J', -2}, {'Z', 0}, {'X', -1},
        {'*', -4}
    }},
    {'T', {
        {'A', 0}, {'R', -1}, {'N', 0}, {'D', -1}, {'C', -1}, {'Q', -1}, {'E', -1}, {'G', -2},
        {'H', -2}, {'I', -1}, {'L', -1}, {'K', -1}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 1},
        {'T', 5}, {'W', -2}, {'Y', -2}, {'V', 0}, {'B', -1}, {'J', -1}, {'Z', -1}, {'X', -1},
        {'*', -4}
    }},
    {'W', {
        {'A', -3}, {'R', -3}, {'N', -4}, {'D', -4}, {'C', -2}, {'Q', -2}, {'E', -3}, {'G', -2},
        {'H', -2}, {'I', -3}, {'L', -2}, {'K', -3}, {'M', -1}, {'F', 1}, {'P', -4}, {'S', -3},
        {'T', -2}, {'W', 11}, {'Y', 2}, {'V', -3}, {'B', -4}, {'J', -2}, {'Z', -2}, {'X', -1},
        {'*', -4}
    }},
    {'Y', {
        {'A', -2}, {'R', -2}, {'N', -2}, {'D', -3}, {'C', -2}, {'Q', -1}, {'E', -2}, {'G', -3},
        {'H', 2}, {'I', -1}, {'L', -1}, {'K', -2}, {'M', -1}, {'F', 3}, {'P', -3}, {'S', -2},
        {'T', -2}, {'W', 2}, {'Y', 7}, {'V', -1}, {'B', -3}, {'J', -1}, {'Z', -2}, {'X', -1},
        {'*', -4}
    }},
    {'V', {
        {'A', 0}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', -1}, {'Q', -2}, {'E', -2}, {'G', -3},
        {'H', -3}, {'I', 3}, {'L', 1}, {'K', -2}, {'M', 1}, {'F', -1}, {'P', -2}, {'S', -2},
        {'T', 0}, {'W', -3}, {'Y', -1}, {'V', 4}, {'B', -3}, {'J', 2}, {'Z', -2}, {'X', -1},
        {'*', -4}
    }},
    {'B', {
        {'A', -2}, {'R', -1}, {'N', 4}, {'D', 4}, {'C', -3}, {'Q', 0}, {'E', 1}, {'G', -1},
        {'H', 0}, {'I', -3}, {'L', -4}, {'K', 0}, {'M', -3}, {'F', -3}, {'P', -2}, {'S', 0},
        {'T', -1}, {'W', -4}, {'Y', -3}, {'V', -3}, {'B', 4}, {'J', -3}, {'Z', 0}, {'X', -1},
        {'*', -4}
    }},
    {'J', {
        {'A', -1}, {'R', -2}, {'N', -3}, {'D', -3}, {'C', -1}, {'Q', -2}, {'E', -3}, {'G', -4},
        {'H', -3}, {'I', 3}, {'L', 3}, {'K', -3}, {'M', 2}, {'F', 0}, {'P', -3}, {'S', -2},
        {'T', -1}, {'W', -2}, {'Y', -1}, {'V', 2}, {'B', -3}, {'J', 3}, {'Z', -3}, {'X', -1},
        {'*', -4}
    }},
    {'Z', {
        {'A', -1}, {'R', 0}, {'N', 0}, {'D', 1}, {'C', -3}, {'Q', 4}, {'E', 4}, {'G', -2},
        {'H', 0}, {'I', -3}, {'L', -3}, {'K', 1}, {'M', -1}, {'F', -3}, {'P', -1}, {'S', 0},
        {'T', -1}, {'W', -2}, {'Y', -2}, {'V', -2}, {'B', 0}, {'J', -3}, {'Z', 4}, {'X', -1},
        {'*', -4}
    }},
    {'X', {
        {'A', -1}, {'R', -1}, {'N', -1}, {'D', -1}, {'C', -1}, {'Q', -1}, {'E', -1}, {'G', -1},
        {'H', -1}, {'I', -1}, {'L', -1}, {'K', -1}, {'M', -1}, {'F', -1}, {'P', -1}, {'S', -1},
        {'T', -1}, {'W', -1}, {'Y', -1}, {'V', -1}, {'B', -1}, {'J', -1}, {'Z', -1}, {'X', -1},
        {'*', -4}
    }},
    {'*', {
        {'A', -4}, {'R', -4}, {'N', -4}, {'D', -4}, {'C', -4}, {'Q', -4}, {'E', -4}, {'G', -4},
        {'H', -4}, {'I', -4}, {'L', -4}, {'K', -4}, {'M', -4}, {'F', -4}, {'P', -4}, {'S', -4},
        {'T', -4}, {'W', -4}, {'Y', -4}, {'V', -4}, {'B', -4}, {'J', -4}, {'Z', -4}, {'X', -4},
        {'*', 1}
    }}
};

struct Result {
    int score;
    std::string seq1;
    std::string seq2;
};

Result global_alignment(const std::string& seq1, const std::string& seq2, int gap_penalty = 5) {
    int m = seq1.length();
    int n = seq2.length();
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<char>> traceback(m + 1, std::vector<char>(n + 1, 0));

    for (int i = 0; i <= m; i++) {
        dp[i][0] = -i * gap_penalty;
    }
    for (int j = 0; j <= n; j++) {
        dp[0][j] = -j * gap_penalty;
    }
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int match_score = blosum62[seq1[i-1]][seq2[j-1]];
            int diagonal = dp[i-1][j-1] + match_score;
            int up = dp[i-1][j] - gap_penalty;
            int left = dp[i][j-1] - gap_penalty;
            
            dp[i][j] = std::max({diagonal, up, left});
            
            if (dp[i][j] == left) {
                traceback[i][j] = 'L';
            } else if (dp[i][j] == up) {
                traceback[i][j] = 'U';
            } else {
                traceback[i][j] = 'D';
            }
        }
    }
    std::string align1, align2;
    int i = m, j = n;
    
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && traceback[i][j] == 'D') {
            align1 = seq1[i-1] + align1;
            align2 = seq2[j-1] + align2;
            i--; j--;
        } else if (i > 0 && (j == 0 || traceback[i][j] == 'U')) {
            align1 = seq1[i-1] + align1;
            align2 = '-' + align2;
            i--;
        } else {
            align1 = '-' + align1;
            align2 = seq2[j-1] + align2;
            j--;
        }
    }
    
    return {dp[m][n], align1, align2};
}

std::vector<std::string> read_input(const std::string& filename) {
    std::vector<std::string> sequences;
    std::string current_seq;
    std::ifstream file(filename);
    std::string line;
    
    while (std::getline(file, line)) {
        if (line[0] == '>') {
            if (!current_seq.empty()) {
                sequences.push_back(current_seq);
                current_seq.clear();
            }
        } else {
            current_seq += line;
        }
    }
    if (!current_seq.empty()) {
        sequences.push_back(current_seq);
    }
    
    return sequences;
}

void output_to_file(string input,string filename){
    std::ofstream out(filename);
    out << input;
    out.close();
}

#include <sstream>

vector<string> split(string s,char c){

    std::stringstream test(s);
    std::string segment;
    std::vector<std::string> seglist;
    while(std::getline(test, segment, c))
    {
       seglist.push_back(segment);
    }
    return seglist;
}

int main(int argc, char** argv) {
    

    if(argc==2 ){
        string cmd = argv[1];
        string file = string(argv[1]);

        string tfname = "1_global_alignment_with_a_fixed_indel_penalty/1_global_alignment_with_a_fixed_indel_penalty/1_global_alignment_with_a_fixed_indel_penalty/Debugging/inputs/input_1.txt";

        auto stripped = split(file,'.')[0];
        auto id = split(stripped,'_').back();
        int fnum = atoi(id.c_str());
    
        // Read sequences from file
        auto sequences = read_input(file);
        if (sequences.size() < 2) {
            cout << "Error: Need at least 2 sequences" << endl;
            return 1;
        }
        auto result = global_alignment(sequences[0], sequences[1]);
        
        // Print results
        cout << result.score << endl;
        cout << result.seq1 << endl;
        cout << result.seq2 << endl;

        string output = to_string(result.score)+"\n";
        output += result.seq1+"\n";
        output += result.seq2;

        output_to_file(output,"sol_q1_t"+to_string(fnum)+".txt");
    }  
    
    return 0;
}