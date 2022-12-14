#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <numeric>
#include <stack>

#include "../base/functions.h"
#include "../base/ThreeDInfoVector.h"

int main(int argc, char **argv) {

    std::cout << "Processing..." << std::endl;
    std::cout << " " << std::endl;

    if (argc > 6) {
        std::cout << "Too many command-line arguments!!" << std::endl;
        std::cout
                << "./ring_mapper_analysis <directory> <structure directory> <number of reads> <i to retrieve> <j to retrieve>"
                << std::endl;
        exit(0);
    } else if (argc < 1) {
        std::cout << "Need command-line arguments for file!!" << std::endl;
        std::cout
                << "./ring_mapper_analysis <directory> <structure directory> <number of reads> <i to retrieve> <j to retrieve>"
                << std::endl;
        exit(0);
    }

    int number_of_reads;
    number_of_reads = std::stoi(argv[3]);
    int i_retrieval;
    i_retrieval = std::stoi(argv[4]);
    int j_retrieval;
    j_retrieval = std::stoi(argv[5]);
    int n = 0;
    int counter = 0;
    int vector_size;
    int num_mutations_int;
    int sum_muts;

    String filename(argv[1]);
    String structure_filename(argv[2]);
    String line;
    String data;
    String num_mutations_str;
    Strings data_temp;
    Strings lines;
    std::ifstream input_file(filename);

    if (!input_file.is_open()) {
        std::cerr << "Could not open the file - '" << filename << "'" << std::endl;
        return EXIT_FAILURE;
    }

    String data_counter_single;
    String data_counter_line;
    Strings data_counter;
    Strings mutation_counter;

    for (int i = 0; i < 4; i++) {
        getline(input_file, data_counter_line);
        data_counter = split(data_counter_line, "\t");
        if (i >= 3) {
            data_counter_single = data_counter[1];
        }
    }

    vector_size = data_counter_single.length();
    std::ofstream ring_mapper_analysis_csv;
    ring_mapper_analysis_csv.open("ring_mapper_results.csv");
    ring_mapper_analysis_csv << "i" << "," << "j" << "," << "chisq" << "," << "a" << "," << "b" << "," << "c" << ","
                             << "d" << std::endl;
    std::ofstream ring_mapper_histogram_csv;
    ring_mapper_histogram_csv.open("ring_mapper_histogram.csv");
    std::ofstream ring_mapper_retrieval_csv;
    ring_mapper_retrieval_csv.open("ring_mapper_temp.csv");
    std::ofstream ring_mapper_histogram_positions_csv;
    ring_mapper_histogram_positions_csv.open("ring_mapper_counter.csv");
    std::ofstream ring_mapper_retrieval_real_csv;
    ring_mapper_retrieval_real_csv.open("ring_mapper_retrieval.csv");
    ring_mapper_retrieval_real_csv << "i" << "," << "j" << std::endl;
    std::vector<std::vector<float>> a;
    std::vector<std::vector<float>> b;
    std::vector<std::vector<float>> c;
    std::vector<std::vector<float>> d;

    for (int i = 0; i < vector_size; i++) {
        std::vector<float> v;
        for (int j = 0; j < vector_size; j++) {
            v.push_back(0.0);
        }
        a.push_back(v);
        b.push_back(v);
        c.push_back(v);
        d.push_back(v);
    }

    double five_percent_cutoff = 0.05 * vector_size;
    std::vector<float> mut_positions;
    std::vector<float> mut_positions_frequency;
    for (int i = 0; i < vector_size; i++) {
        mut_positions.push_back(0.0f);
        mut_positions_frequency.push_back(0);
    }

    /// @brief - ring mapper algorithm
    while (input_file.good()) {
        getline(input_file, line);
        if (line.length() < 2) {
            break;
        }
        data_temp = split(line, "\t");
        if (counter > 2) {
            data = data_temp[1];
            num_mutations_int = std::stoi(data_temp[2]);
            n += 1;
            if (n > number_of_reads) {
                break;
            }
            if (num_mutations_int > 0) {
                if (num_mutations_int < five_percent_cutoff) {
                    for (int i = 0; i < data.length(); i++) {
                        int mut_i = 0;
                        if ((data[i] == 'A') || (data[i] == 'T') || (data[i] == 'C') || (data[i] == 'G')) {
                            mut_i = 1;
                        }
                        for (int j = i + 1; j < data.length(); j++) {
                            int mut_j = 0;
                            if ((data[j] == 'A') || (data[j] == 'T') || (data[j] == 'C') || (data[j] == 'G')) {
                                mut_j = 1;
                            }
                            if (mut_i == 0 && mut_j == 0) {
                                a[i][j]++;
                            } else if (mut_i == 1 && mut_j == 0) {
                                if (num_mutations_int == 1) {
                                    mut_positions[i]++;
                                }
                                b[i][j]++;
                            } else if (mut_i == 0) {
                                c[i][j]++;
                            } else {
                                d[i][j]++;
                            }
                            // creates a temp file
                            if (i == i_retrieval) {
                                if (j == j_retrieval) {
                                    ring_mapper_retrieval_csv << data[i_retrieval] << "," << data[j_retrieval]
                                                              << std::endl;
                                }
                            }

                        }
                    }
                }
            }
        }
        counter++;
    }

    /// @brief - retrieval
    String temp_filename;
    temp_filename = "ring_mapper_temp.csv";
    std::ifstream input_temp_file(temp_filename);
    String temp_line;
    Strings data_temp_vector;
    while (input_temp_file.good()) {
        getline(input_temp_file, temp_line);
        data_temp_vector = split(temp_line, ",");
        String data_i = data_temp_vector[0];
        String data_j = data_temp_vector[1];

        if (data_i == "?" || data_i == ".") {
            data_i = " ";
            data_j = " ";
        }

        if (data_j == "?" || data_j == ".") {
            data_i = " ";
            data_j = " ";
        }

        if (data_i == "A" || data_i == "T" || data_i == "C" || data_i == "G") {
            data_i = "1";
        }

        if (data_j == "A" || data_j == "T" || data_j == "C" || data_j == "G") {
            data_j = "1";
        }

        if (data_i == "0" || data_j == "0" || data_i == "1" || data_j == "1") {
            ring_mapper_retrieval_real_csv << data_i << "," << data_j << std::endl;

        }
    }

    /// @brief - histogram
    sum_muts = std::accumulate(mut_positions.begin(), mut_positions.end(), 0);
    ring_mapper_histogram_csv << "i,histogram_values" << std::endl;
    ring_mapper_histogram_positions_csv << "i, num_single_muts, total_muts: " << sum_muts << std::endl;
    for (int i = 0; i < vector_size; i++) {
        mut_positions_frequency[i] = (mut_positions[i] / sum_muts);
        ring_mapper_histogram_positions_csv << i + 1 << "," << mut_positions[i] << std::endl;
        ring_mapper_histogram_csv << i + 1 << "," << mut_positions_frequency[i] << std::endl;
    }
    std::cout << "Histogram finished..." << std::endl;

    // temp file creation for pairing algorithm
    std::ofstream ring_mapper_positions;
    ring_mapper_positions.open("ring_mapper_positions.csv");

    /// @brief - ring mapper results
    for (int i = 0; i < vector_size; i++) {
        for (int j = i + 1; j < vector_size; j++) {
            float dim = (a[i][j] + b[i][j]) * (c[i][j] + d[i][j]) * (a[i][j] + c[i][j]) * (b[i][j] + d[i][j]);
            n = a[i][j] + b[i][j] + c[i][j] + d[i][j];
            float chisq;
            chisq = n * ((abs(a[i][j] * d[i][j] - b[i][j] * c[i][j]) - 0.5 * n) *
                         (abs(a[i][j] * d[i][j] - b[i][j] * c[i][j]) - 0.5 * n)) / dim;
            if (dim == 0) {
                continue;
            }
            if (chisq > 20) {
                ring_mapper_analysis_csv << i + 1 << "," << j + 1 << "," << chisq << "," << a[i][j] << "," << b[i][j]
                                         << "," << c[i][j] << "," << d[i][j] << std::endl;
                ring_mapper_positions << i + 1 << "," << j + 1 << std::endl;
            }
        }
    }
    std::cout << "Chisq finished..." << std::endl;

    /// @brief - pairing algorithm
    // variables are a clusterfuck
    String analysis_filename = "ring_mapper_results.csv";
    std::ifstream input_analysis_file(analysis_filename);
    std::ifstream input_structure_file(structure_filename);
    String structure_raw;
    String i_pos;
    String j_pos;
    String pos_temp;
    String base_chain;
    String basepair_structure_chain;
    String temp;
    String temp_1;
    String basepair_i_string;
    String basepair_j_string;
    String ring_mapper_positions_line;
    String ring_mapper_input;
    Strings pos_temp_vector;
    Strings i_pos_vector_temp;
    Strings j_pos_vector_temp;
    Strings i_pos_vector_str;
    Strings j_pos_vector_str;
    Strings cooked_structure;
    Strings basepairs;
    Strings basepair_structure;
    Strings ring_mapper_inputs;
    Strings i_positions;
    Strings j_positions;
    Ints i_pos_vector;
    Ints j_pos_vector;

    /// @brief - structure processing
    // initializes a vector that contains the processed (cooked) structure
    for (int i = 0; i < 3; i++) {
        cooked_structure.push_back("");
    }
    // writes each line in the input text file "structure_input.txt" to a spot in cooked_structure
    for (int i = 0; i < 3; i++) {
        getline(input_structure_file, structure_raw);
        cooked_structure[i] = structure_raw;
    }
    // processes the cooked_structure and retrieves proper data
    base_chain = cooked_structure[0];
    basepair_structure_chain = cooked_structure[1];
    // more vectors, this time processing the data into its final form to be ready to plug in the algorithm
    for (int i = 0; i < vector_size; i++) {
        temp = base_chain[i];
        temp_1 = basepair_structure_chain[i];
        basepairs.push_back(temp);
        basepair_structure.push_back(temp_1);
    }

    /// @brief - position processing, retrieves positions of basepairs
    int pos_loop_counter;
    while (ring_mapper_positions.good()) {
        if (pos_loop_counter > 0) {
            getline(input_analysis_file, pos_temp);
            pos_temp_vector = split(pos_temp, ",");
            i_pos = pos_temp_vector[0];
            j_pos = pos_temp_vector[1];
            if (pos_temp.length() < 2) {
                break;
            }
            i_pos_vector_temp.push_back(i_pos);
            j_pos_vector_temp.push_back(j_pos);
            //ring_mapper_pos_csv << i_pos << "," << j_pos << std::endl;
        }
        pos_loop_counter++;
    }

    // temp pos vectors go into str vectors
    for (int i = 0; i < i_pos_vector_temp.size(); i++) {
        i_pos_vector_str.push_back(i_pos_vector_temp[i + 1]);
        j_pos_vector_str.push_back(j_pos_vector_temp[i + 1]);
    }

    /// @brief - build the RNA structure using information vectors
    // stack creation
    std::stack<ThreeDInfoVector> stack_info_vectors;
    ThreeDInfoVector twoDInfoVector = {0, 0, 0};
    ThreeDInfoVector threeDInfoVector = {0, 0, 0}; // vector args: {bracketNumber, basepairNumber, pairNumber}
    ThreeDInfoVector emptyVector = {0, 0, 0};
    ThreeDInfoVector tempUsageVector = {0, 0, 0}; // vector for temporary/transfer usage throughout algo
    std::vector<ThreeDInfoVector> basepair_structure_info_vecs;
    for (int i = 0; i < basepair_structure.size(); i++) {
        basepair_structure_info_vecs.push_back(emptyVector);
    }

    for (int y = 0; y < basepair_structure.size(); y++) {
        // motifs
        if (basepair_structure[y] == ".") {
            if (basepair_structure[y - 1] != ".") {
                twoDInfoVector.add_bracket_number();
                twoDInfoVector.set_basepair_number(0);
            }
            twoDInfoVector.add_basepair_number();
            basepair_structure_info_vecs[y] = twoDInfoVector;
        }
        // helices
        if (basepair_structure[y] == "(") {
            if (basepair_structure[y - 1] != "(") {
                threeDInfoVector.add_bracket_number();
                threeDInfoVector.set_basepair_number(0);
            }
            threeDInfoVector.add_basepair_number();
            threeDInfoVector.set_pair_number(1);
            stack_info_vectors.push(threeDInfoVector);
            basepair_structure_info_vecs[y] = threeDInfoVector;
        }
        if (basepair_structure[y] == ")") {
            tempUsageVector = stack_info_vectors.top();
            tempUsageVector.set_pair_number(2);
            basepair_structure_info_vecs[y] = tempUsageVector;
            stack_info_vectors.pop();
        }
    }

    /// @brief - debugging code, prints pairing info vectors to a csv
    /*
    std::ofstream pairing_debug;
    pairing_debug.open("pairing_debug.csv");

    for (int z = 0; z < basepair_structure_info_vecs.size(); z++) {
        pairing_debug << basepair_structure_info_vecs[z].get_bracket_number() << "," << basepair_structure_info_vecs[z].get_basepair_number()
                      << "," << basepair_structure_info_vecs[z].get_pair_number() << std::endl;
    }
     */

    /// @brief - sort each pair of positions into WC/NC/LR
    String basepair_i;
    String basepair_j;
    String bp_pos_temp;
    Strings bp_pos_temps;
    ThreeDInfoVector vector_i;
    ThreeDInfoVector vector_j;
    std::ofstream ring_mapper_pairmap;
    ring_mapper_pairmap.open("ring_mapper_pairmap.csv");
    std::ifstream input_position_file("ring_mapper_positions.csv");

    while (ring_mapper_positions.good()) {
        getline(input_position_file, bp_pos_temp);
        if (bp_pos_temp.length() < 2) {
            break;
        }
        bp_pos_temps = split(bp_pos_temp, ",");
        basepair_i = bp_pos_temps[0];
        basepair_j = bp_pos_temps[1];

        vector_i = basepair_structure_info_vecs[std::stoi(basepair_i) - 1];
        vector_j = basepair_structure_info_vecs[std::stoi(basepair_j) - 1];

        if (vector_i.get_pair_number() == 0) {
            if (vector_i.get_pair_number() != vector_j.get_pair_number()) {
                ring_mapper_pairmap << basepair_i << "," << basepair_j << "," << "LR" << std::endl;
            } else if (vector_i.get_bracket_number() == vector_j.get_bracket_number()) {
                ring_mapper_pairmap << basepair_i << "," << basepair_j << "," << "NC" << std::endl;
            } else {
                ring_mapper_pairmap << basepair_i << "," << basepair_j << "," << "LR" << std::endl;
            }
        } else {
            if (vector_i.get_bracket_number() == vector_j.get_bracket_number()) {
                if (vector_i.get_basepair_number() == vector_j.get_basepair_number()) {
                    ring_mapper_pairmap << basepair_i << "," << basepair_j << "," << "WC" << std::endl;
                } else {
                    ring_mapper_pairmap << basepair_i << "," << basepair_j << "," << "LOCAL" << std::endl;
                }
            } else {
                ring_mapper_pairmap << basepair_i << "," << basepair_j << "," << "LR" << std::endl;
            }
        }
    }
    std::cout << "Pairing finished..." << std::endl;

    /// @brief - write everything to a single file

    std::ifstream input_ring_mapper_results("ring_mapper_results.csv");
    std::ifstream input_pairmap_file("ring_mapper_pairmap.csv");

    std::ofstream ring_mapper_final_analysis;
    ring_mapper_final_analysis.open("ring_mapper_analysis.csv");

    String ring_mapper_analysis_line;
    Strings ring_mapper_analysis_lines;

    String ring_mapper_pairmap_line;
    Strings ring_mapper_pairmap_lines;

    int line_counter = 0;

    while (ring_mapper_analysis_csv.good()) {
        getline(input_ring_mapper_results, ring_mapper_analysis_line);
        ring_mapper_analysis_lines = split(ring_mapper_analysis_line, ",");

        if (line_counter >= 1) {
            getline(input_pairmap_file, ring_mapper_pairmap_line);
            ring_mapper_pairmap_lines = split(ring_mapper_pairmap_line, ",");
        }

        if (ring_mapper_analysis_line.length() < 3) {
            break;
        }

        if (line_counter < 1) {
            ring_mapper_final_analysis << ring_mapper_analysis_lines[0] << "," << ring_mapper_analysis_lines[1] << ","
                                       << ring_mapper_analysis_lines[2] << "," << ring_mapper_analysis_lines[3] << ","
                                       << ring_mapper_analysis_lines[4] << "," << ring_mapper_analysis_lines[5] << ","
                                       << ring_mapper_analysis_lines[6] << "," << "basepair_type" << std::endl;
        } else {
            ring_mapper_final_analysis << ring_mapper_analysis_lines[0] << "," << ring_mapper_analysis_lines[1] << ","
                                       << ring_mapper_analysis_lines[2] << "," << ring_mapper_analysis_lines[3] << ","
                                       << ring_mapper_analysis_lines[4] << "," << ring_mapper_analysis_lines[5] << ","
                                       << ring_mapper_analysis_lines[6] << "," << ring_mapper_pairmap_lines[2] << std::endl;
        }

        line_counter++;


        /*
        std::cout << ring_mapper_analysis_lines[0] << std::endl;
        std::cout << ring_mapper_analysis_lines[1] << std::endl;
        std::cout << ring_mapper_analysis_lines[2] << std::endl;
        std::cout << ring_mapper_analysis_lines[3] << std::endl;
        std::cout << ring_mapper_analysis_lines[4] << std::endl;
        std::cout << ring_mapper_analysis_lines[5] << std::endl;
        std::cout << ring_mapper_analysis_lines[6] << std::endl;
        */
        //break;
    }

    /*

    //int line_counter;

    Strings analysis_lines_container;
    Strings pairmap_lines_container;

    while (ring_mapper_analysis_csv.good()) {
        getline(input_analysis_file, ring_mapper_analysis_line);
        ring_mapper_analysis_lines = split(ring_mapper_analysis_line, ",");

        getline(input_pairmap_file, ring_mapper_pairmap_line);
        ring_mapper_pairmap_lines = split(ring_mapper_pairmap_line, ",");

        ring_mapper_analysis_lines = split(ring_mapper_analysis_line, ",");
        ring_mapper_pairmap_lines = split(ring_mapper_pairmap_line, ",");

        std::cout << ring_mapper_analysis_lines[0] << std::endl;
        std::cout << ring_mapper_analysis_lines[1] << std::endl;
        std::cout << ring_mapper_analysis_lines[2] << std::endl;
        std::cout << ring_mapper_analysis_lines[3] << std::endl;
        std::cout << ring_mapper_analysis_lines[4] << std::endl;
        std::cout << ring_mapper_analysis_lines[5] << std::endl;
        std::cout << ring_mapper_analysis_lines[6] << std::endl;
        std::cout << ring_mapper_pairmap_lines[6] << std::endl;


        if (line_counter < 1) {
            ring_mapper_final_analysis << ring_mapper_analysis_lines[0] << "," << ring_mapper_analysis_lines[1] << ","
                                       << ring_mapper_analysis_lines[2] << "," << ring_mapper_analysis_lines[3] << ","
                                       << ring_mapper_analysis_lines[4] << "," << ring_mapper_analysis_lines[5] << ","
                                       << ring_mapper_analysis_lines[6] << "," << ring_mapper_analysis_lines[7] << std::endl;
        } else {
            ring_mapper_final_analysis << ring_mapper_analysis_lines[0] << "," << ring_mapper_analysis_lines[1] << ","
                                       << ring_mapper_analysis_lines[2] << "," << ring_mapper_analysis_lines[3] << ","
                                       << ring_mapper_analysis_lines[4] << "," << ring_mapper_analysis_lines[5] << ","
                                       << ring_mapper_analysis_lines[6] << "," << ring_mapper_pairmap_lines[2] << std::endl;
        }

        line_counter++;

    }
    while (ring_mapper_pairmap.good()) {
        getline(input_pairmap_file, ring_mapper_pairmap_line);
        ring_mapper_pairmap_lines = split(ring_mapper_pairmap_line, ",");
    }
*/





    std::cout << "Sorting results..." << std::endl;

    std::cout << " " << std::endl;

    std::cout << "Process finished!" << std::endl;
}
