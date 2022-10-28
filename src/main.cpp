#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include "../base/functions.h"
#include <istream>
#include <ostream>
#include <stdexcept>
#include <numeric>

int main(int argc, char **argv) {
    std::cout << "Processing..." << std::endl;
    if (argc > 4) {
        std::cout << "Too many command-line arguments!!" << std::endl;
        std::cout << "./ring_mapper_analysis <directory> <number of reads> <i to retrieve> <j to involve>" << std::endl;
        exit(0);
    } else if (argc < 1) {
        std::cout << "Need command-line arguments for file!!" << std::endl;
        std::cout << "./ring_mapper_analysis <directory> <number of reads> <i to retrieve> <j to involve>" << std::endl;
        exit(0);
    } /*else if (argc != 4) {
        std::cout << "Argument error!" << std::endl;
        std::cout << "Arguments enter as follows:" << std::endl;
        std::cout << "./ring_mapper_analysis <directory> <number of reads> <i> <j>" << std::endl;
        exit(0);
    }*/
    int number_of_reads;
    number_of_reads = std::stoi(argv[2]);
    int n = 0;
    int counter = 0;
    int vector_size;
    int num_mutations_int;
    int sum_muts;
    String filename(argv[1]);
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
    //std::ofstream ring_mapper_histogram_positions_csv;
    //ring_mapper_histogram_positions_csv.open("ring_mapper_positions.csv");
    std::ofstream debug_output_csv;
    debug_output_csv.open("debug_output.csv");
    /*
    std::ofstream ring_mapper_a_csv;
    ring_mapper_a_csv.open("ring_mapper_a.csv");
    ring_mapper_a_csv << "i,j" << std::endl;
    std::ofstream ring_mapper_b_csv;
    ring_mapper_b_csv.open("ring_mapper_b.csv");
    ring_mapper_b_csv << "i,j" << std::endl;
    std::ofstream ring_mapper_c_csv;
    ring_mapper_c_csv.open("ring_mapper_c.csv");
    ring_mapper_c_csv << "i,j" << std::endl;
    std::ofstream ring_mapper_d_csv;
    ring_mapper_d_csv.open("ring_mapper_d.csv");
    ring_mapper_d_csv << "i,j" << std::endl;
     */
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

    while (input_file.good()) {
        getline(input_file, line);
        if (line.length() < 2) {
            break;
        }
        data_temp = split(line, "\t");
        if (counter > 2) /* trims the 1st 3 lines */ {
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
                        }
                    }
                }
            }
        }
        counter++;
    }

    sum_muts = std::accumulate(mut_positions.begin(), mut_positions.end(), 0);
    ring_mapper_histogram_csv << "i,histogram_values" << std::endl;
    //ring_mapper_histogram_positions_csv << "i, num_single_muts, total_muts: " << sum_muts << std::endl;
    for (int i = 0; i < vector_size; i++) {
        mut_positions_frequency[i] = (mut_positions[i] / sum_muts);
        //ring_mapper_histogram_positions_csv << i + 1 << "," << mut_positions[i] << std::endl;
        ring_mapper_histogram_csv << i + 1 << "," << mut_positions_frequency[i] << std::endl;
    }
    //std::cout << std::accumulate(mut_positions_frequency.begin(), mut_positions_frequency.end(), 0) << std::endl;

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
            }
        }
    }
    std::cout << "Process finished!" << std::endl;
}
