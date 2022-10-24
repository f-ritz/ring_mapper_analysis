#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include "../base/functions.h"
#include <istream>
#include <ostream>
#include <stdexcept>

int main(int argc, char **argv) {
    std::cout << "Processing..." << std::endl;
    if (argc > 2) {
        std::cout << "Too many command-line arguments!!" << std::endl;
        exit(0);
    } else if (argc < 1) {
        std::cout << "Need command-line argument for file!!" << std::endl;
        exit(0);
    }
    int n = 0;
    int counter = 0;
    int vector_size = 0;
    int num_mutations_int = 0;
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
    ring_mapper_analysis_csv << "i" << "," << "j" << "," << "chisq" << "," << "a" << "," << "b" << "," << "c" << "," << "d" << std::endl;
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
    while (input_file.good()) {
        getline(input_file, line);
        if (line.length() < 2) {
            break;
        }
        data_temp = split(line, "\t");
        if (counter > 2) /* trims the 1st 3 lines */ {
            data = data_temp[1];
            num_mutations_str = data_temp[2];
            num_mutations_int = std::stoi(num_mutations_str);
            n += 1;
            if (n > 1000) {
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
                ring_mapper_analysis_csv << i + 1 << "," << j + 1 << "," << chisq << "," << a[i][j] << "," << b[i][j] << "," << c[i][j] << "," << d[i][j] << std::endl;
            }
        }
    }
    std::cout << "Process finished!" << std::endl;
}
