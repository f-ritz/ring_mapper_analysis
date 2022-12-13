//
// Created by Joe Yesselman on 8/18/22.
//

#ifndef RING_MAPPER_ANALYSIS_FUNCTIONS_H
#define RING_MAPPER_ANALYSIS_FUNCTIONS_H

#include <iostream>

using String = std::string;
using Strings = std::vector<std::string>;
using Char = char;
using Chars = std::vector<char>;
using Ints = std::vector<int>;

inline Strings split(const String &org_s, const String &delimiter) {
    String s = org_s;
    String token;
    Strings tokens;
    size_t pos;
    while ((pos = s.find(delimiter)) != String::npos) {
        token = s.substr(0, pos);
        tokens.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    if (s.length() > 0) {
        tokens.push_back(s);
    }
    if (tokens.empty()) {
        tokens.emplace_back();
    }
    return tokens;
}

inline void num_to_int(String s)
{
    // Initialize a variable
    int num = 0;
    int n = s.length();
    // Iterate till length of the string
    for (int i = 0; i < n; i++)
        // Subtract 48 from the current digit
        num = num * 10 + (int(s[i]) - 48);
}




#endif //RING_MAPPER_ANALYSIS_FUNCTIONS_H
