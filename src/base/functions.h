//
// Created by Joe Yesselman on 8/18/22.
//

#ifndef RING_MAPPER_ANALYSIS_FUNCTIONS_H
#define RING_MAPPER_ANALYSIS_FUNCTIONS_H

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



#endif //RING_MAPPER_ANALYSIS_FUNCTIONS_H
