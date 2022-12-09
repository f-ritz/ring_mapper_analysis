//
// Created by Joe Yesselman on 12/7/22.
//

#ifndef RING_MAPPER_ANALYSIS_THREEDINFOVECTOR_H
#define RING_MAPPER_ANALYSIS_THREEDINFOVECTOR_H

#include <vector>

class ThreeDInfoVector {
public:
    inline ThreeDInfoVector() = default;

    inline ThreeDInfoVector(const ThreeDInfoVector &v) = default;

    inline ThreeDInfoVector(const int &chainNumber, const int &chainPosition, const int &nullPos1, const int &nullPos2,
                            const int &nullPos3) : _bracketNumber(chainNumber), _basepairNumber(chainPosition),
                                                   _pairNumber(nullPos1), _depth(nullPos2),
                                                   _hasBeenSearched(nullPos3) {}

    inline explicit ThreeDInfoVector(const std::vector<int> &v) : _bracketNumber(v[0]), _basepairNumber(v[1]),
                                                                  _pairNumber(v[2]), _depth(v[3]),
                                                                  _hasBeenSearched(v[4]) {
    }

public: /// @brief - deletion
    inline ~ThreeDInfoVector() = default;

public: /// @brief - operators
    /// @brief - copy assignment
    inline ThreeDInfoVector &operator=(const ThreeDInfoVector &v) {
        if (this != &v) {
            _bracketNumber = v._bracketNumber;
            _basepairNumber = v._basepairNumber;
            _pairNumber = v._pairNumber;
            _depth = v._depth;
            _hasBeenSearched = v._hasBeenSearched;
        }
        return *this;
    }

    /// @brief - vector + vector
    friend inline ThreeDInfoVector operator+(const ThreeDInfoVector &a, const ThreeDInfoVector &b) {
        return {a._bracketNumber + b._bracketNumber, a._basepairNumber + b._basepairNumber,
                a._pairNumber + b._pairNumber,
                a._depth + b._depth, a._hasBeenSearched + b._hasBeenSearched};
    }

    /// @brief vector + int
    friend inline ThreeDInfoVector operator+(const ThreeDInfoVector &v, const int &t) {
        return {v._bracketNumber + t, v._basepairNumber + t, v._pairNumber + t, v._depth + t, v._hasBeenSearched + t};
    }

    /// @brief int + vector
    friend inline ThreeDInfoVector operator+(const int &t, const ThreeDInfoVector &v) {
        return {t + v._bracketNumber, t + v._basepairNumber, t + v._pairNumber, t + v._depth, t + v._hasBeenSearched};
    }

    /// @brief - vector - vector
    friend inline ThreeDInfoVector operator-(const ThreeDInfoVector &a, const ThreeDInfoVector &b) {
        return {a._bracketNumber - b._bracketNumber, a._basepairNumber - b._basepairNumber,
                a._pairNumber - b._pairNumber,
                a._depth - b._depth, a._hasBeenSearched - b._hasBeenSearched};
    }

    /// @brief vector - int
    friend inline ThreeDInfoVector operator-(const ThreeDInfoVector &v, const int &t) {
        return {v._bracketNumber - t, v._basepairNumber - t, v._pairNumber - t, v._depth - t, v._hasBeenSearched - t};
    }

    /// @brief int - vector
    friend inline ThreeDInfoVector operator-(const int &t, const ThreeDInfoVector &v) {
        return {t - v._bracketNumber, t - v._basepairNumber, t - v._pairNumber, t - v._depth, t - v._hasBeenSearched};
    }

    /// @brief += vector
    inline ThreeDInfoVector &operator+=(const ThreeDInfoVector &v) {
        _bracketNumber += v._bracketNumber;
        _basepairNumber += v._basepairNumber;
        _pairNumber += v._pairNumber;
        _depth += v._depth;
        _hasBeenSearched += v._hasBeenSearched;
        return *this;
    }

    /// @brief += double
    inline ThreeDInfoVector &operator+=(const int &t) {
        _bracketNumber += t;
        _basepairNumber += t;
        _pairNumber += t;
        _depth += t;
        _hasBeenSearched += t;
        return *this;
    }

    /// @brief - getters
    inline int get_bracket_number() const { return _bracketNumber; }

    inline int get_basepair_number() const { return _basepairNumber; }

    inline int get_pair_number() const { return _pairNumber; }

    inline int get_helix_depth() const { return _depth; }

    inline int get_has_been_searched() const { return _hasBeenSearched; }

    /// @brief - setters
    inline void set_bracket_number(const int &bracket_number) { _bracketNumber = bracket_number; }

    inline void set_basepair_number(const int &basepair_number) { _basepairNumber = basepair_number; }

    inline void set_pair_number(const int &pair_number) { _pairNumber = pair_number; }

    inline void set_depth(const int &depth) { _depth = depth; }

    inline void set_has_been_searched(const double &has_been_searched) { _hasBeenSearched = has_been_searched; }

    /// @brief - element addition
    inline void add_bracket_number() { _bracketNumber++; }

    inline void add_basepair_number() { _basepairNumber++; }

    inline void add_pair_number() { _pairNumber++; }

    inline void add_depth() { _depth++; }

    inline void add_has_been_searched() { _hasBeenSearched++; }

    /// @brief - element subtraction
    inline void subtract_bracket_number() { _bracketNumber--; }

    inline void subtract_basepair_number() { _basepairNumber--; }

    inline void subtract_pair_number() { _pairNumber--; }

    inline void subtract_depth() { _depth--; }

    inline void subtract_has_been_searched() { _hasBeenSearched--; }

    /// @brief - storers
    // inline void store_bracket_number() {}
    // maybe unused

private:
    /// @brief - values of the 3D info vector
    int _bracketNumber;
    int _basepairNumber;
    int _pairNumber;
    int _depth;
    int _hasBeenSearched;
};


#endif //RING_MAPPER_ANALYSIS_THREEDINFOVECTOR_H
