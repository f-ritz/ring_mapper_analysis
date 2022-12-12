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

    inline ThreeDInfoVector(const int &bracketNumber, const int &basepairNumber, const int &pairNumber)
            : _bracketNumber(
            bracketNumber), _basepairNumber(basepairNumber),
              _pairNumber(
                      pairNumber) {}

    inline explicit ThreeDInfoVector(const std::vector<int> &v) : _bracketNumber(v[0]), _basepairNumber(v[1]),
                                                                  _pairNumber(v[2]) {
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
        }
        return *this;
    }

    /// @brief - vector + vector
    friend inline ThreeDInfoVector operator+(const ThreeDInfoVector &a, const ThreeDInfoVector &b) {
        return {a._bracketNumber + b._bracketNumber, a._basepairNumber + b._basepairNumber,
                a._pairNumber + b._pairNumber};
    }

    /// @brief vector + int
    friend inline ThreeDInfoVector operator+(const ThreeDInfoVector &v, const int &t) {
        return {v._bracketNumber + t, v._basepairNumber + t, v._pairNumber + t};
    }

    /// @brief int + vector
    friend inline ThreeDInfoVector operator+(const int &t, const ThreeDInfoVector &v) {
        return {t + v._bracketNumber, t + v._basepairNumber, t + v._pairNumber};
    }

    /// @brief - vector - vector
    friend inline ThreeDInfoVector operator-(const ThreeDInfoVector &a, const ThreeDInfoVector &b) {
        return {a._bracketNumber - b._bracketNumber, a._basepairNumber - b._basepairNumber,
                a._pairNumber - b._pairNumber};
    }

    /// @brief vector - int
    friend inline ThreeDInfoVector operator-(const ThreeDInfoVector &v, const int &t) {
        return {v._bracketNumber - t, v._basepairNumber - t, v._pairNumber - t};
    }

    /// @brief int - vector
    friend inline ThreeDInfoVector operator-(const int &t, const ThreeDInfoVector &v) {
        return {t - v._bracketNumber, t - v._basepairNumber, t - v._pairNumber};
    }

    /// @brief += vector
    inline ThreeDInfoVector &operator+=(const ThreeDInfoVector &v) {
        _bracketNumber += v._bracketNumber;
        _basepairNumber += v._basepairNumber;
        _pairNumber += v._pairNumber;
        return *this;
    }

    /// @brief += double
    inline ThreeDInfoVector &operator+=(const int &t) {
        _bracketNumber += t;
        _basepairNumber += t;
        _pairNumber += t;
        return *this;
    }

    /// @brief - getters
    inline int get_bracket_number() const { return _bracketNumber; }

    inline int get_basepair_number() const { return _basepairNumber; }

    inline int get_pair_number() const { return _pairNumber; }

    /// @brief - setters
    inline void set_bracket_number(const int &bracket_number) { _bracketNumber = bracket_number; }

    inline void set_basepair_number(const int &basepair_number) { _basepairNumber = basepair_number; }

    inline void set_pair_number(const int &pair_number) { _pairNumber = pair_number; }

    /// @brief - element addition
    inline void add_bracket_number() { _bracketNumber++; }

    inline void add_basepair_number() { _basepairNumber++; }

    inline void add_pair_number() { _pairNumber++; }

    /// @brief - element subtraction
    inline void subtract_bracket_number() { _bracketNumber--; }

    inline void subtract_basepair_number() { _basepairNumber--; }

    inline void subtract_pair_number() { _pairNumber--; }

    /// @brief - storers
    // inline void store_bracket_number() {}
    // maybe unused

private:
    /// @brief - values of the 3D info vector
    int _bracketNumber;
    int _basepairNumber;
    int _pairNumber;
};


#endif //RING_MAPPER_ANALYSIS_THREEDINFOVECTOR_H
