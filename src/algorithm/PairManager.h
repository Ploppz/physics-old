#pragma once
#include <unordered_map>

struct Pair {
    Pair(int box1, int box2) : first(box1), second(box2) {}
    Pair(const Pair& other) : first(other.first), second(other.second) {}
    /* The std::pair is returned when iterating the std::unordered_map: */
    Pair(std::pair<const Pair, Pair>& pair_pair) : first(pair_pair.first.first), second(pair_pair.first.second) {}
    Pair() : first(-1), second(-1) {}
    int first;
    int second;
    bool operator== (const Pair& other) const { return first == other.first && second == other.second; }
};
namespace std {
    template <>
    struct hash<Pair>
    {
        std::size_t operator()(const Pair& pair) const
        {
            return hash<int>()(pair.first) + hash<int>()(pair.second);
        }
    };
}

class PairManager
{
 public: /* Interface */
    PairManager();
    void register_pair(int box1, int box2);
    void unregister_pair(int box1, int box2);
    void print();
    // Iterating
    std::unordered_map<Pair, Pair>::iterator begin();
    std::unordered_map<Pair, Pair>::iterator end();
 private: /* Structures */
 private: /* Data */
    std::unordered_map<Pair, Pair> pairs; // could be a hash map.. key == value
 private: /* Methods */
    bool exists(Pair pair);
    bool exists(int box1, int box2);
};

