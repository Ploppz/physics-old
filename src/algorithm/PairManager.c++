#include "PairManager.h"
#include <iostream>


PairManager::PairManager() : pairs {}
{
}

void PairManager::register_pair(int box1, int box2)
{
    Pair new_pair(box1, box2);
    if (!pairs.count(new_pair)) {
        pairs[new_pair] = new_pair;
    }
}
void PairManager::unregister_pair(int box1, int box2)
{
    auto found_it = pairs.find(Pair(box1, box2));
    if (found_it != pairs.end())
        pairs.erase(found_it);
}

bool PairManager::exists(Pair pair)
{
    return pairs.count(pair);
}
bool PairManager::exists(int box1, int box2)
{
    return pairs.count(Pair(box1, box2));
}

void PairManager::print()
{
    for (auto it = pairs.begin(); it != pairs.end(); ++ it)
    {
        std::cout << "Pair: " << it->first.first << ", " << it->first.second << std::endl;
    }
}
std::unordered_map<Pair, Pair>::iterator PairManager::begin()
{
    return pairs.begin();
}
std::unordered_map<Pair, Pair>::iterator PairManager::end()
{
    return pairs.end();
}
