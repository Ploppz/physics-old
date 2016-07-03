#include "StatisticsCollection.h"
#include "Statistics.h"
#include "debug.h"

#include <string>

void StatisticsCollection::add_value(const std::string& name, float value)
{
    auto it = statistics_objects.find(name);
    if (it == statistics_objects.end()) {
        // add new Statistics object
        statistics_objects[name] = Statistics(value);
    } else {
        // use existing Statistics object
        it->second.add_value(value);
    }
}
void StatisticsCollection::count(const std::string& name)
{
    DebugBegin();
    auto it = counters.find(name);
    if (it == counters.end()) {
        // add new Statistics object
        dout << "Init counter " << name << newl;
        counters[name] = 0;
    } else {
        // use existing Statistics object
        dout << "Add to counter " << name << newl;
        ++ it->second;
    }
}
std::map<std::string, Statistics>::iterator StatisticsCollection::begin()
{
    return statistics_objects.begin();
}
std::map<std::string, Statistics>::iterator StatisticsCollection::end()
{
    return statistics_objects.end();
}
std::map<std::string, int>::iterator StatisticsCollection::begin_counters()
{
	return counters.begin();
}

std::map<std::string, int>::iterator StatisticsCollection::end_counters()
{
	return counters.end();
}

