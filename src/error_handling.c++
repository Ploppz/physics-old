#include "error_handling.h"
#include <iostream>

using namespace std;

void _runtime_fatal(string message, const string& file, const string& function, int line)
{
    cout << function << " in " << file << ":" << line << "  -  runtime error: " << message << endl;
    exit(1);
}
void logic_fatal(string message)
{
    cout << __FILE__ << ":" << __LINE__ << "  -  logic error: " << message << endl;
    exit(1);
}
