#include "error_handling.h"
#include <iostream>


void runtime_fatal(std::string message)
{
    std::cout << __FILE__ << ":" << __LINE__ << "  -  runtime error: " << message << std::endl;
    exit(1);
}
void logic_fatal(std::string message)
{
    std::cout << __FILE__ << ":" << __LINE__ << "  -  logic error: " << message << std::endl;
    exit(1);
}
