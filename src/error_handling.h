#pragma once

#include <sstream>
#include <string>

#define runtime_fatal(message) _runtime_fatal(message, __FILE__, __FUNCTION__ , __LINE__);
void _runtime_fatal(std::string message, const std::string& file, const std::string& function, int line);
void logic_fatal(std::string message);



class Formatter
{
public:
    Formatter() {}
    ~Formatter() {}

    template <typename Type>
    Formatter & operator << (const Type & value)
    {
        stream_ << value;
        return *this;
    }

    std::string str() const         { return stream_.str(); }
    operator std::string () const   { return stream_.str(); }

    enum ConvertToString 
    {
        to_str
    };
    std::string operator >> (ConvertToString) { return stream_.str(); }

private:
    std::stringstream stream_;

    Formatter(const Formatter &);
    Formatter & operator = (Formatter &);
};
