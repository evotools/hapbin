/*
 * Hapbin: A fast binary implementation EHH, iHS, and XPEHH
 * Copyright (C) 2014-2017 Colin MacLean <cmaclean@illinois.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef ARGPARSE_HPP
#define ARGPARSE_HPP

#include <cstdlib>
#include <string>
#include <tuple>
#include <vector>
#include <iostream>
#include <type_traits>
#include <cstring>
#include <cassert>
#include <utility>
#include <functional>
#include <cstdlib>
#include <climits>
#include "config.h"

bool arithParseErrorCheck(const char* optLong, const char* last)
{
    int err = errno;
    if (err == EINVAL || *last != '\0')
    {
        std::cerr<< "Could not parse value for --" << optLong << std::endl;
        return false;
    }
    if (err == ERANGE)
    {
        std::cerr<< "Value out of range of type for --" << optLong << std::endl;
        return false;
    }
    return true;
}

void removeBackslashes(char* in)
{
    std::size_t len = strlen(in);
    int i = 0, pos=0;
    while(i < len)
    {
        if (in[i] == '\\')
        {
            in[pos] = in[i+1];
            pos++;
            i += 2;
        }
        else
        {
            in[pos] = in[i];
            pos++;
            i++;
        }
    }
    in[pos] = '\0';
}

template<typename T>
T parseArg(char *arg, const char* optLong, bool *ok);

template<>
int parseArg(char *arg, const char* optLong, bool *ok)
{
    char *tmp;
    int ret = strtol(arg, &tmp, 0);
    *ok = arithParseErrorCheck(optLong, tmp);
    return ret;
}

template<>
unsigned int parseArg(char* arg, const char* optLong, bool *ok)
{
    char *tmp;
    unsigned int ret = strtoul(arg, &tmp, 0);
    *ok = arithParseErrorCheck(optLong, tmp);
    return ret;
}

template<>
long int parseArg(char *arg, const char* optLong, bool *ok)
{
    char *tmp;
    long ret = strtol(arg, &tmp, 0);
    *ok = arithParseErrorCheck(optLong, tmp);
    return ret;
}

template<>
long long int parseArg(char *arg, const char* optLong, bool *ok)
{
    char *tmp;
    long long ret = strtoll(arg, &tmp, 0);
    *ok = arithParseErrorCheck(optLong, tmp);
    return ret;
}

template<>
unsigned long int parseArg(char* arg, const char* optLong, bool *ok)
{
    char *tmp;
    unsigned long ret = strtoul(arg, &tmp, 0);
    *ok = arithParseErrorCheck(optLong, tmp);
    return ret;
}

template<>
unsigned long long int parseArg(char* arg, const char* optLong, bool *ok)
{
    char *tmp;
    unsigned long long ret = strtoull(arg, &tmp, 0);
    *ok = arithParseErrorCheck(optLong, tmp);
    return ret;
}

template<>
const char* parseArg(char* arg, const char* optLong, bool *ok)
{
    *ok = true;
    if (arg[0] == '-')
    {
        std::cerr << "Error: --" << optLong << " expects a value." << std::endl;
        *ok = false;
    }
    removeBackslashes(arg);
    return arg;
}

template<>
float parseArg(char* arg, const char* optLong, bool *ok)
{
    char *tmp;
    float ret = strtof(arg, &tmp);
    *ok = arithParseErrorCheck(optLong, tmp);
    return ret;
}

template<>
double parseArg(char* arg, const char* optLong, bool *ok)
{
    char *tmp;
    double ret;
    try {
        ret = std::stod(arg);
        *ok=true;
    } catch (const std::invalid_argument& ia) {
        std::cout << "Invalid argument. Could not parse: " << arg << std::endl;
	    *ok=false;
    }
    return ret;
}

template<>
std::string parseArg(char* arg, const char* optLong, bool *ok)
{
    *ok = true;
    if (arg[0] == '-')
    {
        std::cerr << "Error: --" << optLong << " expects a value." << std::endl;
        *ok = false;
    }
    removeBackslashes(arg);
    return std::string(arg);
}

class ArgumentBase
{
public:
    ArgumentBase(char opt, const char* optLong, const char* desc, bool multi, bool required)
        : m_opt(opt), m_optLong(optLong), m_desc(desc), m_found(false), m_multi(multi), m_required(required)
        {}
    const char* desc() const { return m_desc; }
    char opt() const { return m_opt; }
    const char* optLong() const { return m_optLong; }
    void found() { m_found = true; }
    bool wasFound() const { return m_found; }
    bool isMulti() const { return m_multi; }
    bool isRequired() const { return m_required; }
    virtual int parse(int argc, char** argv) = 0;

    static const char NO_SHORT_OPT;
protected:
    const char* m_optLong;
    char m_opt;
    const char* m_desc;
    bool m_found;
    bool m_multi;
    bool m_required;
};

const char ArgumentBase::NO_SHORT_OPT = '\0';

template<typename T>
class Argument : public ArgumentBase
{
public:
    Argument(char opt, const char* optLong, const char* desc, bool multi, bool required, T defaultValue)
        : ArgumentBase(opt, optLong, desc, multi, required), m_vals{}, m_val{defaultValue} {}
    virtual int parse(int argc, char** argv)
    {
        int arglength = strlen(argv[0]);
        int optlength = strlen(optLong());
        if (arglength == 2)
        {
            if(argv[0][0] == '-' && argv[0][1] == m_opt)
            {
                return doParse(argv[1]);
            }
        }
        if (arglength == optlength + 2)
        {
            if (strncmp("--", argv[0], 2) == 0 && strcmp(optLong(), &argv[0][2]) == 0)
            {
                return doParse(argv[1]);
            }
        }
        return 0;
    }
    T value() { return m_val; }
    std::vector<T> values() const { return m_vals; }
    bool isFlag() const { return false; }
    bool toggleFlag() {}
protected:
    int doParse(char* v)
    {
        bool ok;
        if (isMulti())
        {
            if(wasFound())
            {

                T val = parseArg<T>(v, optLong(), &ok);
                if (!ok)
                    return INT_MAX;
                m_vals.push_back(val);
                m_val = val;
            }
            else
            {
                T val = parseArg<T>(v, optLong(), &ok);
                if (!ok)
                    return INT_MAX;
                m_vals.push_back(val);
                m_val = val;
                found();
            }
            return 2;
        } else {
            T val = parseArg<T>(v, optLong(), &ok);
            if (!ok)
                return INT_MAX;
            m_val = val;
            found();
            return 2;
        }
    }
    std::vector<T> m_vals;
    T m_val;
};

template<>
class Argument<bool> : public ArgumentBase
{
public:
    Argument(char opt, const char* optLong, const char* desc, bool flag, bool defaultValue)
        : ArgumentBase(opt, optLong, desc, false, false), m_isFlag(flag), m_val(defaultValue)
        {
            //flags should default to false
            assert(!(flag && defaultValue));
        }
    virtual int parse(int argc, char** argv)
    {
        int arglength = strlen(argv[0]);
        int optlength = strlen(optLong());
        if(arglength == 2 && argv[0][0] == '-' && argv[0][1] == m_opt)
        {
            m_val = true;
            found();
            return 1;
        }
        if (isFlag())
        {
            if (arglength == optlength + 2)
            {
                if (strncmp("--", argv[0], 2) == 0 && strcmp(optLong(), &argv[0][2]) == 0)
                {
                    m_val = true;
                    found();
                    return 1;
                }
            }
        }
        else
        {
            if (arglength == optlength + 5)
            {
                if (strncmp("--no-", argv[0], 4) == 0 && strcmp(optLong(), &argv[0][5]) == 0)
                {
                    m_val = false;
                    found();
                    return 1;
                }
            }
            else if (arglength == optlength + 2)
            {
                if (strncmp("--", argv[0], 2) == 0 && strcmp(optLong(), &argv[0][2]) == 0)
                {
                    m_val = true;
                    found();
                    return 1;
                }
            }
        }
        return 0;
    }
    bool value() const { return m_val; }
    bool isFlag() const { return m_isFlag; }
    bool toggleFlag(bool val) { m_val = val; }
protected:
    bool m_val;
    bool m_defaultValue;
    bool m_isFlag;
};

class BareArgument : public ArgumentBase
{
public:
    BareArgument(bool multi, bool required) : ArgumentBase('\0', NULL, "", multi, required) {}

    virtual int parse(int argc, const char** argv)
    {

    }
protected:
    std::vector<const char*> m_vals;
};

class ArgParse
{
public:
    ArgParse(const std::vector<ArgumentBase*>& args, const std::string& usageSummary = "")
        : m_current_argument(1), m_arguments(args), m_usageSummary(usageSummary)
    {}

    bool parseArguments(int argc, char** argv)
    {
        for(int currArg = 1; currArg < argc;)
        {
            bool found = false;
            for (ArgumentBase* b : m_arguments)
            {
                int advance = b->parse(argc-currArg, &argv[currArg]);
                if (advance == INT_MAX)
                    return false;
                if (advance > 0)
                {
                    currArg += advance;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                std::cerr << "Error: Unknown argument " << argv[currArg] << std::endl;
                return false;
            }
        }
        bool missingArg = false;
        for (ArgumentBase *b : m_arguments)
        {
            if (b->isRequired() && !b->wasFound())
            {
                std::cerr << "Error: Must specify --" << b->optLong() << " argument." << std::endl;
                missingArg = true;
            }
        }
        if (missingArg)
        {
            showHelp();
            return false;
        }
        return true;
    }
    void showHelp() const
    {
        std::cout << m_usageSummary << std::endl << std::endl;
        for (ArgumentBase* b : m_arguments)
        {
            std::cout << "\t";
            if (b->opt() != ArgumentBase::NO_SHORT_OPT)
                std::cout << "-" << b->opt() << ",";
            std::cout << "--" << b->optLong() << "\t\t\t" << b->desc() << std::endl;
        }
    }
    void showVersion() const
    {
        std::cout << "Version: " << VERSION_SHORT << " (" << VERSION << ") " << std::endl;
    }
protected:
    int m_current_argument;
    std::string m_usageSummary;
    std::vector<ArgumentBase*> m_arguments;
};

template<typename... Args>
void parseArguments(Argument<Args>&... args)
{
}

#endif /* ARGPARSE_HPP */
