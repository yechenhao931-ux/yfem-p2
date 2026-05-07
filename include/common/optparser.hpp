// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef OPTPARSER
#define OPTPARSER

#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
/** Class for parsing command-line options.

    The class is initialized with argc and argv, and new options are added with
    the AddOption method. Currently options of type bool, int, double, char*,
    mfem::Array<int>, and mfem::Vector are supported.

    See the MFEM examples for sample use.
*/
class OptionsParser
{
public:
    enum OptionType
    {
        INT,
        DOUBLE,
        CHAR,
        STRING,
        STD_STRING,
        ENABLE,
        DISABLE,
        ARRAY,
        VECTOR
    };

private:
    struct Option
    {
        OptionType type;
        void *var_ptr;
        const char *short_name;
        const char *long_name;
        const char *description;
        bool required;

        Option() = default;

        Option(OptionType type_, void *var_ptr_, const char *short_name_,
               const char *long_name_, const char *description_, bool req)
            : type(type_), var_ptr(var_ptr_), short_name(short_name_),
              long_name(long_name_), description(description_), required(req) {}
    };

    int argc;
    char **argv;
    std::vector<Option> options;
    std::vector<int> option_check;
    // error_type can be:
    //  0 - no error
    //  1 - print help message
    //  2 - unrecognized option at argv[error_idx]
    //  3 - missing argument for the last option argv[argc-1]
    //  4 - option with index error_idx is specified multiple times
    //  5 - invalid argument in argv[error_idx] for option in argv[error_idx-1]
    //  6 - required option with index error_idx is missing
    int error_type, error_idx;

    static void WriteValue(const Option &opt, std::ostream &out);

public:
    /// Construct a command line option parser with 'argc_' and 'argv_'.
    OptionsParser(int argc_, char *argv_[])
        : argc(argc_), argv(argv_)
    {
        error_type = error_idx = 0;
    }

    /** @brief Add a boolean option and set 'var' to receive the value.
        Enable/disable tags are used to set the bool to true/false
        respectively. */
    void AddOption(bool *var, const char *enable_short_name,
                   const char *enable_long_name, const char *disable_short_name,
                   const char *disable_long_name, const char *description,
                   bool required = false)
    {
        options.push_back(Option(ENABLE, var, enable_short_name, enable_long_name,
                                 description, required));
        options.push_back(Option(DISABLE, var, disable_short_name, disable_long_name,
                                 description, required));
    }

    /// Add an integer option and set 'var' to receive the value.
    void AddOption(int *var, const char *short_name, const char *long_name,
                   const char *description, bool required = false)
    {
        options.push_back(Option(INT, var, short_name, long_name, description,
                                 required));
    }

    /// Add a double option and set 'var' to receive the value.
    void AddOption(double *var, const char *short_name, const char *long_name,
                   const char *description, bool required = false)
    {
        options.push_back(Option(DOUBLE, var, short_name, long_name, description,
                                 required));
    }

    /// Add an integer option and set 'var' to receive the value.
    void AddOption(char *var, const char *short_name, const char *long_name,
                   const char *description, bool required = false)
    {
        options.push_back(Option(CHAR, var, short_name, long_name, description,
                                 required));
    }

    /// Add a string (char*) option and set 'var' to receive the value.
    void AddOption(const char **var, const char *short_name,
                   const char *long_name, const char *description,
                   bool required = false)
    {
        options.push_back(Option(STRING, var, short_name, long_name, description,
                                 required));
    }

    /// Add a string (std::string) option and set 'var' to receive the value.
    void AddOption(std::string *var, const char *short_name,
                   const char *long_name, const char *description,
                   bool required = false)
    {
        options.push_back(Option(STD_STRING, var, short_name, long_name, description,
                                 required));
    }

    /** Add an integer array (separated by spaces) option and set 'var' to
        receive the values. */
    void AddOption(std::vector<int> *var, const char *short_name,
                   const char *long_name, const char *description,
                   bool required = false)
    {
        options.push_back(Option(ARRAY, var, short_name, long_name, description,
                                 required));
    }

    /** Add a vector (doubles separated by spaces) option and set 'var' to
        receive the values. */
    void AddOption(std::vector<double> *var, const char *short_name,
                   const char *long_name, const char *description,
                   bool required = false)
    {
        options.push_back(Option(VECTOR, var, short_name, long_name, description,
                                 required));
    }

    /** @brief Parse the command-line options.
        Note that this function expects all the options provided through the
        command line to have a corresponding AddOption. In particular, this
        function cannot be used for partial parsing. */
    void Parse();

    /// Parse the command line options, and exit with an error if the options
    /// cannot be parsed successfully. The selected options are printed to the
    /// given stream (defaulting to mfem::out).
    void ParseCheck(std::ostream &out = std::cout);

    /// Return true if the command line options were parsed successfully.
    bool Good() const { return (error_type == 0); }

    /// Return true if we are flagged to print the help message.
    bool Help() const { return (error_type == 1); }

    /// Print the options
    void PrintOptions(std::ostream &out) const;

    /// Print the error message
    void PrintError(std::ostream &out) const;

    /// Print the help message
    void PrintHelp(std::ostream &out) const;

    /// Print the usage message
    void PrintUsage(std::ostream &out) const;

    int isValidAsInt(char *s)
    {
        if (s == NULL || *s == '\0')
        {
            return 0; // Empty string
        }

        if (*s == '+' || *s == '-')
        {
            ++s;
        }

        if (*s == '\0')
        {
            return 0; // sign character only
        }

        while (*s)
        {
            if (!isdigit(*s))
            {
                return 0;
            }
            ++s;
        }

        return 1;
    }

    int isValidAsDouble(char *s)
    {
        // A valid floating point number for atof using the "C" locale is formed by
        // - an optional sign character (+ or -),
        // - followed by a sequence of digits, optionally containing a decimal-point
        //   character (.),
        // - optionally followed by an exponent part (an e or E character followed by
        //   an optional sign and a sequence of digits).

        if (s == NULL || *s == '\0')
        {
            return 0; // Empty string
        }

        if (*s == '+' || *s == '-')
        {
            ++s;
        }

        if (*s == '\0')
        {
            return 0; // sign character only
        }

        while (*s)
        {
            if (!isdigit(*s))
            {
                break;
            }
            ++s;
        }

        if (*s == '\0')
        {
            return 1; // s = "123"
        }

        if (*s == '.')
        {
            ++s;
            while (*s)
            {
                if (!isdigit(*s))
                {
                    break;
                }
                ++s;
            }
            if (*s == '\0')
            {
                return 1; // this is a fixed point double s = "123." or "123.45"
            }
        }

        if (*s == 'e' || *s == 'E')
        {
            ++s;
            return isValidAsInt(s);
        }
        else
        {
            return 0; // we have encounter a wrong character
        }
    }

    int isValidAsChar(char *s)
    {
        if(s[0]>=0 && s[0] <= 128){
            return 1;
        }

        return 0;
    }

    void parseArray(char *str, std::vector<int> &var)
    {
        var.resize(0);
        std::stringstream input(str);
        int val;
        while (input >> val)
        {
            var.push_back(val);
        }
    }

    void parseVector(char *str, std::vector<double> &var)
    {
        int nentries = 0;
        double val;
        {
            std::stringstream input(str);
            while (input >> val)
            {
                ++nentries;
            }
        }

        var.resize(nentries);
        {
            nentries = 0;
            std::stringstream input(str);
            while (input >> val)
            {
                var[nentries++] = val;
            }
        }
    }
};



#endif
