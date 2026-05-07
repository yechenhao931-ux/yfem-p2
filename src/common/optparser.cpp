#include <common/optparser.hpp>
#include <cctype>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <type_traits>
#include <initializer_list>
#include <sstream>
#include <string>
using namespace std;

void OptionsParser::Parse()
{
   option_check.resize(options.size());
   option_check.assign(option_check.size(),0);
   for (int i = 1; i < argc; )
   {
      if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
      {
         // print help message
         error_type = 1;
         return;
      }

      for (int j = 0; true; j++)
      {
         if (j >= options.size())
         {
            // unrecognized option
            error_type = 2;
            error_idx = i;
            return;
         }

         if (strcmp(argv[i], options[j].short_name) == 0 ||
             strcmp(argv[i], options[j].long_name) == 0)
         {
            OptionType type = options[j].type;

            if ( option_check[j] )
            {
               error_type = 4;
               error_idx = j;
               return;
            }
            option_check[j] = 1;

            i++;
            if (type != ENABLE && type != DISABLE && i >= argc)
            {
               // missing argument
               error_type = 3;
               error_idx = j;
               return;
            }

            int isValid = 1;
            switch (options[j].type)
            {
               case INT:
                  isValid = isValidAsInt(argv[i]);
                  *(int *)(options[j].var_ptr) = atoi(argv[i++]);
                  break;
               case DOUBLE:
                  isValid = isValidAsDouble(argv[i]);
                  *(double *)(options[j].var_ptr) = atof(argv[i++]);
                  break;
               case CHAR:
                  isValid = isValidAsChar(argv[i]);
                  *(char *)(options[j].var_ptr) = argv[i++][0];
                  break;
               case STRING:
                  *(const char **)(options[j].var_ptr) = argv[i++];
                  break;
               case STD_STRING:
                  *(std::string *)(options[j].var_ptr) = argv[i++];
                  break;
               case ENABLE:
                  *(bool *)(options[j].var_ptr) = true;
                  option_check[j+1] = 1;  // Do not allow the DISABLE Option
                  break;
               case DISABLE:
                  *(bool *)(options[j].var_ptr) = false;
                  option_check[j-1] = 1;  // Do not allow the ENABLE Option
                  break;
               case ARRAY:
                  parseArray(argv[i++], *(std::vector<int> *)(options[j].var_ptr) );
                  break;
               case VECTOR:
                  parseVector(argv[i++], *(std::vector<double> *)(options[j].var_ptr) );
                  break;
            }

            if (!isValid)
            {
               error_type = 5;
               error_idx = i;
               return;
            }

            break;
         }
      }
   }

   // check for missing required options
   for (int i = 0; i < options.size(); i++)
      if (options[i].required &&
          (option_check[i] == 0 ||
           (options[i].type == ENABLE && option_check[++i] == 0)))
      {
         error_type = 6; // required option missing
         error_idx = i; // for a boolean option i is the index of DISABLE
         return;
      }

   error_type = 0;
}

void OptionsParser::ParseCheck(std::ostream &os)
{
   Parse();
   int my_rank = 0;
   if (!Good())
   {
      if (my_rank == 0) { PrintUsage(os); }
      std::exit(1);
   }
   if (my_rank == 0) { PrintOptions(os); }
}

void OptionsParser::WriteValue(const Option &opt, std::ostream &os)
{
   switch (opt.type)
   {
      case INT:
         os << *(int *)(opt.var_ptr);
         break;

      case DOUBLE:
         os << *(double *)(opt.var_ptr);
         break;

      case STRING:
         os << *(const char **)(opt.var_ptr);
         break;

      case STD_STRING:
         os << *(std::string *)(opt.var_ptr);
         break;

      case ARRAY:
      {
         std::vector<int> &list = *(std::vector<int>*)(opt.var_ptr);
         os << '\'';
         if (list.size() > 0)
         {
            os << list[0];
         }
         for (int i = 1; i < list.size(); i++)
         {
            os << ' ' << list[i];
         }
         os << '\'';
         break;
      }

      case VECTOR:
      {
         std::vector<double> &list = *(std::vector<double>*)(opt.var_ptr);
         os << '\'';
         if (list.size() > 0)
         {
            os << list[0];
         }
         for (int i = 1; i < list.size(); i++)
         {
            os << ' ' << list[i];
         }
         os << '\'';
         break;
      }

      default: // provide a default to suppress warning
         break;
   }
}

void OptionsParser::PrintOptions(std::ostream &os) const
{
   static const char *indent = "   ";

   os << "Options used:\n";
   for (int j = 0; j < options.size(); j++)
   {
      OptionType type = options[j].type;

      os << indent;
      if (type == ENABLE)
      {
         if (*(bool *)(options[j].var_ptr) == true)
         {
            os << options[j].long_name;
         }
         else
         {
            os << options[j+1].long_name;
         }
         j++;
      }
      else
      {
         os << options[j].long_name << " ";
         WriteValue(options[j], os);
      }
      os << '\n';
   }
}

void OptionsParser::PrintError(std::ostream &os) const
{
   static const char *line_sep = "";

   os << line_sep;
   switch (error_type)
   {
      case 2:
         os << "Unrecognized option: " << argv[error_idx] << '\n' << line_sep;
         break;

      case 3:
         os << "Missing argument for the last option: " << argv[argc-1]
            << '\n' << line_sep;
         break;

      case 4:
         if (options[error_idx].type == ENABLE )
            os << "Option " << options[error_idx].long_name << " or "
               << options[error_idx + 1].long_name
               << " provided multiple times\n" << line_sep;
         else if (options[error_idx].type == DISABLE)
            os << "Option " << options[error_idx - 1].long_name << " or "
               << options[error_idx].long_name
               << " provided multiple times\n" << line_sep;
         else
            os << "Option " << options[error_idx].long_name
               << " provided multiple times\n" << line_sep;
         break;

      case 5:
         os << "Wrong option format: " << argv[error_idx - 1] << " "
            << argv[error_idx] << '\n' << line_sep;
         break;

      case 6:
         os << "Missing required option: " << options[error_idx].long_name
            << '\n' << line_sep;
         break;
   }
   os << std::endl;
}

void OptionsParser::PrintHelp(std::ostream &os) const
{
   static const char *indent = "   ";
   static const char *seprtr = ", ";
   static const char *descr_sep = "\n\t";
   static const char *line_sep = "";
   static const char *types[] = { " <int>", " <double>", " <string>",
                                  " <string>", "", "", " '<int>...'",
                                  " '<double>...'"
                                };

   os << indent << "-h" << seprtr << "--help" << descr_sep
      << "Print this help message and exit.\n" << line_sep;
   for (int j = 0; j < options.size(); j++)
   {
      OptionType type = options[j].type;

      os << indent << options[j].short_name << types[type]
         << seprtr << options[j].long_name << types[type]
         << seprtr;
      if (options[j].required)
      {
         os << "(required)";
      }
      else
      {
         if (type == ENABLE)
         {
            j++;
            os << options[j].short_name << types[type] << seprtr
               << options[j].long_name << types[type] << seprtr
               << "current option: ";
            if (*(bool *)(options[j].var_ptr) == true)
            {
               os << options[j-1].long_name;
            }
            else
            {
               os << options[j].long_name;
            }
         }
         else
         {
            os << "current value: ";
            WriteValue(options[j], os);
         }
      }
      os << descr_sep;

      if (options[j].description)
      {
         os << options[j].description << '\n';
      }
      os << line_sep;
   }
}

void OptionsParser::PrintUsage(std::ostream &os) const
{
   static const char *line_sep = "";

   PrintError(os);
   os << "Usage: " << argv[0] << " [options] ...\n" << line_sep
      << "Options:\n" << line_sep;
   PrintHelp(os);
}
