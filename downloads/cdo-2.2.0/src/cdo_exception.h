#ifndef CDO_EXCEPTION_H
#define CDO_EXCEPTION_H

#include <string>
#include <stdexcept>
struct CdoException : std::logic_error
{
  CdoException(std::string p_msg,std::string p_file, std::string p_line) : std::logic_error(p_msg), file(p_file), line(p_line) {}
  std::string file;
  std::string line;
};

#endif
