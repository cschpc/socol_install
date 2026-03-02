/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef LISTBUF_H
#define LISTBUF_H

#include <stdio.h>
#include <sys/stat.h>
#include <string>
#include <vector>

class ListBuffer
{
public:
  std::vector<char> buffer;
  std::string name;

  int
  read(FILE *fp, const char *_name)
  {
    int filedes = fileno(fp);
    struct stat buf;
    const size_t filesize = (fstat(filedes, &buf) == 0) ? (size_t) buf.st_size : 0;
    if (filesize == 0)
      {
        fprintf(stderr, "ListBuffer: empty stream: %s\n", _name);
        return -1;
      }

    this->buffer.resize(filesize);
    if (fread(buffer.data(), 1, filesize, fp) != filesize)
      {
        fprintf(stderr, "ListBuffer: read failed on %s!\n", _name);
        return -1;
      }

    if (_name) this->name = _name;

    return 0;
  }
};

#endif
