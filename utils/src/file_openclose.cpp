/*
 * gcc -c file_openclose.c -Wall // to object file
 * Opening and closing files with control of existence
 * Possibly universaly useful
 * v1.1 Nov 2019
 * Author: J. Janek
 */

#include <cstdio>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <string>
#include "file_openclose.hpp"
#include "general_utils.hpp"

int my_fopen_r(FILE **fr, char *filename, char *mode)
/* Open a file to read with existence controll */
{
  FILE *file;
  errno = 0;
  if ((file = fopen(filename, mode)) == NULL)
  {
    print_warning(0, "File '" + std::string(filename) + "' cannot be opened.\n");
    if (errno != 0)
      perror("ERROR: ");
    return 1;
  }
  *fr = file;
  return 0;
}

int my_fclose(FILE **f, char *filename)
/* Close a file with error controll */
{
  errno = 0;
  if (fclose(*f) == EOF)
  {
    print_warning(0, "Error while closing '" + std::string(filename) + "'.\n");
    if (errno != 0)
      perror("ERROR: ");
    return 1;
  }
  print_warning(2, "File '"+ std::string(filename) + "' was successfuly closed.\n");
  *f = nullptr;
  return 0;
}

int my_fopen_w(FILE **fw, char *filename, char *mode)
/* Open a file to write, firstly asking if you want to
 * overwrite the existing file (if exists)
 */
{
  FILE *file, *backup;
  char backupfile[100];
  char c;

  if ((file = fopen(filename, "r")) != NULL)
  {
    strcpy(backupfile, filename);
    strcat(backupfile, "~");
    print_warning(2, "File '"+ std::string(filename) + "' already exists. Creating a backup '"+ std::string(backupfile) + "'.\n");
    if ((backup = fopen(backupfile, "w")) == NULL)
    {
      print_warning(1, "Cannot create backupfile '" + std::string(backupfile) + "'!\n");
    }
    else
    {
      while (c = fgetc(file), !feof(file))
      {
        fputc(c, backup);
      }
      fclose(backup);
    }
    if (my_fclose(&file, filename))
      return 1;
  }

  errno = 0;
  if ((file = fopen(filename, mode)) == NULL)
  {
    print_warning(0, "Error while opening '"+ std::string(filename) + "'.\n");
    if (errno != 0)
      perror("ERROR: ");
    return 2;
  }
  *fw = file;
  file = nullptr;
  print_warning(2, "File '"+ std::string(filename) + "' was successfuly opened to write.\n");
  return 0;
}

// copy files
int my_copy(char *srcname, char *destname)
{
  FILE *src, *dest, *backup;
  char backupfile[100];
  char c;

  if ((src = fopen(srcname, "r")) != NULL)
  {
    if ((dest = fopen(destname, "r")) != NULL)
    {
      strcpy(backupfile, destname);
      strcat(backupfile, "~");
      print_warning(2, "File '"+ std::string(destname) + "' already exists. Creating a backup '"+ std::string(backupfile) + "'.\n");
      if ((backup = fopen(backupfile, "w")) == NULL)
      {
        print_warning(1, "Cannot create backupfile '" + std::string(backupfile) + "'\n");
      }
      else
      {
        while (c = fgetc(dest), !feof(dest))
        {
          fputc(c, backup);
        }
        fclose(backup);
      }
      if (fclose(dest) != 0)
        return 1;
    }
    if ((dest = fopen(destname, "w")) != NULL)
    {
      while (c = fgetc(src), !feof(src))
      {
        fputc(c, dest);
      }
      if (fclose(dest) != 0)
        return 1;
    }
    else
    {
      return 2;
    }
    fclose(src);
  }
  else
  {
    return 3;
  }
  return 0;
}
