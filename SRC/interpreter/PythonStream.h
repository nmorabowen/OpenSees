/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision$
// $Date$

#ifndef _PythonStream
#define _PythonStream

#include <StandardStream.h>
#include <Python.h>
#include <sstream>
#include <cstdio>   // Ladruno WP-258: fputs to stderr off the GIL

class PythonStream : public StandardStream {
public:
  PythonStream(int indentSize = 2, bool echo = true, bool standard_echo = false);


  ~PythonStream();

  void setError(PyObject *err);
  int setFile(const char *fileName, openMode mode, bool echo);
  OPS_Stream &operator<<(char c);
  OPS_Stream &operator<<(unsigned char c);
  OPS_Stream &operator<<(signed char c);
  OPS_Stream &operator<<(const char *s);
  OPS_Stream &operator<<(const unsigned char *s);
  OPS_Stream &operator<<(const signed char *s);
  OPS_Stream &operator<<(int n);
  OPS_Stream &operator<<(unsigned int n);
  OPS_Stream &operator<<(long n);
  OPS_Stream &operator<<(unsigned long n);
  OPS_Stream &operator<<(short n);
  OPS_Stream &operator<<(unsigned short n);
  OPS_Stream &operator<<(bool b);
  OPS_Stream &operator<<(double n);
  OPS_Stream &operator<<(float n);
  OPS_Stream &operator<<(const void *p);

private:
  template<class T>void err_out(T err);
  PyObject *error;
  std::string msg;
  bool echoApplication;
};

template<class T>void PythonStream::err_out(T err)
{
  std::stringstream ss;
  ss << err;
  // Ladruno WP-258: a LOCAL string (the `msg` member was one buffer shared by
  // every thread that prints), and no call into CPython from a thread that does
  // not hold the GIL. An opserr from an OpenMP worker of a threaded element or
  // particle loop reached PySys_FormatStderr with no Python thread state and
  // dereferenced null -- the WP-107 "ManzariDafalias IntScheme 1 segfaults on a
  // worker" (it faulted at the first plastic-branch warning; the elastic branch
  // never prints). Taking the GIL there would deadlock: the main thread holds it
  // for the whole parallel region. Such a thread writes to the C stderr instead
  // (the same fd 2 that sys.stderr wraps in a normal session); the thread that
  // holds the GIL -- every serial print -- keeps the Python path unchanged.
  const std::string text = ss.str();
  if (!PyGILState_Check()) {
    fputs(text.c_str(), stderr);
    fflush(stderr);
    return;
  }
  // Ladruno: the already-formatted message must be passed as a %s ARGUMENT, not as
  // the format string itself. Any literal '%' in an opserr message (e.g. "50% of
  // model mass", "exceeds -maxAddedMass cap 5%") was otherwise consumed as a bogus
  // printf conversion by PySys_FormatStderr and silently dropped/garbled under
  // openseespy (the Tcl StandardStream path was unaffected). Classic format-string bug.
  PySys_FormatStderr("%s", text.c_str());
}


#endif
