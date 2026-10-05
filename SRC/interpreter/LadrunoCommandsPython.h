// Ladruno WP-168 (apeGmsh#1490): Ladruno_registerCommands for the Python engine
// (PythonWrapper: opensees.pyd and openseesmp.pyd). PythonWrapper::
// addOpenSeesCommands calls it once; it registers every row of
// LadrunoCommandTable.h whose `dl` column is not LADRUNO_NONE.
//
// The bridge is generated, one instantiation per row, and is exactly the shape
// every hand-written Py_ops_Ladruno* bridge had before WP-168:
//
//     wrapper->resetCommandLine(PyTuple_Size(args), 1, args);
//     if (OPS_X() < 0) { opserr<<(void*)0; return NULL; }
//     return wrapper->getResults();
//
// Include this from PythonWrapper.cpp only.

#ifndef LadrunoCommandsPython_h
#define LadrunoCommandsPython_h

#include "PythonWrapper.h"
#include "OpenSeesCommands.h"
#include <OPS_Globals.h>

namespace ladruno_python_commands {

// The wrapper that registered the commands. PythonWrapper.cpp's own bridges use
// its file-static `wrapper`, which the constructor sets to `this`; the hook is
// called from that same object's addOpenSeesCommands, so this holds the same
// pointer.
inline PythonWrapper*& host()
{
    static PythonWrapper* theHost = 0;
    return theHost;
}

template <int (*F)(void)>
PyObject* command(PyObject* self, PyObject* args)
{
    host()->resetCommandLine(PyTuple_Size(args), 1, args);
    if (F() < 0) {
        opserr << (void*)0;
        return NULL;
    }
    return host()->getResults();
}

template <int (*F)(void)>
void add(PythonWrapper* w, const char* name)
{
    if constexpr (F != nullptr)
        w->addCommand(name, &command<F>);
}

}  // namespace ladruno_python_commands

inline void Ladruno_registerCommands(PythonWrapper* w)
{
    ladruno_python_commands::host() = w;
#define LADRUNO_NONE nullptr
#define LADRUNO_COMMAND(name, dl, classic) \
    ladruno_python_commands::add<dl>(w, name);
#include "LadrunoCommandTable.h"
#undef LADRUNO_COMMAND
#undef LADRUNO_NONE
}

#endif
