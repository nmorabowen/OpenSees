// Ladruno WP-168 (apeGmsh#1490): Ladruno_registerCommands for classic Tcl, the
// engine OpenSees.exe, OpenSeesSP and OpenSeesMP run (SRC/tcl/commands.cpp).
// OpenSeesAppInit calls it once; it registers every row of
// SRC/interpreter/LadrunoCommandTable.h whose `classic` column is not
// LADRUNO_NONE, with the same Tcl_CreateCommand call (NULL ClientData, no
// delete proc) the hand-written registrations used.
//
// The classic bridges themselves stay in commands.cpp: several bind this
// engine's own state (see the table's header comment). Include this from
// commands.cpp only, after every bridge the table names is declared.

#ifndef LadrunoCommandsClassicTcl_h
#define LadrunoCommandsClassicTcl_h

#include <tcl.h>

static void Ladruno_registerCommands(Tcl_Interp* interp)
{
    struct Row {
        const char* name;
        Tcl_CmdProc* proc;
    };
    static const Row rows[] = {
#define LADRUNO_NONE nullptr
#define LADRUNO_COMMAND(name, dl, classic) {name, classic},
#include "../interpreter/LadrunoCommandTable.h"
#undef LADRUNO_COMMAND
#undef LADRUNO_NONE
    };
    for (const Row& r : rows)
        if (r.proc != nullptr)
            Tcl_CreateCommand(interp, r.name, r.proc,
                              (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
}

#endif
