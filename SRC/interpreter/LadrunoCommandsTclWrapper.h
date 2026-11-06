// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

// Ladruno WP-168 (apeGmsh#1490): Ladruno_registerCommands for the DL Tcl engine
// (TclWrapper). TclWrapper::addOpenSeesCommands calls it once; it registers
// every row of LadrunoCommandTable.h whose `dl` column is not LADRUNO_NONE.
//
// The bridge is generated, one instantiation per row, and is exactly the shape
// every hand-written Tcl_ops_Ladruno* bridge had before WP-168:
//
//     wrapper->resetCommandLine(argc, 1, argv);
//     if (OPS_X() < 0) return TCL_ERROR;
//     return TCL_OK;
//
// Include this from TclWrapper.cpp only.

#ifndef LadrunoCommandsTclWrapper_h
#define LadrunoCommandsTclWrapper_h

#include "TclWrapper.h"
#include "OpenSeesCommands.h"

namespace ladruno_tclwrapper_commands {

// The wrapper that registered the commands. TclWrapper.cpp's own bridges use its
// file-static `wrapper`, which the constructor sets to `this`; the hook is called
// from that same object's addOpenSeesCommands, so this holds the same pointer.
inline TclWrapper*& host()
{
    static TclWrapper* theHost = 0;
    return theHost;
}

template <int (*F)(void)>
int command(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
    host()->resetCommandLine(argc, 1, argv);
    if (F() < 0) return TCL_ERROR;
    return TCL_OK;
}

template <int (*F)(void)>
void add(TclWrapper* w, Tcl_Interp* interp, const char* name)
{
    if constexpr (F != nullptr)
        w->addCommand(interp, name, &command<F>);
}

}  // namespace ladruno_tclwrapper_commands

inline void Ladruno_registerCommands(TclWrapper* w, Tcl_Interp* interp)
{
    ladruno_tclwrapper_commands::host() = w;
#define LADRUNO_NONE nullptr
#define LADRUNO_COMMAND(name, dl, classic) \
    ladruno_tclwrapper_commands::add<dl>(w, interp, name);
#include "LadrunoCommandTable.h"
#undef LADRUNO_COMMAND
#undef LADRUNO_NONE
}

#endif
