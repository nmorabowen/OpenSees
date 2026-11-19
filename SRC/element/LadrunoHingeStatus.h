/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

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

// Ladruno (WP-235) — the `hingeStatus` element response, shared by
// LadrunoForceBeamColumn, LadrunoDispBeamColumn and LadrunoIMKBeam (2-D and 3-D).
//
// `hingeSummary` answers "what rotation and moment", and reads 0 both for an
// unloaded hinge and for no hinge at all. Damage needs a different token per
// element and hinge type. `hingeStatus` answers the other two questions with one
// token and one layout everywhere:
//
//     kindZ_i, damageZ_i, kindZ_j, damageZ_j, kindY_i, damageY_i, kindY_j, damageY_j
//
// The slot order is hingeSummary's. A 2-D element has the same 8 columns; its Y
// slots are always absent, so a reader needs no per-dimension branch.
//
//   kind    0 = no hinge in this slot
//           1 = a uniaxial law owns this slot
//           2 = a biaxial law owns this slot. One law serves BOTH the Z and the Y
//               slot of that end, so both slots read kind 2 and the SAME damage.
//               It is one damage value, not two: never sum or count them apart.
//
//   damage  the scalar the law returns for its `damage` response, unchanged.
//           NaN when the slot is absent, when the law has no `damage` response
//           (Steel01, Bilin, ...), or when that response is not one scalar
//           (snap Bilinear returns three). NaN means "not defined"; 0 always
//           means a law reported zero damage.
//
// Each element fills the vector from its live material pointers on every call.
// Nothing is cached, because recvSelf can replace a hinge material.

#ifndef LadrunoHingeStatus_h
#define LadrunoHingeStatus_h

#include <limits>

#include <DummyStream.h>
#include <Information.h>
#include <NDMaterial.h>
#include <OPS_Stream.h>
#include <Response.h>
#include <UniaxialMaterial.h>
#include <Vector.h>

namespace LadrunoHingeStatus {

enum { numColumns = 8 };   // size of the response vector
enum { kindAbsent = 0, kindUniaxial = 1, kindBiaxial = 2 };
enum { axisZ = 0, axisY = 1 };
enum { endI = 0, endJ = 1 };

// Response id, free on all six element classes.
enum { responseID = 36 };

inline double notDefined()
{
  return std::numeric_limits<double>::quiet_NaN();
}

// The names reach COLUMN_MAP/COMP_NAMES of a .ladruno file through these tags.
inline void tagColumns(OPS_Stream &output)
{
  output.tag("ResponseType", "kindZ_i");
  output.tag("ResponseType", "damageZ_i");
  output.tag("ResponseType", "kindZ_j");
  output.tag("ResponseType", "damageZ_j");
  output.tag("ResponseType", "kindY_i");
  output.tag("ResponseType", "damageY_i");
  output.tag("ResponseType", "kindY_j");
  output.tag("ResponseType", "damageY_j");
}

// Ask a material for its own `damage` response, the one `hingeI damage` and
// `hinge 0 damage` reach. The probe gets a DummyStream: a material's setResponse
// may emit tags, and those must not become columns of this response.
template <class Material>
inline double damageOf(Material *material)
{
  const char *argv[1] = {"damage"};
  DummyStream silent;
  Response *response = material->setResponse(argv, 1, silent);
  if (response == 0)
    return notDefined();

  double damage = notDefined();
  if (response->getResponse() >= 0) {
    const Vector &data = response->getInformation().getData();
    if (data.Size() == 1)
      damage = data(0);
  }
  delete response;
  return damage;
}

inline void setSlot(Vector &status, int axis, int end, int kind, double damage)
{
  const int column = 4 * axis + 2 * end;
  status(column) = kind;
  status(column + 1) = damage;
}

// Every slot absent.
inline void reset(Vector &status)
{
  for (int axis = axisZ; axis <= axisY; ++axis)
    for (int end = endI; end <= endJ; ++end)
      setSlot(status, axis, end, kindAbsent, notDefined());
}

// Fill one end. A biaxial law takes precedence over uniaxial laws at the same
// end, as it does in hingeSummary. Pass 0 for whatever the element lacks.
inline void fillEnd(Vector &status, int end, UniaxialMaterial *hingeZ,
                    UniaxialMaterial *hingeY, NDMaterial *hingeBiaxial)
{
  if (hingeBiaxial != 0) {
    const double damage = damageOf(hingeBiaxial);
    setSlot(status, axisZ, end, kindBiaxial, damage);
    setSlot(status, axisY, end, kindBiaxial, damage);
    return;
  }
  if (hingeZ != 0)
    setSlot(status, axisZ, end, kindUniaxial, damageOf(hingeZ));
  if (hingeY != 0)
    setSlot(status, axisY, end, kindUniaxial, damageOf(hingeY));
}

} // namespace LadrunoHingeStatus

#endif
