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
                                                                        
// $Revision: 1.9 $
// $Date: 2010-06-01 23:46:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/actor/MovableObject.cpp,v $
                                                                        
                                                                        
// Written: fmk
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the implementation of MovableObject.
//
// What: "@(#) MovableObject.C, revA"

#include <MovableObject.h>
#include <OPS_Globals.h>
#include <profiler/ProfilerMacros.h>   // Ladruno P4: TaggedObject live-component census

MovableObject::MovableObject(int cTag, int dTag)
:classTag(cTag), dbTag(dTag)
{
  OPS_PROFILE_CENSUS_BORN(classTag);   // Ladruno P4
}


MovableObject::MovableObject(int theTag)
:classTag(theTag), dbTag(0)
{
  OPS_PROFILE_CENSUS_BORN(classTag);   // Ladruno P4
}

MovableObject::~MovableObject()
{
  OPS_PROFILE_CENSUS_DIED(classTag);   // Ladruno P4
}




int 
MovableObject::getClassTag(void) const
{
    return classTag;
}

static char unknownClassType[] = {"UnknownMovableObject"};

const char *
MovableObject::getClassType(void) const
{
  return unknownClassType;
}


int 
MovableObject::getDbTag(void) const
{
    return dbTag;
}

void
MovableObject::setDbTag(int newTag)
{
  dbTag = newTag;
}

int
MovableObject::setParameter(const char **argv, int argc, Parameter &param)
{
  return -1;
}

int
MovableObject::updateParameter(int parameterID, Information &info)
{
  return 0;
}

int
MovableObject::activateParameter(int parameterID)
{
  return 0;
}

int 
MovableObject::setVariable(const char *variable, Information &theInfo)
{
  return -1;
}

int 
MovableObject::getVariable(const char *variable, Information &theInfo)
{
  return -1;
}
