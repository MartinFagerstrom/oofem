/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef bodyload_h
#define bodyload_h

#include "load.h"

#include "classtype.h"
#include "bcgeomtype.h"

namespace oofem {

/**
 * Class implementing element body load, acting over whole element volume (e.g., the dead weight).
 * Body load is usually attribute of one or more elements.
 *
 * This base body load class only defines the common services common to all derived classes.
 * Derived classes need to implement services declared by base Load class.
 */
class BodyLoad : public Load
{
public:
    /**
     * Constructor. Creates Body Load object with given number, belonging to given domain.
     * @param i Load number
     * @param d Domain to which new load belongs.
     */
    BodyLoad(int i, Domain *d) : Load(i, d) { }

    /**
     * Returns receiver's load geometry type.
     * @return BodyLoadBGT.
     */
    virtual bcGeomType giveBCGeoType() const { return BodyLoadBGT; }

    virtual const char *giveClassName() const { return "BodyLoad"; }
    virtual classType giveClassID() const { return BodyLoadClass; }
};
} // end namespace oofem
#endif // bodyload_h








