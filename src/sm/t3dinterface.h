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

#ifndef t3dinterface_h
#define t3dinterface_h

#include "mesherinterface.h"

namespace oofem {
#define BMF_FILENAME "t3d.bmf"

class TimeStep;

/**
 * This class represents the interface to t3d mesh generation package.
 * This interface is primarily responsible for two main tasks:
 * - to create input mesher file, containing all information including the mesh density information
 *   based on information from remeshing criteria.
 * - possibly to launch the mesher and transform its output to oofem input (using t3d2oofem)
 */
class T3DInterface : public MesherInterface
{
public:
    /// Constructor
    T3DInterface(Domain *d) : MesherInterface(d) { }
    /// Destructor
    virtual ~T3DInterface() { }

    virtual returnCode createMesh(TimeStep *tStep, int domainNumber, int domainSerNum, Domain **dNew);

protected:
    /// Creates the mesher input, containing the required mesh density information.
    int createInput(Domain *d, TimeStep *stepN);
};
} // end namespace oofem
#endif // t3dinterface_h
