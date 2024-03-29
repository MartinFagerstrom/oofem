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

#ifndef freeminterface_h
#define freeminterface_h

#include "mesherinterface.h"
#include "flotarry.h"

namespace oofem {
class TimeStep;

/**
 * This class represents the interface to freem mesh generation package.
 * This interface is primarily responsible for two main tasks:
 * - To create input mesher file, containing all information including the mesh density information
 *   based on the remeshing criteria.
 * - Possibly to launch the mesher and transform its output to oofem input.
 */
class FreemInterface : public MesherInterface
{
public:
    /// Constructor
    FreemInterface(Domain *d) : MesherInterface(d) { }
    /// Destructor
    virtual ~FreemInterface() { }

    virtual returnCode createMesh(TimeStep *tStep, int domainNumber, int domainSerNum, Domain **dNew);

protected:
    /// Creates the mesher input, containing the required mesh density information.
    int createInput(Domain *d, TimeStep *stepN);
    /// Service for smoothing the densities for freem.
    void smoothNodalDensities(Domain *d, FloatArray &nodalDensities, TimeStep *stepN);
};
} // end namespace oofem
#endif // freeminterface_h
