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

#ifndef deformationtheorymaterial_h
#define deformationtheorymaterial_h

#include "structuralmaterial.h"

namespace oofem {

/**
 * This class implements an abstract base material, which behaves
 * according to deformation theory.
 */
class DeformationTheoryMaterial : public StructuralMaterial
{
public:
    DeformationTheoryMaterial(int n, Domain *d) : StructuralMaterial(n, d) { }
    virtual ~DeformationTheoryMaterial()  { }

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour() { return 1; }
    virtual const char *giveClassName() const { return "DeformationTheoryMaterial"; }
    virtual classType giveClassID() const { return DeformationTheoryMaterialClass; }
};
} // end namespace oofem
#endif // deformationtheorymaterial_h
