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

#ifndef nonlocmatstiffinterface_h
#define nonlocmatstiffinterface_h

#include "interface.h"
#include "nonlocalmaterialext.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
class SparseMtrx;
class GaussPoint;
class TimeStep;

/**
 * Class Nonlocal Material Stiffness Interface. This is only abstract class.
 * This interface allows material model to add those services required to
 * compute and assemble nonlocal contribution to stiffness matrix.
 */
class NonlocalMaterialStiffnessInterface : public Interface
{
public:
    /// Constructor.
    NonlocalMaterialStiffnessInterface() : Interface() { }

    /// Computes and adds IP contributions to destination matrix.
    virtual void NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s,
                                                                      GaussPoint *gp, TimeStep *atTime) = 0;
    /**
     * Returns integration list of receiver. Contains localIntegrationRecord structures, containing
     * references to integration points and their weights that influence to nonlocal average in
     * receiver's associated integration point.
     */
    virtual dynaList< localIntegrationRecord > *NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp) = 0;

#ifdef __OOFEG
    /**
     * Plots the sparse structure of stiffness contribution.
     */
    virtual void NonlocalMaterialStiffnessInterface_showSparseMtrxStructure(GaussPoint *gp, oofegGraphicContext &gc, TimeStep *) { };
#endif
};
} // end namespace oofem
#endif // nonlocmatstiffinterface_h
