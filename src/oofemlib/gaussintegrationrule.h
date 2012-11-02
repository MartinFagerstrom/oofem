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

#ifndef gaussintegrationrule_h
#define gaussintegrationrule_h

#include "integrationrule.h"
#include "element.h"

namespace oofem {
/**
 * Class representing Gaussian-quadrature integration rule.
 * The number of integration points and their coordinates and integration weights depends on
 * integration rule type (rule for integration in 1d, 2d, 3d) and required accuracy.
 * 
 * Tasks:
 * - Returning number of integration points used
 * - Returning requested integration point
 * - Updating itself
 * - Saving and restoring context
 * 
 * @see GaussPoint
 */
class GaussIntegrationRule : public IntegrationRule
{
public:
    /**
     * Constructor.
     * @param n Number associated with receiver.
     * @param e Element associated with receiver.
     * @param startIndx First component, for which rule applies.
     * @param endIndx Last component, for which rule applies.
     * @param dynamic Flag indicating that receiver can change.
     */
    GaussIntegrationRule(int n, Element *e, int startIndx, int endIndx, bool dynamic = false);
    GaussIntegrationRule(int n, Element *e);
    /// Destructor
    virtual ~GaussIntegrationRule();

    virtual classType giveClassID() const { return GaussIntegrationRuleClass; }
    virtual const char *giveClassName() const { return "GaussIntegrationRule"; }
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }

    /**
     * Returns requred number of integration points to exactly integrate
     * polynomial of order approxOrder on given domain.
     * When approxOrder is too large and is not supported by implementation
     * method returns -1.
     */
    int getRequiredNumberOfIntegrationPoints(integrationDomain dType, int approxOrder);

protected:
    virtual int SetUpPointsOnLine(int, MaterialMode, GaussPoint * * *);
    virtual int SetUpPointsOnTriangle(int, MaterialMode, GaussPoint * * *);
    virtual int SetUpPointsOnSquare(int, MaterialMode, GaussPoint * * *);
    virtual int SetUpPointsOnCube(int, MaterialMode, GaussPoint * * *);
    virtual int SetUpPointsOnTetrahedra(int, MaterialMode, GaussPoint * * *);
    virtual int SetUpPointsOnWedge(int, MaterialMode, GaussPoint * * *);
    virtual int SetUpPointsOn2DEmbeddedLine(int nPoints, MaterialMode mode, GaussPoint ***,
                                    const FloatArray **coords);
};
} // end namespace oofem
#endif // gaussintegrationrule_h
