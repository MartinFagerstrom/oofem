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

#ifndef macrolspace_h
#define macrolspace_h

#include "lspace.h"
#include "sparsemtrx.h"
#include "engngm.h"
#include "structengngmodel.h"
//#include "micromaterial.h"
#include "metastep.h"
#include "oofem_limits.h"
//#include "nlstructuralelement.h"
//#include "fei3dhexalin.h"

namespace oofem {
class MicroMaterial;

/**
 * This class implements a macroelement. It is derived from eight-node brick element.
 * The stiffness matrix is computed from underlying RVE and is condensed to 24 DoFs to corner nodes.
 *
 * @author Vit Smilauer
 */
class MacroLSpace : public LSpace
{
public:
    MacroLSpace(int n, Domain *d);
    virtual ~MacroLSpace();

    virtual const char *giveClassName() const { return "MacroLSpace"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    /// Related to setting the boundary conditions of micro problem.
    virtual void changeMicroBoundaryConditions(TimeStep *tStep);

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);

    /**
     * Evaluates shape function at a given pointnodal representation of real internal forces obtained from microProblem.
     * @param answer Array of shape function values at given node.
     * @param coords Coordinates of nodes defining the interpolation geometry.
     * @param gcoords Global coordinates of point of interest.
     */
    virtual void evalInterpolation(FloatArray &answer, const FloatArray **coords, const FloatArray &gcoords);

    virtual void updateYourself(TimeStep *tStep);

protected:
    /// Array containing the node mapping from microscale (which microMasterNodes corresponds to which macroNode)
    IntArray microMasterNodes;
    IntArray microBoundaryNodes;
    IntArray microDOFs;
    bool firstCall;
    MicroMaterial *microMaterial;
    Domain *microDomain;
    EngngModel *microEngngModel;
    /// Information of iteration number.
    int iteration;
    /// Stores node number on the boundary in the triplets.
    IntArray microBoundaryDofManager;
    FloatMatrix stiffMatrix;
    /// Process with external file for the storage of stiffness matrix 0-None, 1-read, 2-write.
    int stiffMatrxFileNoneReadingWriting;
    /// Array containg the force vector from nodes (if condensation is skipped, use this vector).
    FloatArray internalMacroForcesVector;
    /// Last time step when stiffness matrix was assembled.
    TimeStep *lastStiffMatrixTimeStep;
};
} // end namespace oofem
#endif //macrolspace_h
