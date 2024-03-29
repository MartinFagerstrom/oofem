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

#ifndef eigenvaluedynamic_h
#define eigenvaluedynamic_h

#include "engngm.h"
#include "sparsegeneigenvalsystemnm.h"
#include "sparsemtrx.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "geneigvalsolvertype.h"

namespace oofem {
/**
 * This class implements way for examining eigenvalues and eigenvectors in
 * dynamic problems.
 *
 * Solution of this problem is base on equation in the form of: @f$ K\cdot y=w M\cdot y @f$
 * Currently eigenvalue problem is solved using subspace iteration.
 * Tasks:
 * - Assembling the governing equation in the form @f$ K\cdot y=wM\cdot y@f$.
 * - Creating Numerical method for @f$ K\cdot y=wM\cdot y@f$.
 * - Interfacing Numerical method to Elements.
 */
class EigenValueDynamic : public EngngModel
{
private:
    SparseMtrx *stiffnessMatrix;
    SparseMtrx *massMatrix;
    SparseMtrxType sparseMtrxType;
    FloatMatrix eigVec;
    FloatArray eigVal;
    int numberOfRequiredEigenValues;
    int activeVector;
    int restoreFlag;
    /// Relative tolerance.
    double rtolv;
    /// Numerical method used to solve the problem.
    SparseGeneralEigenValueSystemNM *nMethod;
    GenEigvalSolverType solverType;

public:
    EigenValueDynamic(int i, EngngModel *_master = NULL) : EngngModel(i, _master)
    {
        stiffnessMatrix = NULL;
        massMatrix = NULL;
        numberOfSteps = 1;
        ndomains = 1;
        nMethod = NULL;
    }
    virtual ~EigenValueDynamic() {
        delete  stiffnessMatrix;
        delete massMatrix;
        if ( nMethod ) { delete nMethod; } }

    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void terminate(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);

    virtual double giveUnknownComponent(EquationID eid, ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    virtual double giveUnknownComponent(UnknownType ut, ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);
    virtual void setActiveVector(int i) { activeVector = i; }
    virtual int resolveCorrespondingEigenStepNumber(void *obj);

#ifdef __SLEPC_MODULE
    virtual void initPetscContexts();
#endif

    /**
     * DOF printing routine. Called by DofManagers to print Dof specific part.
     * Dof class provides component printing routines, but emodel is responsible
     * for what will be printed at DOF level.
     * @param stream Output stream.
     * @param iDof Dof to be processed.
     * @param atTime Solution step.
     */
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

    // identification
    virtual const char *giveClassName() const { return "EigenValueDynamic"; }
    virtual classType giveClassID() const { return EigenValueDynamicClass; }
};
} // end namespace oofem
#endif // eigenvaluedynamic_h
