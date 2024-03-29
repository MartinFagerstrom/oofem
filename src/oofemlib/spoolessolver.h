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
#ifndef spoolessolver_h
#define spoolessolver_h

#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "flotarry.h"

#include "spoolesinterface.h"

namespace oofem {
class Domain;
class EngngModel;
class FloatMatrix;

/**
 * Implements the solution of linear system of equation in the form @f$ A\cdot x = b @f$ using solvers
 * from SPOOLES library. Can work with only SPOOLES sparse matrix implementation.
 */
class SpoolesSolver : public SparseLinearSystemNM
{
private:
#ifdef __SPOOLES_MODULE
    /// Last mapped LHS matrix
    SparseMtrx *Lhs;
    /// Last mapped matrix version
    SparseMtrx :: SparseMtrxVersionType lhsVersion;
    int msglvl;
    FILE *msgFile;
    int msgFileCloseFlag;

    FrontMtx *frontmtx;
    IV *oldToNewIV, *newToOldIV;
    ETree *frontETree;
    IVL *adjIVL, *symbfacIVL;
    SubMtxManager *mtxmanager;
    Graph *graph;
#endif

public:
    /**
     * Constructor.
     * @param i Solver index.
     * @param d Domain which solver belongs to.
     * @param m Engineering model which solver belongs to.
     */
    SpoolesSolver(int i, Domain *d, EngngModel *m);

    ///Destructor
    virtual ~SpoolesSolver();

    /**
     * Solves the given linear system by LDL^T factorization.
     */
    virtual NM_Status solve(SparseMtrx *A, FloatArray *b, FloatArray *x);

    /// Initializes receiver from given record. Empty implementation.
    virtual IRResultType initializeFrom(InputRecord *ir);

    // identification
    virtual const char *giveClassName() const { return "SpoolesSolver"; }
    virtual classType giveClassID() const { return SpoolesSolverClass; }
    virtual LinSystSolverType giveLinSystSolverType() const { return ST_Spooles; }
};
} // end namespace oofem
#endif // spoolessolver_h
