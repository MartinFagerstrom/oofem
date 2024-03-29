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


#ifndef ilucomprowprecond_h
#define ilucomprowprecond_h

#include "flotarry.h"
#include "intarray.h"
//#include "compcol.h"
#include "dyncomprow.h"
#include "precond.h"

namespace oofem {
/**
 * Implemantation of ILU (Incomplete LU) Preconditioner for compressed row sparse matrices.
 * Fill - up supported.
 */
class CompRow_ILUPreconditioner : public Preconditioner
{
private:
    DynCompRow A;

    double drop_tol;
    int part_fill;

public:
    /// Constructor. Initializes the the receiver (constructs the precontioning matrix M) of given matrix.
    CompRow_ILUPreconditioner(const SparseMtrx &A, InputRecord &attributes);
    /// Constructor. The user should call initializeFrom and init services in this given order to ensure consistency.
    CompRow_ILUPreconditioner() : Preconditioner() { }
    /// Destructor
    ~CompRow_ILUPreconditioner(void) { };

    /**
     * Initializes the receiver (constructs the precontioning matrix M) of given matrix.
     * @param a Sparse matrix to be preconditioned
     */
    virtual void init(const SparseMtrx &a);

    //void initialize (const CompCol &A);
    void initialize(const DynCompRow &A);

    /// Solves the linear system
    virtual void solve(const FloatArray &x, FloatArray &y) const;
    /// Solves transposed system
    virtual void trans_solve(const FloatArray &x, FloatArray &y) const;

    /// returns the preconditioner name
    virtual const char *giveClassName() const { return "ILUT"; }
    /// Initializes receiver from given record. Empty implementation.
    virtual IRResultType initializeFrom(InputRecord *ir);


protected:
    void qsortCol(IntArray &, FloatArray &, int l, int r);
    int  qsortColPartition(IntArray &, FloatArray &, int l, int r);
};
} // end namespace oofem
#endif // ilucomprowprecond_h
