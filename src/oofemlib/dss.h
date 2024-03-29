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

#ifndef dss_h
#define dss_h

#ifdef __DSS_MODULE

#include "sparsemtrx.h"
#include "intarray.h"

/* DSS module lives in global namespace, not in oofem namespace */
class DSSolver;
struct SparseMatrixF;

namespace oofem {
/**
 * Interface to Direct Sparse Solver written by R.Vonracek.
 * This class represent the sparse matrix interface to DSS library. It allows to build internal structure,
 * assemble the DSS sparse matrix, and to factorize and back substitution operations.
 */
class DSSMatrix : public SparseMtrx
{
public:
    /// Possible storage schemes and factorization types
    enum dssType { sym_LDL, sym_LL, unsym_LU };

protected:
    /// Pointer to SparseMatrixF class rep
    SparseMatrixF *_sm;
    /// pointer to DSSolver class (representation of Direct Sparse Solver in DSS lib)
    DSSolver *_dss;
    /// Flag indicating whether factorized.
    bool isFactorized;
    /// type of storage & factorization
    dssType _type;

    /// implements 0-based access
    double operator()(int i, int j) const;
    /// implements 0-based access
    double &operator()(int i, int j);

public:
    /**
     * Constructor. Before any operation an internal profile must be built.
     * @param t Storage type
     * @param n Size of row and columns of square matrix
     * @see buildInternalStructure
     */
    DSSMatrix(dssType t, int n);
    /**
     * Constructor. Before any operation an internal profile must be built.
     * @param t Storage type
     * @see buildInternalStructure
     */
    DSSMatrix(dssType t);
    /// Copy constructor
    DSSMatrix(const DSSMatrix &S);
    /// Destructor
    ~DSSMatrix();

    // Overloaded methods
    SparseMtrx *GiveCopy() const;
    void times(const FloatArray &x, FloatArray &answer) const;
    virtual void times(double x);
    int buildInternalStructure(EngngModel *, int, EquationID, const UnknownNumberingScheme & s);
    int assemble(const IntArray &loc, const FloatMatrix &mat);
    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);
    bool canBeFactorized() const { return true; }
    virtual SparseMtrx *factorized();
    void solve(FloatArray *b, FloatArray *x);
    void zero();
    double &at(int i, int j);
    double at(int i, int j) const;
    SparseMtrxType giveType() const { return SMT_SymCompCol; }
    bool isAsymmetric() const { return false; }
};
} // end namespace oofem
#endif // ifdef __DSS_MODULE
#endif // dss_h

