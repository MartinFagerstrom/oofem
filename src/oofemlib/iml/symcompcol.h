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

// Inspired by SPARSELib++
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             ********   ***                                 SparseLib++    */
/*          *******  **  ***       ***      ***               v. 1.5c        */
/*           *****      ***     ******** ********                            */
/*            *****    ***     ******** ********              R. Pozo        */
/*       **  *******  ***   **   ***      ***                 K. Remington   */
/*        ********   ********                                 A. Lumsdaine   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                     SparseLib++ : Sparse Matrix Library                   */
/*                                                                           */
/*               National Institute of Standards and Technology              */
/*                        University of Notre Dame                           */
/*              Authors: R. Pozo, K. Remington, A. Lumsdaine                 */
/*                                                                           */
/*                                 NOTICE                                    */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that the above notice appear in all copies and supporting        */
/* documentation.                                                            */
/*                                                                           */
/* Neither the Institutions (National Institute of Standards and Technology, */
/* University of Notre Dame) nor the Authors make any representations about  */
/* the suitability of this software for any purpose.  This software is       */
/* provided ``as is'' without expressed or implied warranty.                 */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*        Compressed column sparse matrix  (O-based, Fortran)            */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#ifndef symcompcol_h
#define symcompcol_h

#include "compcol.h"

namespace oofem {
/**
 * Implementation of symmetric sparse matrix stored using compressed column/row storage.
 * Only the lower part is stored.
 */
class SymCompCol : public CompCol
{
public:
    /**
     * Constructor.
     * Before any operation an internal profile must be built.
     * @param n Size of matrix
     * @see buildInternalStructure
     */
    SymCompCol(int n);
    /**
     * Constructor.
     * Before any operation an internal profile must be built.
     * @see buildInternalStructure
     */
    SymCompCol();
    /// Copy constructor
    SymCompCol(const SymCompCol &S);
    /// Destructor
    ~SymCompCol() { }

    // Overloaded methods
    SparseMtrx *GiveCopy() const;
    void times(const FloatArray &x, FloatArray &answer) const;
    void timesT(const FloatArray &x, FloatArray &answer) const { this->times(x,answer); }
    virtual void times(double x);
    int buildInternalStructure(EngngModel *, int, EquationID, const UnknownNumberingScheme &);
    int assemble(const IntArray &loc, const FloatMatrix &mat);
    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);
    bool canBeFactorized() const { return false; }
    virtual void zero();
    virtual double &at(int i, int j);
    virtual double at(int i, int j) const;
    SparseMtrxType giveType() const { return SMT_SymCompCol; }
    bool isAntisymmetric() const { return false; }


    const double &val(int i) const { return val_(i); }
    const int &row_ind(int i) const { return rowind_(i); }
    const int &col_ptr(int i) const { return colptr_(i); }
    int dim(int i) const { return dim_ [ i ]; };

protected:
    /*******************************/
    /*  Access and info functions  */
    /*******************************/

    double &val(int i) { return val_(i); }
    int &row_ind(int i) { return rowind_(i); }
    int &col_ptr(int i) { return colptr_(i); }

    int size(int i) const { return dim_ [ i ]; };
    int NumNonzeros() const { return nz_; };
    int base() const { return base_; }

    /***********************************/
    /*  General access function (slow) */
    /***********************************/
    /// implements 0-based access
    double operator()(int i, int j) const;
    /// implements 0-based access
    double &operator()(int i, int j);
};
} // end namespace oofem
#endif // symcompcol_h
