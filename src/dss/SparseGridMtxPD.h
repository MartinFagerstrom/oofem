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

#ifndef _SPARSEGRIDMTXPD_H__
#define _SPARSEGRIDMTXPD_H__

#include "SparseGridMtxLDL.h"

DSS_NAMESPASE_BEGIN

/**
 * @author: Richard Vondracek
 */

class SparseGridMtxPD
{
public:
    // Allocates new space according to bskl and reads old matrix with respect
    // to permutation blockP
    SparseGridMtxPD(SparseMatrixF &sm, long block_size, Ordering *block_order, long fixed_blocks, MathTracer *MT);

    // Allocates new space according to bskl and reads old matrix with respect
    // to permutation blockP
    SparseGridMtxPD(SparseMatrixF &sm, long block_size, Ordering *block_order, Ordering *node_order, long fixed_blocks, MathTracer *MT);

    SparseGridMtxPD(SparseGridMtx *pMatrix, long fixed_blocks);

    virtual ~SparseGridMtxPD();

private:
    SparseGridMtx *m_pMatrix;
    long m_lFixed_blocks;

public:
    // This functions computes LDL' decomposition of first (nblocks-fixed_bn) columns
    // The resulting submatrix will contain the schur complement A22 - A21 * A11^(-1) * A12
    void SchurComplementFactorization();
    void SolveA11(double *x);
    void Sub_A21_A11inv(double *x);
    void Sub_A11inv_A12(double *x);

    void WriteCondensedMatrixA22(double *a, Ordering *mcn, IntArrayList *lncn);

    SparseGridMtx *Matrix();
}; //class SparseGridMtx

DSS_NAMESPASE_END

#endif // _SPARSEGRIDMTX_H__
