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

// Recasted by Ladislav Svoboda

#ifndef fei3dhexaquad_h
#define fei3dhexaquad_h

#include "feinterpol3d.h"

namespace oofem {
/**
 * Class representing implementation of quadratic hexahedra interpolation class.
 *  ***  numbering like T3d  ***            ***  numbering like T3d  ***
 *  ***  numbering of nodes  ***            ***  numbering of surfs  ***
 *
 *              zeta
 *        1      ^  9         2
 *         +-----|--+--------+                     +-----------------+
 *        /|     |          /|                    /|                /|
 *       / |     |         / |                   / |               / |
 *    12+  |     o      10+  |                  /  |     1        /  |
 *     /   |     |       /   |                 /   |             /   |
 *  4 /  17+  11 |    3 /    +18              /    |       (3)  /    |
 *   +--------+--------+     |               +-----------------+     |
 *   |     |     |     |     |               |     |           |     |
 *   |     |     +-----|--o------> eta       | (6) |           |  4  |
 *   |     |    /   13 |     |               |     |           |     |
 *   |   5 +---/----+--|-----+ 6             |     +-----------|-----+
 * 20+    /   o        +19  /                |    /   5        |    /
 *   |   /   /         |   /                 |   /             |   /
 *   |16+   /          |  +14                |  /       (2)    |  /
 *   | /   /           | /                   | /               | /
 *   |/   L ksi        |/                    |/                |/
 *   +--------+--------+                     +-----------------+
 *  8         15        7
 */
class FEI3dHexaQuad : public FEInterpolation3d
{
public:
    FEI3dHexaQuad() : FEInterpolation3d(2) { }

    virtual double giveCharacteristicLength(const FEICellGeometry &cellgeo) const;

    // Bulk
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo);

    // Edge
    virtual void edgeEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void edgeEvaldNdx(FloatMatrix &answer, int iedge,
                              const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void edgeLocal2global(FloatArray &answer, int iedge,
                                  const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords,
                                                  const FEICellGeometry &cellgeo);
    virtual void computeLocalEdgeMapping(IntArray &edgeNodes, int iedge);

    // Surface
    virtual void surfaceEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    //virtual void surfaceEvaldNdx (FloatMatrix&answer, int isurf, const FloatArray& lcoords, const FEICellGeometry& cellgeo);
    virtual void surfaceLocal2global(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void computeLocalSurfaceMapping(IntArray &nodes, int iSurf);
    void computeGlobalSurfaceMapping(IntArray &edgeNodes, IntArray &elemNodes, int iedge);

    virtual void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo);

protected:
    double edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo);
    void giveLocalDerivative(FloatMatrix &dN, const FloatArray &lcoords);
};
} // end namespace oofem
#endif
