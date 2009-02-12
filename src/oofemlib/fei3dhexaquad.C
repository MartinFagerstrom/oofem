/* $Header: /home/cvs/bp/oofem/sm/src/qspace.C,v 1.3.4.1 2004/04/05 15:19:47 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   file fei3dhexaquad.C

//
//   recasted by L. Svoboda
//
#include "fei3dhexaquad.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "node.h"
#include "mathfem.h"
#include "fei2dquadlin.h"


void
FEI3dHexaQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, double time)
{
    double u, v, w;
    answer.resize(20);

    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);
    /*
     *  ***  numbering like T3d  ***            ***  numbering like T3d  ***
     *  ***  numbering of nodes  ***            ***  numbering of surfs  ***
     *				          				         
     *             dzeta                      	                            
     *        1      ^  9         2           	                            
     *         +-----|--+--------+                	   +-----------------+  
     *        /|     |          /|                	  /|                /|  
     *       / |     |         / |                	 / |               / |  
     *    12+  |     o      10+  |                  /  |     1        /  |  
     *     /   |     |       /   |                 /   |             /   |  
     *  4 /  17+  11 |    3 /    +18              /    |        3   /    |  
     *   +--------+--------+     |               +-----------------+     |  
     *   |     |     |     |     |               |     |           |     |  
     *   |     |     +-----|--o------> eta       |  6  |           |  4  |  
     *   |     |    /   13 |     |               |     |           |     |  
     *   |   5 +---/----+--|-----+ 6             |     +-----------|-----+  
     * 20+    /   o        +19  /                |    /   5        |    /   
     *   |   /   /         |   /                 |   /             |   /    
     *   |16+   /          |  +14                |  /        2     |  /     
     *   | /   /           | /                   | /               | /      
     *   |/   L ksi        |/                    |/                |/       
     *   +--------+--------+                     +-----------------+        
     *  8         15        7                                               
     */
    answer.at(1)  = 0.125 * ( 1.0 - u ) * ( 1.0 - v ) * ( 1.0 + w ) * ( -u - v + w - 2.0 ); 
    answer.at(2)  = 0.125 * ( 1.0 - u ) * ( 1.0 + v ) * ( 1.0 + w ) * ( -u + v + w - 2.0 ); 
    answer.at(3)  = 0.125 * ( 1.0 + u ) * ( 1.0 + v ) * ( 1.0 + w ) * ( u + v + w - 2.0 ); 
    answer.at(4)  = 0.125 * ( 1.0 + u ) * ( 1.0 - v ) * ( 1.0 + w ) * ( u - v + w - 2.0 ); 
    answer.at(5)  = 0.125 * ( 1.0 - u ) * ( 1.0 - v ) * ( 1.0 - w ) * ( -u - v - w - 2.0 ); 
    answer.at(6)  = 0.125 * ( 1.0 - u ) * ( 1.0 + v ) * ( 1.0 - w ) * ( -u + v - w - 2.0 ); 
    answer.at(7)  = 0.125 * ( 1.0 + u ) * ( 1.0 + v ) * ( 1.0 - w ) * ( u + v - w - 2.0 ); 
    answer.at(8)  = 0.125 * ( 1.0 + u ) * ( 1.0 - v ) * ( 1.0 - w ) * ( u - v - w - 2.0 ); 
                                                                                          
    answer.at(9)  = 0.25 * ( 1.0 - v * v ) * ( 1.0 - u ) * ( 1.0 + w ); 
    answer.at(10)  = 0.25 * ( 1.0 - u * u ) * ( 1.0 + v ) * ( 1.0 + w ); 
    answer.at(11)  = 0.25 * ( 1.0 - v * v ) * ( 1.0 + u ) * ( 1.0 + w ); 
    answer.at(12)  = 0.25 * ( 1.0 - u * u ) * ( 1.0 - v ) * ( 1.0 + w ); 
                                                                        
    answer.at(13)  = 0.25 * ( 1.0 - v * v ) * ( 1.0 - u ) * ( 1.0 - w ); 
    answer.at(14)  = 0.25 * ( 1.0 - u * u ) * ( 1.0 + v ) * ( 1.0 - w ); 
    answer.at(15)  = 0.25 * ( 1.0 - v * v ) * ( 1.0 + u ) * ( 1.0 - w ); 
    answer.at(16)  = 0.25 * ( 1.0 - u * u ) * ( 1.0 - v ) * ( 1.0 - w ); 
                                                                        
    answer.at(17)  = 0.25 * ( 1.0 - u ) * ( 1.0 - v ) * ( 1.0 - w * w ); 
    answer.at(18)  = 0.25 * ( 1.0 - u ) * ( 1.0 + v ) * ( 1.0 - w * w ); 
    answer.at(19)  = 0.25 * ( 1.0 + u ) * ( 1.0 + v ) * ( 1.0 - w * w ); 
    answer.at(20)  = 0.25 * ( 1.0 + u ) * ( 1.0 - v ) * ( 1.0 - w * w ); 

    return;
}

void
FEI3dHexaQuad :: evaldNdx(FloatMatrix &answer, const FloatArray **coords, const FloatArray &lcoords, double time)
{
    int i;
    FloatMatrix jacobianMatrix(3, 3), inv(3, 3);
    FloatArray dx(20), dy(20), dz(20);
    double u, v, w;

    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);

    this->giveJacobianMatrixAt(jacobianMatrix, coords, lcoords);
    inv.beInverseOf(jacobianMatrix);

    this->giveDerivativeKsi(dx, u, v, w);
    this->giveDerivativeEta(dy, u, v, w);
    this->giveDerivativeDzeta(dz, u, v, w);

    answer.resize(20, 3);

    for ( i = 1; i <= 20; i++ ) {
        answer.at(i, 1) = dx.at(i) * inv.at(1, 1) + dy.at(i) * inv.at(1, 2) + dz.at(i) * inv.at(1, 3);
        answer.at(i, 2) = dx.at(i) * inv.at(2, 1) + dy.at(i) * inv.at(2, 2) + dz.at(i) * inv.at(2, 3);
        answer.at(i, 3) = dx.at(i) * inv.at(3, 1) + dy.at(i) * inv.at(3, 2) + dz.at(i) * inv.at(3, 3);
    }
}

void
FEI3dHexaQuad :: local2global(FloatArray &answer, const FloatArray **coords, const FloatArray &lcoords, double time)
{
    int i;
    FloatArray n(20);

    this->evalN(n, lcoords, time);   // inspirace

    answer.resize(3);
    answer.zero();

    for ( i = 1; i <= 20; i++ ) {
        answer.at(1) += n.at(i) * coords [ i - 1 ]->at(1);
        answer.at(2) += n.at(i) * coords [ i - 1 ]->at(2);
        answer.at(3) += n.at(i) * coords [ i - 1 ]->at(3);
    }
}

// #define POINT_TOL 1.e-3
//
int
FEI3dHexaQuad :: global2local(FloatArray &answer, const FloatArray **coords, const FloatArray &lcoords, double time)
{ 
  OOFEM_ERROR("FEI3dHexaQuad :: global2local not implemented");
  return 1; 
}

double
FEI3dHexaQuad :: giveTransformationJacobian(const FloatArray **coords, const FloatArray &lcoords, double time)
{
    FloatMatrix jacobianMatrix(3, 3);

    this->giveJacobianMatrixAt(jacobianMatrix, coords, lcoords);
    return jacobianMatrix.giveDeterminant();
}

void FEI3dHexaQuad :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords, double time)
{ OOFEM_ERROR("FEI3dHexaQuad :: edgeEvalN not implemented"); }
void FEI3dHexaQuad :: edgeEvaldNdx(FloatMatrix &answer, int iedge, const FloatArray **coords, const FloatArray &lcoords, double time)
{ OOFEM_ERROR("FEI3dHexaQuad :: edgeEvaldNdx not implemented"); }
void FEI3dHexaQuad :: edgeLocal2global(FloatArray &answer, int iedge, const FloatArray **coords, const FloatArray &lcoords, double time)
{ OOFEM_ERROR("FEI3dHexaQuad :: edgeLocal2global not implemented"); }
double FEI3dHexaQuad :: edgeGiveTransformationJacobian(int iedge, const FloatArray **coords, const FloatArray &lcoords, double time)
{ OOFEM_ERROR("FEI3dHexaQuad :: edgeGiveTransformationJacobian not implemented");
  return 0.0; }

void
FEI3dHexaQuad :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{ OOFEM_ERROR("FEI3dHexaQuad :: computeLocalEdgeMapping not implemented"); }

void
FEI3dHexaQuad :: surfaceEvalN(FloatArray &answer, const FloatArray &lcoords, double time)
{
    double ksi, eta;
    answer.resize(8);

    ksi = lcoords.at(1);
    eta = lcoords.at(2);

    answer.at(1) = ( 1. + ksi ) * ( 1. + eta ) * 0.25 * ( ksi + eta - 1. );
    answer.at(2) = ( 1. - ksi ) * ( 1. + eta ) * 0.25 * ( -ksi + eta - 1. );
    answer.at(3) = ( 1. - ksi ) * ( 1. - eta ) * 0.25 * ( -ksi - eta - 1. );
    answer.at(4) = ( 1. + ksi ) * ( 1. - eta ) * 0.25 * ( ksi - eta - 1. );
    answer.at(5) = 0.5 * ( 1. - ksi * ksi ) * ( 1. + eta );
    answer.at(6) = 0.5 * ( 1. - ksi ) * ( 1. - eta * eta );
    answer.at(7) = 0.5 * ( 1. - ksi * ksi ) * ( 1. - eta );
    answer.at(8) = 0.5 * ( 1. + ksi ) * ( 1. - eta * eta );

    return;
}


void
FEI3dHexaQuad :: surfaceLocal2global(FloatArray &answer, int isurf,
                                     const FloatArray **coords, const FloatArray &lcoords, double time)
{
    IntArray nodes(8);
    FloatArray n;

    computeLocalSurfaceMapping(nodes, isurf);

    this->surfaceEvalN(n, lcoords, time); // inspirace

    answer.resize(3);
    answer.zero();

    for ( int i = 1; i <= 8; i++ ) {
        answer.at(1) += n.at(i) * coords [ nodes.at(i) - 1 ]->at(1);
        answer.at(2) += n.at(i) * coords [ nodes.at(i) - 1 ]->at(2);
        answer.at(3) += n.at(i) * coords [ nodes.at(i) - 1 ]->at(3);
    }
}



double
FEI3dHexaQuad :: surfaceGiveTransformationJacobian(int isurf, const FloatArray **coords, const FloatArray &lcoords,
                                                   double time)
{
    // only plane surface is supported !!!

    // FEI2dQuadLin cannot be used without modifying coordinates of surface nodes
    // therefore local calculation is done without using FEI2dQuadLin :: giveJacobianMatrixAt

    int n1, n2, n3, n4, n;
    IntArray snodes;
    FloatMatrix jacobianMatrix(2, 2);
    FloatArray sn [ 4 ];
    FloatArray const *nod1, *nod, *nod2;
    double length, k;
    int i;

    this->computeLocalSurfaceMapping(snodes, isurf);

    // check whether all surface nodes are in plane
    n1 = snodes.at(1) - 1;
    n2 = snodes.at(2) - 1;
    n3 = snodes.at(3) - 1;
    n4 = snodes.at(4) - 1;

    // get the normal using nodes 1 2 3
    FloatArray a(3), b(3), c(3);
    for ( i = 1; i <= 3; i++ ) {
        b.at(i) = coords [ n2 ]->at(i) - coords [ n1 ]->at(i);
        a.at(i) = coords [ n3 ]->at(i) - coords [ n1 ]->at(i);
    }

    n = n4;
    c.beVectorProductOf(b, a);
    length = sqrt( dotProduct(c, c, 3) );

    if ( length < 1.0e-10 ) {
        // try nodes 1 3 4
        for ( i = 1; i <= 3; i++ ) {
            b.at(i) = coords [ n4 ]->at(i) - coords [ n1 ]->at(i);
        }

        n = n2;
        c.beVectorProductOf(a, b);
        length = sqrt( dotProduct(c, c, 3) );

        if ( length < 1.0e-10 ) {
            OOFEM_ERROR("FEI3dHexaQuad :: surfaceGiveTransformationJacobian: degenerated surface");
        }
    }

    // c is normed vector
    c.times(1.0 / length);
    for ( i = 1; i <= 3; i++ ) {
        b.at(i) = coords [ n ]->at(i) - coords [ n1 ]->at(i);
    }

    // check distance of the 4th node n
    if ( fabs( dotProduct(b, c, 3) ) > 1.0e-6 ) {
        OOFEM_ERROR("FEI3dHexaQuad :: surfaceGiveTransformationJacobian: not planar surface");
    }

    // check whether edges are direct
    for ( i = 1; i <= 4; i++ ) {
        nod1 = coords [ snodes.at(i) - 1 ];
        nod2 = coords [ snodes.at(i % 4 + 1) - 1 ];
        nod  = coords [ snodes.at(i + 4) - 1 ];

        k =           ( nod2->at(1) - nod1->at(1) ) / ( nod->at(1) - nod1->at(1) );
        if ( fabs( k - ( nod2->at(2) - nod1->at(2) ) / ( nod->at(2) - nod1->at(2) ) ) > 1.0e-6 ||
             fabs( k - ( nod2->at(3) - nod1->at(3) ) / ( nod->at(3) - nod1->at(3) ) ) > 1.0e-6    ) {
            OOFEM_ERROR("surfaceGiveTransformationJacobian: not direct edge");
        }
    }

    // map nodes to the surface (x,y) plane
    // get x and y unit vectors a and b (a is already computed as n3-n1)
    length = sqrt( dotProduct(a, a, 3) );
    a.times(1.0 / length);
    b.beVectorProductOf(c, a);

    for ( i = 0; i < 4; i++ ) {
        sn [ i ].resize(2);
    }

    sn [ 0 ].at(1) = 0.0;
    sn [ 0 ].at(2) = 0.0;

    sn [ 2 ].at(1) = length;
    sn [ 2 ].at(2) = 0.0;

    for ( i = 1; i <= 3; i++ ) {
        c.at(i) = coords [ n2 ]->at(i) - coords [ n1 ]->at(i);
    }

    sn [ 1 ].at(1) = dotProduct(c, a, 3);
    sn [ 1 ].at(2) = dotProduct(c, b, 3);

    for ( i = 1; i <= 3; i++ ) {
        c.at(i) = coords [ n4 ]->at(i) - coords [ n1 ]->at(i);
    }

    sn [ 3 ].at(1) = dotProduct(c, a, 3);
    sn [ 3 ].at(2) = dotProduct(c, b, 3);

    double ksi, eta, x, y;
    FloatArray nx(4), ny(4);

    jacobianMatrix.resize(2, 2);
    jacobianMatrix.zero();

    ksi = lcoords.at(1);
    eta = lcoords.at(2);

    //FEI2dQuadLin interpolation(1, 2);
    //interpolation.giveDerivativeKsi(nx, eta);
    //interpolation.giveDerivativeEta(ny, ksi);

    nx.at(1) =  0.25 * ( 1. + eta );
    nx.at(2) = -0.25 * ( 1. + eta );
    nx.at(3) = -0.25 * ( 1. - eta );
    nx.at(4) =  0.25 * ( 1. - eta );

    ny.at(1) =  0.25 * ( 1. + ksi );
    ny.at(2) =  0.25 * ( 1. - ksi );
    ny.at(3) = -0.25 * ( 1. - ksi );
    ny.at(4) = -0.25 * ( 1. + ksi );


    for ( i = 1; i <= 4; i++ ) {
        x = sn [ i - 1 ].at(1);
        y = sn [ i - 1 ].at(2);

        jacobianMatrix.at(1, 1) += nx.at(i) * x;
        jacobianMatrix.at(1, 2) += nx.at(i) * y;
        jacobianMatrix.at(2, 1) += ny.at(i) * x;
        jacobianMatrix.at(2, 2) += ny.at(i) * y;
    }

    return jacobianMatrix.giveDeterminant();
}


void
FEI3dHexaQuad :: computeLocalSurfaceMapping(IntArray &nodes, int isurf)
{
    nodes.resize(8);

    // the actual numbering  has a positive normal pointing outwards from the element  - (LSpace compatible)
    //

    if      ( isurf == 1 ) { // surface 1 - nodes   3 2 1 4  10  9 12 11
        nodes.at(1) =  3;
        nodes.at(2) =  2;
        nodes.at(3) =  1;
        nodes.at(4) =  4;
        nodes.at(5) = 10;
        nodes.at(6) =  9;
        nodes.at(7) = 12;
        nodes.at(8) = 11;
    } else if ( isurf == 2 ) { // surface 2 - nodes   7 8 5 6  15 16 13 14
        nodes.at(1) =  7;
        nodes.at(2) =  8;
        nodes.at(3) =  5;
        nodes.at(4) =  6;
        nodes.at(5) = 15;
        nodes.at(6) = 16;
        nodes.at(7) = 13;
        nodes.at(8) = 14;
    } else if ( isurf == 3 ) { // surface 3 - nodes   2 6 5 1  18 13 17  9
        nodes.at(1) =  2;
        nodes.at(2) =  6;
        nodes.at(3) =  5;
        nodes.at(4) =  1;
        nodes.at(5) = 18;
        nodes.at(6) = 13;
        nodes.at(7) = 17;
        nodes.at(8) =  9;
    } else if ( isurf == 4 )     { // surface 4 - nodes   3 7 6 2  19 14 18 10
        nodes.at(1) =  3;
        nodes.at(2) =  7;
        nodes.at(3) =  6;
        nodes.at(4) =  2;
        nodes.at(5) = 19;
        nodes.at(6) = 14;
        nodes.at(7) = 18;
        nodes.at(8) = 10;
    } else if ( isurf == 5 )     { // surface 5 - nodes   3 4 8 7  11 20 15 19
        nodes.at(1) =  3;
        nodes.at(2) =  4;
        nodes.at(3) =  8;
        nodes.at(4) =  7;
        nodes.at(5) = 11;
        nodes.at(6) = 20;
        nodes.at(7) = 15;
        nodes.at(8) = 19;
    } else if ( isurf == 6 ) { // surface 6 - nodes   4 1 5 8  12 17 16 20
        nodes.at(1) =  4;
        nodes.at(2) =  1;
        nodes.at(3) =  5;
        nodes.at(4) =  8;
        nodes.at(5) = 12;
        nodes.at(6) = 17;
        nodes.at(7) = 16;
        nodes.at(8) = 20;
    } else   {
        OOFEM_ERROR2("FEI3dHexaQuad :: computeLocalSurfaceMapping: wrong surface number (%d)", isurf);
    }

    /*
    // this commented numbering is symmetrical with respect to local coordinate axes

    if      ( isurf == 1 ) { // surface 1 - nodes   3 2 1 4  10  9 12 11
        nodes.at(1) =  3;
        nodes.at(2) =  2;
        nodes.at(3) =  1;
        nodes.at(4) =  4;
        nodes.at(5) = 10;
        nodes.at(6) =  9;
        nodes.at(7) = 12;
        nodes.at(8) = 11;
    } else if (iSurf == 2) { // surface 2 - nodes   7 6 5 8  14 13 16 15
      nodes.at(1) =  7;  nodes.at(2) =  6;  nodes.at(3) =  5;  nodes.at(4) =  8;
       nodes.at(5) = 14;  nodes.at(6) = 13;  nodes.at(7) = 16;  nodes.at(8) = 15;
    } else if (iSurf == 3) { // surface 3 - nodes   2 1 5 6   9 17 13 18
      nodes.at(1) =  2;  nodes.at(2) =  1;  nodes.at(3) =  5;  nodes.at(4) =  6;
      nodes.at(5) =  9;  nodes.at(6) = 17;  nodes.at(7) = 13;  nodes.at(8) = 18;
    } else if ( isurf == 4 )     { // surface 4 - nodes   3 7 6 2  19 14 18 10
        nodes.at(1) =  3;
        nodes.at(2) =  7;
        nodes.at(3) =  6;
        nodes.at(4) =  2;
        nodes.at(5) = 19;
        nodes.at(6) = 14;
        nodes.at(7) = 18;
        nodes.at(8) = 10;
    } else if ( isurf == 5 )     { // surface 5 - nodes   3 4 8 7  11 20 15 19
        nodes.at(1) =  3;
        nodes.at(2) =  4;
        nodes.at(3) =  8;
        nodes.at(4) =  7;
        nodes.at(5) = 11;
        nodes.at(6) = 20;
        nodes.at(7) = 15;
        nodes.at(8) = 19;
    } else if (iSurf == 6) { // surface 6 - nodes   4 8 5 1  20 16 17 12
      nodes.at(1) =  4;  nodes.at(2) =  8;  nodes.at(3) =  5;  nodes.at(4) =  1;
      nodes.at(5) = 20;  nodes.at(6) = 16;  nodes.at(7) = 17;  nodes.at(8) = 12;
    } else   {
      OOFEM_ERROR2("FEI3dHexaQuad :: computeLocalSurfaceMapping: wrong surface number (%d)", isurf);
    }
    */
}


void
FEI3dHexaQuad :: computeGlobalSurfaceMapping(IntArray &surfNodes, IntArray &elemNodes, int iSurf)
{
    IntArray nodes;
    surfNodes.resize(8);

    computeLocalSurfaceMapping(nodes, iSurf);

    for ( int i = 1; i <= 8; i++ ) {
        surfNodes.at(i) = elemNodes.at( nodes.at(i) );
    }
}


void
FEI3dHexaQuad :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray **coords, const FloatArray &lcoords)
// Returns the jacobian matrix  J (x,y,z)/(ksi,eta,dzeta)  of the receiver.
// Computes it if it does not exist yet.
{
    int i;
    double x, y, z, u, v, w;
    FloatArray dx(20), dy(20), dz(20);

    jacobianMatrix.resize(3, 3);
    jacobianMatrix.zero();

    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);

    this->giveDerivativeKsi(dx, u, v, w);
    this->giveDerivativeEta(dy, u, v, w);
    this->giveDerivativeDzeta(dz, u, v, w);

    for ( i = 1; i <= 20; i++ ) {
        x = coords [ i - 1 ]->at(1);
        y = coords [ i - 1 ]->at(2);
        z = coords [ i - 1 ]->at(3);

        jacobianMatrix.at(1, 1) += dx.at(i) * x;
        jacobianMatrix.at(1, 2) += dx.at(i) * y;
        jacobianMatrix.at(1, 3) += dx.at(i) * z;
        jacobianMatrix.at(2, 1) += dy.at(i) * x;
        jacobianMatrix.at(2, 2) += dy.at(i) * y;
        jacobianMatrix.at(2, 3) += dy.at(i) * z;
        jacobianMatrix.at(3, 1) += dz.at(i) * x;
        jacobianMatrix.at(3, 2) += dz.at(i) * y;
        jacobianMatrix.at(3, 3) += dz.at(i) * z;
    }
}


void
FEI3dHexaQuad :: giveDerivativeKsi(FloatArray &dx, double u, double v, double w)
{
    dx.at(1)  = 0.125 * ( 1.0 - v ) * ( 1.0 + w ) * ( 2.0 * u + v - w + 1.0 );
    dx.at(2)  = 0.125 * ( 1.0 + v ) * ( 1.0 + w ) * ( 2.0 * u - v - w + 1.0 );
    dx.at(3)  = 0.125 * ( 1.0 + v ) * ( 1.0 + w ) * ( 2.0 * u + v + w - 1.0 );
    dx.at(4)  = 0.125 * ( 1.0 - v ) * ( 1.0 + w ) * ( 2.0 * u - v + w - 1.0 );
    dx.at(5)  = 0.125 * ( 1.0 - v ) * ( 1.0 - w ) * ( 2.0 * u + v + w + 1.0 );
    dx.at(6)  = 0.125 * ( 1.0 + v ) * ( 1.0 - w ) * ( 2.0 * u - v + w + 1.0 );
    dx.at(7)  = 0.125 * ( 1.0 + v ) * ( 1.0 - w ) * ( 2.0 * u + v - w - 1.0 );
    dx.at(8)  = 0.125 * ( 1.0 - v ) * ( 1.0 - w ) * ( 2.0 * u - v - w - 1.0 );

    dx.at(9)  = -0.25 * ( 1.0 - v * v ) * ( 1.0 + w );
    dx.at(10)  =  -0.5 * u * ( 1.0 + v ) * ( 1.0 + w );
    dx.at(11)  =  0.25 * ( 1.0 - v * v ) * ( 1.0 + w );
    dx.at(12)  =  -0.5 * u * ( 1.0 - v ) * ( 1.0 + w );

    dx.at(13)  = -0.25 * ( 1.0 - v * v ) * ( 1.0 - w );
    dx.at(14)  =  -0.5 * u * ( 1.0 + v ) * ( 1.0 - w );
    dx.at(15)  =  0.25 * ( 1.0 - v * v ) * ( 1.0 - w );
    dx.at(16)  =  -0.5 * u * ( 1.0 - v ) * ( 1.0 - w );

    dx.at(17)  = -0.25 * ( 1.0 - v ) * ( 1.0 - w * w );
    dx.at(18)  = -0.25 * ( 1.0 + v ) * ( 1.0 - w * w );
    dx.at(19)  =  0.25 * ( 1.0 + v ) * ( 1.0 - w * w );
    dx.at(20)  =  0.25 * ( 1.0 - v ) * ( 1.0 - w * w );
}

void
FEI3dHexaQuad :: giveDerivativeEta(FloatArray &dy, double u, double v, double w)
{
    dy.at(1)  = 0.125 * ( 1.0 - u ) * ( 1.0 + w ) * ( 2.0 * v + u - w + 1.0 );
    dy.at(2)  = 0.125 * ( 1.0 - u ) * ( 1.0 + w ) * ( 2.0 * v - u + w - 1.0 );
    dy.at(3)  = 0.125 * ( 1.0 + u ) * ( 1.0 + w ) * ( 2.0 * v + u + w - 1.0 );
    dy.at(4)  = 0.125 * ( 1.0 + u ) * ( 1.0 + w ) * ( 2.0 * v - u - w + 1.0 );
    dy.at(5)  = 0.125 * ( 1.0 - u ) * ( 1.0 - w ) * ( 2.0 * v + u + w + 1.0 );
    dy.at(6)  = 0.125 * ( 1.0 - u ) * ( 1.0 - w ) * ( 2.0 * v - u - w - 1.0 );
    dy.at(7)  = 0.125 * ( 1.0 + u ) * ( 1.0 - w ) * ( 2.0 * v + u - w - 1.0 );
    dy.at(8)  = 0.125 * ( 1.0 + u ) * ( 1.0 - w ) * ( 2.0 * v - u + w + 1.0 );

    dy.at(9)  =  -0.5 * v * ( 1.0 - u ) * ( 1.0 + w );
    dy.at(10)  =  0.25 * ( 1.0 - u * u ) * ( 1.0 + w );
    dy.at(11)  =  -0.5 * v * ( 1.0 + u ) * ( 1.0 + w );
    dy.at(12)  = -0.25 * ( 1.0 - u * u ) * ( 1.0 + w );

    dy.at(13)  =  -0.5 * v * ( 1.0 - u ) * ( 1.0 - w );
    dy.at(14)  =  0.25 * ( 1.0 - u * u ) * ( 1.0 - w );
    dy.at(15)  =  -0.5 * v * ( 1.0 + u ) * ( 1.0 - w );
    dy.at(16)  = -0.25 * ( 1.0 - u * u ) * ( 1.0 - w );

    dy.at(17)  = -0.25 * ( 1.0 - u ) * ( 1.0 - w * w );
    dy.at(18)  =  0.25 * ( 1.0 - u ) * ( 1.0 - w * w );
    dy.at(19)  =  0.25 * ( 1.0 + u ) * ( 1.0 - w * w );
    dy.at(20)  = -0.25 * ( 1.0 + u ) * ( 1.0 - w * w );
}

void
FEI3dHexaQuad :: giveDerivativeDzeta(FloatArray &dz, double u, double v, double w)
{
    dz.at(1)  = 0.125 * ( 1.0 - u ) * ( 1.0 - v ) * ( 2.0 * w - u - v - 1.0 );
    dz.at(2)  = 0.125 * ( 1.0 - u ) * ( 1.0 + v ) * ( 2.0 * w - u + v - 1.0 );
    dz.at(3)  = 0.125 * ( 1.0 + u ) * ( 1.0 + v ) * ( 2.0 * w + u + v - 1.0 );
    dz.at(4)  = 0.125 * ( 1.0 + u ) * ( 1.0 - v ) * ( 2.0 * w + u - v - 1.0 );
    dz.at(5)  = 0.125 * ( 1.0 - u ) * ( 1.0 - v ) * ( 2.0 * w + u + v + 1.0 );
    dz.at(6)  = 0.125 * ( 1.0 - u ) * ( 1.0 + v ) * ( 2.0 * w + u - v + 1.0 );
    dz.at(7)  = 0.125 * ( 1.0 + u ) * ( 1.0 + v ) * ( 2.0 * w - u - v + 1.0 );
    dz.at(8)  = 0.125 * ( 1.0 + u ) * ( 1.0 - v ) * ( 2.0 * w - u + v + 1.0 );

    dz.at(9)  =  0.25 * ( 1.0 - v * v ) * ( 1.0 - u );
    dz.at(10)  =  0.25 * ( 1.0 - u * u ) * ( 1.0 + v );
    dz.at(11)  =  0.25 * ( 1.0 - v * v ) * ( 1.0 + u );
    dz.at(12)  =  0.25 * ( 1.0 - u * u ) * ( 1.0 - v );

    dz.at(13)  = -0.25 * ( 1.0 - v * v ) * ( 1.0 - u );
    dz.at(14)  = -0.25 * ( 1.0 - u * u ) * ( 1.0 + v );
    dz.at(15)  = -0.25 * ( 1.0 - v * v ) * ( 1.0 + u );
    dz.at(16)  = -0.25 * ( 1.0 - u * u ) * ( 1.0 - v );

    dz.at(17)  =  -0.5 * ( 1.0 - u ) * ( 1.0 - v ) * w;
    dz.at(18)  =  -0.5 * ( 1.0 - u ) * ( 1.0 + v ) * w;
    dz.at(19)  =  -0.5 * ( 1.0 + u ) * ( 1.0 + v ) * w;
    dz.at(20)  =  -0.5 * ( 1.0 + u ) * ( 1.0 - v ) * w;
}