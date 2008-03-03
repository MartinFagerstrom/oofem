/* $Header: /home/cvs/bp/oofem/tm/src/quadaxisym1_ht.C,v 1.1 2003/04/14 16:01:40 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2002   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

#include "quadaxisym1_ht.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "structuralms.h"
#include "load.h"
#ifndef __MAKEDEPEND
#include <math.h>
#include <stdio.h>
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "conTable.h"
#endif


QuadAxisym1_ht :: QuadAxisym1_ht (int n, Domain* aDomain, ElementMode em)
: Quad1_ht (n,aDomain, em)
   // Constructor.
{}

QuadAxisym1_ht :: ~QuadAxisym1_ht ()
// Destructor
{}

double
QuadAxisym1_ht :: computeVolumeAround (GaussPoint* aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
   double       determinant,weight,volume;
  FloatMatrix jacobianMtrx;
  this->computeJacobianMatrix(jacobianMtrx,aGaussPoint);
   determinant = fabs (jacobianMtrx.giveDeterminant());
   weight      = aGaussPoint -> giveWeight();
   volume      = determinant * weight * this->computeRadiusAt(aGaussPoint) ;
  
   return volume;
}

double
QuadAxisym1_ht :: computeEdgeVolumeAround(GaussPoint* gp, int iEdge)
{ 
 double dx,dy, length, radius;
  Node   *nodeA,*nodeB ;
  int aNode=0, bNode=0;
 FloatMatrix n;
 
  if (iEdge == 1) { // edge between nodes 1 2
  aNode = 1;
  bNode = 2;
  } else if (iEdge == 2) { // edge between nodes 2 3
  aNode = 2;
  bNode = 3;
  } else if (iEdge == 3) { // edge between nodes 3 4
  aNode = 3;
  bNode = 4;
  } else if (iEdge == 4) { // edge between nodes 4 1
  aNode = 4;
  bNode = 1;
  } else {
  _error ("computeEdgeVolumeAround: wrong egde number");
  }
 
  nodeA   = this->giveNode(aNode) ;
  nodeB   = this->giveNode(bNode) ;
  
  dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1) ;
  dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2) ;
  length = sqrt(dx*dx + dy*dy);
 this->computeEgdeNMatrixAt (n, gp);
 radius = n.at(1, 1)*nodeA->giveCoordinate(1)+n.at(1, 2)*nodeB->giveCoordinate(1);
 
  return 0.5 * length * gp -> giveWeight() * radius ;
 
}

double
QuadAxisym1_ht :: computeRadiusAt (GaussPoint* gp)
{
 double r = 0.0;
 int i;
 FloatMatrix n;

 this->computeNmatrixAt (n, gp->giveCoordinates());
 for (i=1; i<=4; i++) r+= n.at(1,i)*this->giveNode(i)->giveCoordinate(1);
 return r;
}