/* $Header: /home/cvs/bp/oofem/sm/src/quad1planestrain.C,v 1.1.4.1 2004/04/05 15:19:47 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



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

#include "quad1planestrain.h"
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
#ifndef __MAKEDEPEND
#include <math.h>
#include <stdio.h>
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "oofegutils.h"
#include "conTable.h"
#endif


Quad1PlaneStrain :: Quad1PlaneStrain (int n, Domain* aDomain)
             : StructuralElement (n,aDomain), ZZNodalRecoveryModelInterface(), SPRNodalRecoveryModelInterface(),
        SpatialLocalizerInterface()
   , DirectErrorIndicatorRCInterface(), EIPrimaryUnknownMapperInterface(),
   HuertaErrorEstimatorInterface(), HuertaRemeshingCriteriaInterface()
// Constructor.
{
   numberOfDofMans  = 4 ;
   jacobianMatrix = NULL ;
}

Quad1PlaneStrain :: ~Quad1PlaneStrain ()
// Destructor
{
 delete jacobianMatrix;
}

void
Quad1PlaneStrain :: giveDerivativeKsi (FloatArray& answer, double eta)
{
 answer.resize(4);

 answer.at(1) =  0.25*(1.+eta);
 answer.at(2) = -0.25*(1.+eta);
 answer.at(3) = -0.25*(1.-eta);
 answer.at(4) =  0.25*(1.-eta);

}

void
Quad1PlaneStrain :: giveDerivativeEta (FloatArray& answer, double ksi)
{
 answer.resize(4);
 
 answer.at(1) =  0.25*(1.+ksi);
 answer.at(2) =  0.25*(1.-ksi);
 answer.at(3) = -0.25*(1.-ksi);
 answer.at(4) = -0.25*(1.+ksi);

}

void
Quad1PlaneStrain :: computeBmatrixAt (GaussPoint *aGaussPoint, FloatMatrix& answer, int li, int ui)
//
// Returns the [4x8] strain-displacement matrix {B} of the receiver,
// evaluated at aGaussPoint.
// (epsilon_x,epsilon_y,epsilon_z,gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3,u4,v4)
{
 int i;
 double ksi,eta;
 FloatArray   nx,ny;
 FloatMatrix  jacMtrx, inv ;
#ifdef  Quad1PlaneStrain_reducedShearIntegration
 GaussPoint  *helpGaussPoint;
 FloatArray  *coord;
#endif

 ksi = aGaussPoint -> giveCoordinate(1);
 eta = aGaussPoint -> giveCoordinate(2);
 
 this->giveDerivativeKsi (nx, eta);
 this->giveDerivativeEta (ny, ksi);
 
 this->computeJacobianMatrixAt(jacMtrx, aGaussPoint);
 inv.beInverseOf (jacMtrx);
 
 answer.resize(4,8);
 answer.zero();
 
 for (i=1;i<=4;i++){
  answer.at(1,2*i-1) = nx.at(i)*inv.at(1,1) + ny.at(i)*inv.at(1,2);
  answer.at(2,2*i-0) = nx.at(i)*inv.at(2,1) + ny.at(i)*inv.at(2,2);
 }

#ifdef  Quad1PlaneStrain_reducedShearIntegration

 coord = new FloatArray(2);
 coord -> at(1)    = 0.0;
 coord -> at(2)    = 0.0;
 helpGaussPoint     = new GaussPoint(this,1,coord,4.0,_PlaneStress);

 this->computeJacobianMatrixAt(jacMtrx, helpGaussPoint)
 inv.beInverseOf (jacMtrx);
 
 eta = helpGaussPoint -> giveCoordinate(1);
 ksi = helpGaussPoint -> giveCoordinate(2);
 
 this->giveDerivativeKsi(nx, eta);
 this->giveDerivativeEta(ny, ksi);

#endif
 
 for (i=1;i<=4;i++){
  answer.at(4,2*i-1) = nx.at(i)*inv.at(2,1) + ny.at(i)*inv.at(2,2);
  answer.at(4,2*i-0) = nx.at(i)*inv.at(1,1) + ny.at(i)*inv.at(1,2);
 }

#ifdef Quad1PlaneStrain_reducedShearIntegration
 delete helpGaussPoint;
#endif
 return ;
}



void
Quad1PlaneStrain :: computeGaussPoints ()
// Sets up the array containing the four Gauss points of the receiver.
{
   numberOfIntegrationRules = 1 ;
  integrationRulesArray = new IntegrationRule*[1];
  integrationRulesArray[0] = new GaussIntegrationRule (1,this, 1, 3);
  integrationRulesArray[0]->setUpIntegrationPoints (_Square, numberOfGaussPoints, _PlaneStrain);

}

void
Quad1PlaneStrain :: computeNmatrixAt (GaussPoint* aGaussPoint, FloatMatrix& answer) 
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at aGaussPoint.
{
 int          i;
   double       ksi,eta;
 // FloatArray*  n;
 // FloatMatrix* answer;

   ksi = aGaussPoint -> giveCoordinate(1) ;
   eta = aGaussPoint -> giveCoordinate(2) ;

 FloatArray n (4);

   n.at(1) = (1. + ksi) * (1. + eta) * 0.25 ;
   n.at(2) = (1. - ksi) * (1. + eta) * 0.25 ;
   n.at(3) = (1. - ksi) * (1. - eta) * 0.25 ;
   n.at(4) = (1. + ksi) * (1. - eta) * 0.25 ;

 answer.resize(2,8);
 answer.zero();

 for (i=1;i<=4;i++){
  answer.at(1,2*i-1) = n.at(i);
  answer.at(2,2*i-0) = n.at(i);
 }

 return ;
}



void
Quad1PlaneStrain :: computeEgdeNMatrixAt (FloatMatrix& answer, GaussPoint* aGaussPoint) 
{
  /*
  
  computes interpolation matrix for element edge.
  we assemble locally this matrix for only nonzero 
  shape functions. 
  (for example only two nonzero shape functions for 2 dofs are
  necessary for linear plane stress tringle edge).
  These nonzero shape functions are then mapped to 
  global element functions.
  
  Using mapping technique will allow to assemble shape functions 
  without regarding particular side
  */
 
  double ksi,n1,n2;
 answer.resize (2,4);
  answer.zero();

  ksi = aGaussPoint -> giveCoordinate(1) ;
  n1  = (1. - ksi) * 0.5 ;
  n2  = (1. + ksi) * 0.5 ;

  answer.at(1,1) = n1 ;
  answer.at(1,3) = n2 ;
  answer.at(2,2) = n1 ;
  answer.at(2,4) = n2 ;
  
  return  ;
}


void
Quad1PlaneStrain :: giveEdgeDofMapping (IntArray& answer, int iEdge) const
{
  /* 
   provides dof mapping of local edge dofs (only nonzero are taken into account)
   to global element dofs
  */

  answer.resize(4);
  if (iEdge == 1)  { // edge between nodes 1,2
  answer.at(1) = 1;
  answer.at(2) = 2;
  answer.at(3) = 3;
  answer.at(4) = 4;
  } else if (iEdge == 2) { // edge between nodes 2 3
  answer.at(1) = 3;
  answer.at(2) = 4;
  answer.at(3) = 5;
  answer.at(4) = 6;
  } else if (iEdge == 3) {// edge between nodes 3 4
  answer.at(1) = 5;
  answer.at(2) = 6;
  answer.at(3) = 7;
  answer.at(4) = 8;
  } else if (iEdge == 4) {// edge between nodes 4 1
  answer.at(1) = 7;
  answer.at(2) = 8;
  answer.at(3) = 1;
  answer.at(4) = 2;
  } else {
  _error ("giveEdgeDofMapping: wrong edge number");
  }

  return ;
}

double
Quad1PlaneStrain ::   computeEdgeVolumeAround (GaussPoint* aGaussPoint, int iEdge) 
{
  double dx,dy, length ;
  Node   *nodeA,*nodeB ;
  int aNode=0, bNode=0;

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
  return 0.5 * length * aGaussPoint -> giveWeight() ;

}

void      
Quad1PlaneStrain::computeEdgeIpGlobalCoords (FloatArray& answer, GaussPoint* gp, int iEdge) 
{

  double n1, n2, ksi = gp->giveCoordinate(1) ;
  Node   *nodeA,*nodeB ;
  int aNode=0, bNode=0;

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
  _error ("computeEdgeIpGlobalCoords: wrong egde number");
  }

 n1  = (1. - ksi) * 0.5 ;
 n2  = (1. + ksi) * 0.5 ;

  nodeA   = this->giveNode(aNode) ;
  nodeB   = this->giveNode(bNode) ;
  
 answer.resize (2);
 answer.at(1) = n1*nodeA->giveCoordinate(1)+n2*nodeB->giveCoordinate(1);
 answer.at(2) = n1*nodeA->giveCoordinate(2)+n2*nodeB->giveCoordinate(2);
}


int
Quad1PlaneStrain :: computeLoadLEToLRotationMatrix (FloatMatrix& answer, int iEdge, GaussPoint* gp) 
{
  // returns transformation matrix from
  // edge local coordinate system
  // to element local c.s 
  // (same as global c.s in this case) 
  //
  // i.e. f(element local) = T * f(edge local)
  //
  double dx,dy, length ;
  Node   *nodeA,*nodeB ;
  int aNode=0, bNode=0;

 answer.resize(2,2);
 answer.zero();

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

  answer.at(1,1) = dx/length;
  answer.at(1,2) = -dy/length;
  answer.at(2,1) = answer.at(1,2);
  answer.at(2,2) = dx/length;

  return 1;
}


double
Quad1PlaneStrain :: computeVolumeAround (GaussPoint* aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
   double       determinant,weight,thickness,volume;
  FloatMatrix jacMtrx;

  this->computeJacobianMatrixAt(jacMtrx, aGaussPoint);
   determinant = fabs (jacMtrx.giveDeterminant());
   weight      = aGaussPoint -> giveWeight();
   thickness   = this -> giveCrossSection() -> give('t');
   volume      = determinant * weight * thickness ;

   return volume;
}

void
Quad1PlaneStrain :: computeJacobianMatrixAt (FloatMatrix& answer, GaussPoint* aGaussPoint)
// Returns the jacobian matrix  J (x,y)/(ksi,eta)  of the receiver.
// Computes it if it does not exist yet.
{
 int        i;
 double     ksi,eta,x,y;
 FloatArray nx,ny;

 answer.resize(2,2); answer.zero();

   ksi = aGaussPoint -> giveCoordinate(1);
   eta = aGaussPoint -> giveCoordinate(2);

 this->giveDerivativeKsi(nx, eta);
 this->giveDerivativeEta(ny, ksi);

 for (i=1;i<=4;i++){
  x = this -> giveNode(i) -> giveCoordinate(1);
  y = this -> giveNode(i) -> giveCoordinate(2);

  answer.at(1,1) += nx.at(i)*x;
  answer.at(1,2) += nx.at(i)*y;
  answer.at(2,1) += ny.at(i)*x;
  answer.at(2,2) += ny.at(i)*y;
 }
}

IRResultType
Quad1PlaneStrain :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 this->StructuralElement :: initializeFrom (ir);
 numberOfGaussPoints = 4;
 IR_GIVE_OPTIONAL_FIELD (ir, numberOfGaussPoints, IFT_Quad1PlaneStrain_nip, "nip"); // Macro
 
 if (!((numberOfGaussPoints == 4) ||
   (numberOfGaussPoints == 9) ||
   (numberOfGaussPoints == 16)))
   numberOfGaussPoints = 4;

 this -> computeGaussPoints();
 return IRRT_OK;
}


double
Quad1PlaneStrain :: giveCharacteristicLenght (GaussPoint* gp, const FloatArray& normalToCrackPlane)
//
// returns receivers characteristic length in gp (for some material models)
// for crack formed in plane with normal normalToCrackPlane.
//
{
 // return this -> giveLenghtInDir(normalToCrackPlane) / sqrt (this->numberOfGaussPoints);
 return this -> giveLenghtInDir(normalToCrackPlane);
}

void
Quad1PlaneStrain ::   giveDofManDofIDMask  (int inode, EquationID, IntArray& answer) const {
// returns DofId mask array for inode element node.
// DofId mask array determines the dof ordering requsted from node.
// DofId mask array contains the DofID constants (defined in cltypes.h)
// describing physical meaning of particular DOFs.
 //IntArray* answer = new IntArray (2);
 answer.resize (2);

 answer.at(1) = D_u;
 answer.at(2) = D_v;

 return ;
}


/*
int
Quad1PlaneStrain :: computeGtoNRotationMatrix (FloatMatrix &answer)
   // returns transformation matrix from global coordinate set to
   // nodal coordinate set
{
 FloatMatrix *triplet;
 int i,flag=0,ii;
 
 for (i=1; i<= numberOfNodes; i++)
  flag += this->giveNode(i)->hasLocalCS ();
 if (flag == 0) {answer.beEmptyMtrx(); return 0 ;}
  
 answer.resize (this->computeNumberOfDofs(),this->computeNumberOfDofs());
 answer.zero();

 // loop over nodes
 for (i=1; i<= numberOfNodes; i++) {
  ii = (i-1)*2+1 ;
  if (this->giveNode(i)->hasLocalCS ()) {
   triplet = this->giveNode(i)->giveLocalCoordinateTriplet();
   answer.at(ii,ii)     = triplet->at(1,1);
   answer.at(ii,ii+1)   = triplet->at(1,2);
   answer.at(ii+1,ii)   = triplet->at(2,1);
   answer.at(ii+1,ii+1) = triplet->at(2,2);
  } else {
   // no transformation - unit matrix as 
   // transformation submatrix for node i
   answer.at(ii,ii)     = 1.0;
   answer.at(ii+1,ii+1) = 1.0;
  }
 }

 return 1;
}
*/

int
Quad1PlaneStrain :: computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) 
{
 double ksi, eta, n1, n2, n3, n4;
 
 ksi = lcoords.at(1) ;
 eta = lcoords.at(2) ;
 
 n1 = (1. + ksi) * (1. + eta) * 0.25 ;
 n2 = (1. - ksi) * (1. + eta) * 0.25 ;
 n3 = (1. - ksi) * (1. - eta) * 0.25 ;
 n4 = (1. + ksi) * (1. - eta) * 0.25 ;
 
 answer.resize(2);
 answer.at(1) = n1*this->giveNode(1)->giveCoordinate(1) + n2*this->giveNode(2)->giveCoordinate(1)+
  n3 * this->giveNode(3)->giveCoordinate(1) + n4 * this->giveNode(4)->giveCoordinate(1);
 answer.at(2) = n1*this->giveNode(1)->giveCoordinate(2) + n2*this->giveNode(2)->giveCoordinate(2)+
  n3 * this->giveNode(3)->giveCoordinate(2) + n4 * this->giveNode(4)->giveCoordinate(2);

 return 1;
}


Interface* 
Quad1PlaneStrain :: giveInterface (InterfaceType interface)
{
 if (interface == ZZNodalRecoveryModelInterfaceType) return (ZZNodalRecoveryModelInterface*) this;
 else if (interface == SPRNodalRecoveryModelInterfaceType) return (SPRNodalRecoveryModelInterface*) this;
 else if (interface == SpatialLocalizerInterfaceType) return (SpatialLocalizerInterface*) this;
 else if (interface == DirectErrorIndicatorRCInterfaceType) return (DirectErrorIndicatorRCInterface*) this;
 else if (interface == EIPrimaryUnknownMapperInterfaceType) return (EIPrimaryUnknownMapperInterface*) this;
 else if (interface == HuertaErrorEstimatorInterfaceType) return (HuertaErrorEstimatorInterface*)this;
 else if (interface == HuertaRemeshingCriteriaInterfaceType) return (HuertaRemeshingCriteriaInterface*) this;
 return NULL;
}


int 
Quad1PlaneStrain :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
 if ((type == IST_StressTensor)||(type ==IST_StrainTensor)) return 4;
 
 GaussPoint *gp = integrationRulesArray[0]-> getIntegrationPoint(0) ;
 return this->giveIPValueSize (type, gp);
}

void
Quad1PlaneStrain :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx  
(FloatMatrix& answer, GaussPoint* aGaussPoint, InternalStateType type)
{
  // evaluates N matrix (interpolation estimated stress matrix)
  // according to Zienkiewicz & Zhu paper
  // N(nsigma, nsigma*nnodes)
  // Definition : sigmaVector = N * nodalSigmaVector
   double ksi, eta, n1, n2, n3, n4;

   ksi = aGaussPoint -> giveCoordinate(1) ;
   eta = aGaussPoint -> giveCoordinate(2) ;

   n1 = (1. + ksi) * (1. + eta) * 0.25 ;
   n2 = (1. - ksi) * (1. + eta) * 0.25 ;
   n3 = (1. - ksi) * (1. - eta) * 0.25 ;
   n4 = (1. + ksi) * (1. - eta) * 0.25 ;

  if (this->giveIPValueSize(type, aGaussPoint)) answer.resize(1,4) ;
  else return;

  answer.zero();

   answer.at(1,1) = n1 ;
   answer.at(1,2) = n2 ;
   answer.at(1,3) = n3 ;
   answer.at(1,4) = n4 ;

   return ;
}

void
Quad1PlaneStrain::HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId, 
                                  IntArray &localNodeIdArray, IntArray &globalNodeIdArray, 
                                  HuertaErrorEstimatorInterface::SetupMode sMode, TimeStep* tStep,
                                  int &localNodeId, int &localElemId, int &localBcId,
                                  IntArray &controlNode, IntArray &controlDof,
                                  HuertaErrorEstimator::AnalysisMode aMode)
{
 Element *element = this->HuertaErrorEstimatorI_giveElement();
 int inode, nodes = 4, iside, sides = 4, nd1, nd2;
 FloatArray *corner[4], midSide[4], midNode, cor[4];
 double x = 0.0, y = 0.0;

 static int sideNode[4][2] = {{1, 2}, {2, 3}, {3, 4}, {4, 1}};

 if(sMode == HuertaErrorEstimatorInterface::NodeMode || 
   (sMode == HuertaErrorEstimatorInterface::BCMode && aMode == HuertaErrorEstimator::HEE_linear)){
  for(inode = 0; inode < nodes; inode++){
   corner[inode] = element -> giveNode(inode + 1) -> giveCoordinates();
   if(corner[inode] -> giveSize() != 3){
    cor[inode].resize(3);
    cor[inode].at(1) = corner[inode] -> at(1);
    cor[inode].at(2) = corner[inode] -> at(2);
    cor[inode].at(3) = 0.0;

    corner[inode] = &(cor[inode]);
   }

   x += corner[inode] -> at(1);
   y += corner[inode] -> at(2);
  }
  
  for(iside = 0; iside < sides; iside++){
   midSide[iside].resize(3);

   nd1 = sideNode[iside][0] - 1;
   nd2 = sideNode[iside][1] - 1;

   midSide[iside].at(1) = (corner[nd1] -> at(1) + corner[nd2] -> at(1)) / 2.0;
   midSide[iside].at(2) = (corner[nd1] -> at(2) + corner[nd2] -> at(2)) / 2.0;
   midSide[iside].at(3) = 0.0;
  }
  
  midNode.resize(3);

  midNode.at(1) = x / nodes;
  midNode.at(2) = y / nodes;
  midNode.at(3) = 0.0;
 }

 this -> setupRefinedElementProblem2D (element, refinedElement, level, nodeId, localNodeIdArray, globalNodeIdArray, 
                    sMode, tStep, nodes, corner, midSide, midNode, 
                    localNodeId, localElemId, localBcId,
                    controlNode, controlDof, aMode, "Quad1PlaneStrain");
}



#ifdef __OOFEG
#include "rcm2.h"
#define TR_LENGHT_REDUCT 0.3333

void Quad1PlaneStrain :: drawRawGeometry (oofegGraphicContext &gc)
{
  WCRec p[4];
  GraphicObj *go;

 if (!gc.testElementGraphicActivity(this)) return; 

  EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
  EASValsSetColor(gc.getElementColor());
  EASValsSetEdgeColor(gc.getElementEdgeColor());
  EASValsSetEdgeFlag(TRUE);
  EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
 EASValsSetFillStyle (FILL_HOLLOW);
  p[0].x = (FPNum) this->giveNode(1)->giveCoordinate(1);
  p[0].y = (FPNum) this->giveNode(1)->giveCoordinate(2);
  p[0].z = 0.;
  p[1].x = (FPNum) this->giveNode(2)->giveCoordinate(1);
  p[1].y = (FPNum) this->giveNode(2)->giveCoordinate(2);
  p[1].z = 0.;
  p[2].x = (FPNum) this->giveNode(3)->giveCoordinate(1);
  p[2].y = (FPNum) this->giveNode(3)->giveCoordinate(2);
  p[2].z = 0.;
  p[3].x = (FPNum) this->giveNode(4)->giveCoordinate(1);
  p[3].y = (FPNum) this->giveNode(4)->giveCoordinate(2);
  p[3].z = 0.;
   
  go =  CreateQuad3D(p);
  EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
  EGAttachObject(go, (EObjectP) this);
  EMAddGraphicsToModel(ESIModel(), go);
}


void Quad1PlaneStrain :: drawDeformedGeometry (oofegGraphicContext &gc, UnknownType type)
{
  WCRec p[4];
  GraphicObj *go;
  TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
 double defScale = gc.getDefScale();

 if (!gc.testElementGraphicActivity(this)) return; 

  EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
  EASValsSetColor(gc.getDeformedElementColor());
  EASValsSetEdgeColor(gc.getElementEdgeColor());
  EASValsSetEdgeFlag(TRUE);
  EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
 EASValsSetFillStyle (FILL_HOLLOW);
  p[0].x = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
  p[0].y = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
  p[0].z = 0.;
  p[1].x = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
  p[1].y = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
  p[1].z = 0.;
  p[2].x = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
  p[2].y = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
  p[2].z = 0.;
  p[3].x = (FPNum) this->giveNode(4)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
  p[3].y = (FPNum) this->giveNode(4)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
  p[3].z = 0.;
   
  go =  CreateQuad3D(p);
  EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
  EMAddGraphicsToModel(ESIModel(), go);
}



void Quad1PlaneStrain :: drawScalar   (oofegGraphicContext& context)
{
 int i, indx, result = 0;
  WCRec p[4];
  GraphicObj *tr;
 TimeStep* tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
 FloatArray v[4];
 double s[4], defScale;
 IntArray map;

 if (!context.testElementGraphicActivity(this)) return; 

 EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
 if (context.giveIntVarMode() == ISM_recovered) {
  for (i=1; i<=4; i++) 
   result+= this->giveInternalStateAtNode (v[i-1], context.giveIntVarType(), context.giveIntVarMode(), i, tStep);
  if (result !=4) return;

  this->giveIntVarCompFullIndx (map, context.giveIntVarType());
  if ((indx = map.at(context.giveIntVarIndx())) == 0) return;
  for (i=1; i<=4; i++) s[i-1]=v[i-1].at(indx);
  
  if (context.getScalarAlgo() == SA_ISO_SURF) {
   
   for (i=0; i< 4; i++) {
    if (context.getInternalVarsDefGeoFlag()) {
     // use deformed geometry
     defScale = context.getDefScale();
     p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
     p[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
     p[i].z = 0.;
     
    } else {
     p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
     p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
     p[i].z = 0.;
    }
   }
   
   //EASValsSetColor(gc.getYieldPlotColor(ratio));
   context.updateFringeTableMinMax (s, 4);
   tr =  CreateQuadWD3D(p, s[0], s[1], s[2], s[3]);
   EGWithMaskChangeAttributes(LAYER_MASK, tr);
   EMAddGraphicsToModel(ESIModel(), tr);
   
  } else if ((context.getScalarAlgo() == SA_ZPROFILE)||(context.getScalarAlgo() == SA_COLORZPROFILE)) {
   double landScale= context.getLandScale();
   
   for (i=0; i< 4; i++) {
    if (context.getInternalVarsDefGeoFlag()) {
     // use deformed geometry
     defScale = context.getDefScale();
     p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
     p[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
     p[i].z = s[i]*landScale;
     
    } else {
     p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
     p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
     p[i].z = s[i]*landScale;
    }
// this fixes a bug in ELIXIR
    if(fabs(s[i]) < 1.0e-6)s[i] = 1.0e-6;
   }

   if (context.getScalarAlgo() == SA_ZPROFILE) {
    EASValsSetColor(context.getDeformedElementColor());
    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    tr =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
   } else {
    context.updateFringeTableMinMax (s, 4);
    tr =  CreateQuadWD3D(p, s[0], s[1], s[2], s[3]);
    EGWithMaskChangeAttributes(LAYER_MASK, tr);
   }
   EMAddGraphicsToModel(ESIModel(), tr);
  } 
 } else if (context.giveIntVarMode() == ISM_local) {
  if (numberOfGaussPoints != 4) return;
  int ip;
  GaussPoint* gp;
  IntArray ind(4);
  FloatArray* gpCoords;
   WCRec pp[9];
  
   for (i=0; i< 4; i++) {
    if (context.getInternalVarsDefGeoFlag()) {
     // use deformed geometry
     defScale = context.getDefScale();
     pp[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
     pp[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
     pp[i].z = 0.;
     
    } else {
     pp[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
     pp[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
     pp[i].z = 0.;
    }
   }

  for (i=0; i<3; i++) {
   pp[i+4].x = 0.5*(pp[i].x+pp[i+1].x); pp[i+4].y = 0.5*(pp[i].y+pp[i+1].y); pp[i+4].z = 0.5*(pp[i].z+pp[i+1].z);
  }
  pp[7].x = 0.5*(pp[3].x+pp[0].x); pp[7].y = 0.5*(pp[3].y+pp[0].y); pp[7].z = 0.5*(pp[3].z+pp[0].z);

  pp[8].x = 0.25*(pp[0].x+pp[1].x+pp[2].x+pp[3].x); 
  pp[8].y = 0.25*(pp[0].y+pp[1].y+pp[2].y+pp[3].y);
  pp[8].z = 0.25*(pp[0].z+pp[1].z+pp[2].z+pp[3].z);
  
   for (ip=1 ; ip<= numberOfGaussPoints ; ip++) {
    gp = integrationRulesArray[0]-> getIntegrationPoint(ip-1) ;
   gpCoords = gp->giveCoordinates();
   if ((gpCoords->at(1) >0.)&&(gpCoords->at(2) >0.)) {
    ind.at(1) = 0; ind.at(2) = 4; ind.at(3) = 8; ind.at(4) = 7;
   } else if ((gpCoords->at(1) < 0.)&&(gpCoords->at(2) >0.)) {
    ind.at(1) = 4; ind.at(2) = 1; ind.at(3) = 5; ind.at(4) = 8;
   } else if ((gpCoords->at(1) < 0.)&&(gpCoords->at(2) <0.)) {
    ind.at(1) = 5; ind.at(2) = 2; ind.at(3) = 6; ind.at(4) = 8;
   } else {
    ind.at(1) = 6; ind.at(2) = 3; ind.at(3) = 7; ind.at(4) = 8;
   }
   
   if (giveIPValue (v[0], gp, context.giveIntVarType(), tStep) == 0 ) return;
   this->giveIntVarCompFullIndx (map, context.giveIntVarType());
   if ((indx = map.at(context.giveIntVarIndx())) == 0) return;
   for (i=1; i<=4; i++) s[i-1]=v[0].at(indx);
   
   for (i=0; i<4; i++) {
    p[i].x = pp[ind.at(i+1)].x;
    p[i].y = pp[ind.at(i+1)].y;
    p[i].z = pp[ind.at(i+1)].z;
   }
   
   context.updateFringeTableMinMax (s, 4);
   tr =  CreateQuadWD3D(p, s[0], s[1], s[2], s[3]);
   EGWithMaskChangeAttributes(LAYER_MASK, tr);
   EMAddGraphicsToModel(ESIModel(), tr);
  }
 }
 return;
}




void
Quad1PlaneStrain::drawSpecial(oofegGraphicContext& gc)
{
 int i;
  WCRec l[2];
  GraphicObj *tr;
  StructuralMaterial *mat = (StructuralMaterial*) this->giveMaterial();
  GaussPoint* gp;
 TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
 double defScale = gc.getDefScale();
 FloatArray crackStatuses, cf;

 if (!gc.testElementGraphicActivity(this)) return; 
 if (gc.giveIntVarType() == IST_CrackState) {
  // ask if any active crack exist
   int igp, crackStatus;
   double ax,ay,bx,by,norm,xc,yc,length;
   FloatArray crackDir;
  FloatArray gpglobalcoords;
  
   for (igp=1 ; igp<= numberOfGaussPoints ; igp++) {
    gp = integrationRulesArray[0]-> getIntegrationPoint(igp-1);
    
    if (mat->giveIPValue(cf, gp, CrackedFlag, tStep) == 0) return;
   if ((int) cf.at(1) == 0) return;
    if (mat->giveIPValue(crackDir, gp, CrackDirs, tStep)) {
    mat->giveIPValue(crackStatuses, gp, CrackStatuses, tStep);
    for (i = 1; i <= 3; i++) {
     crackStatus = (int) crackStatuses.at(i);
     if ((crackStatus != pscm_NONE) && (crackStatus != pscm_CLOSED)) {
      // draw a crack
      // this element is 2d element in x-y plane
      // 
      // compute perpendicular line to normal in xy plane
      ax = crackDir.at(i);
      ay = crackDir.at(3+i);
      if (fabs(ax) > 1.e-6) {
       by = 1.;
       bx = -ay*by/ax;
       norm = sqrt(bx*bx+by*by);
       bx = bx/norm;   // normalize to obtain unit vector
       by = by/norm;
      } else { by = 0.0; bx = 1.0;}
      // obtain gp global coordinates 
      if (gc.getInternalVarsDefGeoFlag()) {
       double ksi, eta, n1, n2, n3, n4;
       ksi = gp->giveCoordinate(1);
       eta = gp->giveCoordinate(2);

       n1 = (1. + ksi) * (1. + eta) * 0.25 ;
       n2 = (1. - ksi) * (1. + eta) * 0.25 ;
       n3 = (1. - ksi) * (1. - eta) * 0.25 ;
       n4 = (1. + ksi) * (1. - eta) * 0.25 ;

       gpglobalcoords.resize (2);
       
       gpglobalcoords.at(1) = (n1 * this->giveNode(1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale) + 
                   n2 * this->giveNode(2)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale) +
                   n3 * this->giveNode(3)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale) + 
                   n4 * this->giveNode(4)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale));
       gpglobalcoords.at(2) = (n1*this->giveNode(1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale) + 
                   n2*this->giveNode(2)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale) +
                   n3 * this->giveNode(3)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale)+ 
                   n4 * this->giveNode(4)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale));
      } else {
       computeGlobalCoordinates (gpglobalcoords, *(gp->giveCoordinates()));
      }
      xc = gpglobalcoords.at(1);
      yc = gpglobalcoords.at(2);
      length = sqrt(computeVolumeAround(gp)/this->giveCrossSection()->give('t'))/2.;
      l[0].x = (FPNum) xc+bx*length;
      l[0].y = (FPNum) yc+by*length;
      l[0].z = 0.;
      l[1].x = (FPNum) xc-bx*length;
      l[1].y = (FPNum) yc-by*length;
      l[1].z = 0.;
      EASValsSetLayer(OOFEG_CRACK_PATTERN_LAYER);
      EASValsSetLineWidth(OOFEG_CRACK_PATTERN_WIDTH);
      if ((crackStatus == pscm_SOFTENING)||(crackStatus == pscm_OPEN))
       EASValsSetColor(gc.getActiveCrackColor());
      else 
       EASValsSetColor(gc.getCrackPatternColor());
      tr = CreateLine3D (l);
      EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
      EMAddGraphicsToModel(ESIModel(), tr);     
     }
    }
   }
  }
 }
}

#endif
 
void 
Quad1PlaneStrain::SPRNodalRecoveryMI_giveSPRAssemblyPoints (IntArray& pap)
{
 pap.resize(4);
 for (int i=1; i <5; i++) pap.at(i) = this->giveNode(i)->giveNumber();
}

void
Quad1PlaneStrain::SPRNodalRecoveryMI_giveDofMansDeterminedByPatch (IntArray& answer, int pap)
{
 int found = 0;
 answer.resize(1);

 for (int i=1; i<5; i++) if (pap == this->giveNode(i)->giveNumber()) found = 1;
 if (found) answer.at(1) = pap;
 else _error ("SPRNodalRecoveryMI_giveDofMansDeterminedByPatch: node unknown");
}

int
Quad1PlaneStrain::SPRNodalRecoveryMI_giveNumberOfIP ()
{ return this->giveDefaultIntegrationRulePtr()->getNumberOfIntegrationPoints(); }


void
Quad1PlaneStrain::SPRNodalRecoveryMI_computeIPGlobalCoordinates (FloatArray& coords, GaussPoint* gp)
{
 this->computeGlobalCoordinates(coords, *gp->giveCoordinates());
}

SPRPatchType
Quad1PlaneStrain::SPRNodalRecoveryMI_givePatchType()
{
 return SPRPatchType_2dxy;
}


#define POINT_TOL 1.e-6

int
Quad1PlaneStrain::computeLocalCoordinates (FloatArray& answer, const FloatArray& coords)
{
   Node    *node1,*node2,*node3,*node4;
   double  x1,x2,x3,x4,y1,y2,y3,y4,a1,a2,a3,a4,b1,b2,b3,b4;
  double a,b,c, ksi1, ksi2, ksi3, eta1=0.0, eta2=0.0, denom;
  int nroot, xind = 1, yind = 2;

  answer.resize(2);

   node1 = this -> giveNode(1) ;
   node2 = this -> giveNode(2) ;
   node3 = this -> giveNode(3) ;
   node4 = this -> giveNode(4) ;

   x1 = node1 -> giveCoordinate(1) ;
   x2 = node2 -> giveCoordinate(1) ;
   x3 = node3 -> giveCoordinate(1) ;
   x4 = node4 -> giveCoordinate(1) ;

   y1 = node1 -> giveCoordinate(2) ;
   y2 = node2 -> giveCoordinate(2) ;
   y3 = node3 -> giveCoordinate(2) ;
   y4 = node4 -> giveCoordinate(2) ;

  a1 = x1+x2+x3+x4;
  a2 = x1-x2-x3+x4;
  a3 = x1+x2-x3-x4;
  a4 = x1-x2+x3-x4;

  b1 = y1+y2+y3+y4;
  b2 = y1-y2-y3+y4;
  b3 = y1+y2-y3-y4;
  b4 = y1-y2+y3-y4;

  a = a2*b4-b2*a4;
  b = a1*b4+a2*b3-a3*b2-b1*a4-b4*4.0*coords.at(1)+a4*4.0*coords.at(2);
  c = a1*b3-a3*b1-4.0*coords.at(1)*b3+4.0*coords.at(2)*a3;

  // solve quadratic equation for ksi
  cubic (0.0, a, b, c, &ksi1, &ksi2, &ksi3, &nroot);

  if (nroot == 0) return 0;
  if(nroot){
   denom = (b3+ksi1*b4);
   if (fabs(denom) <= 1.0e-10) 
    eta1 = (4.0*coords.at(xind)-a1-ksi1*a2) / (a3+ksi1*a4);
   else eta1 = (4.0*coords.at(yind)-b1-ksi1*b2) / denom;
  }
  if(nroot > 1){
   double diff_ksi1, diff_eta1, diff_ksi2, diff_eta2, diff1, diff2;

   denom = b3+ksi2*b4;
   if (fabs(denom) <= 1.0e-10) 
    eta2 = (4.0*coords.at(xind)-a1-ksi2*a2) / (a3+ksi2*a4);
   else eta2 = (4.0*coords.at(yind)-b1-ksi2*b2) / denom;
   
   // choose the one which seems to be closer to the parametric space (square <-1;1>x<-1;1>)
   diff_ksi1 = 0.0;
   if(ksi1 > 1.0)diff_ksi1 = ksi1 - 1.0;
   if(ksi1 < -1.0)diff_ksi1 = ksi1 + 1.0;

   diff_eta1 = 0.0;
   if(eta1 > 1.0)diff_eta1 = eta1 - 1.0;
   if(eta1 < -1.0)diff_eta1 = eta1 + 1.0;
   
   diff_ksi2 = 0.0;
   if(ksi2 > 1.0)diff_ksi2 = ksi2 - 1.0;
   if(ksi2 < -1.0)diff_ksi2 = ksi2 + 1.0;

   diff_eta2 = 0.0;
   if(eta2 > 1.0)diff_eta2 = eta2 - 1.0;
   if(eta2 < -1.0)diff_eta2 = eta2 + 1.0;
   
   diff1 = diff_ksi1 * diff_ksi1 + diff_eta1 * diff_eta1;
   diff2 = diff_ksi2 * diff_ksi2 + diff_eta2 * diff_eta2;

   // ksi2, eta2 seems to be closer
   if(diff1 > diff2){
    ksi1 = ksi2;
    eta1 = eta2;
   }
  }

  answer.at(1) = ksi1;
  answer.at(2) = eta1;

  // test if inside
  for (int i=1; i<=2; i++) {
   if (answer.at(i)<(-1.-POINT_TOL)) return 0;
   if (answer.at(i)>(1.+POINT_TOL)) return 0;
  }

  return 1;
}

int 
Quad1PlaneStrain::SpatialLocalizerI_containsPoint (const FloatArray& coords) {

 int result;
 FloatArray lcoords;
 result = this->computeLocalCoordinates (lcoords, coords);

 return result;
}

double 
Quad1PlaneStrain::SpatialLocalizerI_giveDistanceFromParametricCenter (const FloatArray& coords)
{
 FloatArray lcoords(2), gcoords;
 double dist;
 int size, gsize;

 lcoords.at(1) = lcoords.at(2) = 0.0;
 this -> computeGlobalCoordinates (gcoords, lcoords);

 if((size = coords.giveSize()) < (gsize = gcoords.giveSize()))
  _error("SpatialLocalizerI_giveDistanceFromParametricCenter: coordinates size mismatch");

 if(size == gsize){
  dist = coords.distance(gcoords);
 }
 else{
  FloatArray helpCoords = coords;

  helpCoords.resize(gsize);
  dist = helpCoords.distance(gcoords);
 }

 return dist;
}


double 
Quad1PlaneStrain::DirectErrorIndicatorRCI_giveCharacteristicSize () {

 int i;
 IntegrationRule* iRule;
 GaussPoint* gp ;
 double volume = 0.0;

 iRule = integrationRulesArray[giveDefaultIntegrationRule()];
 for (i=0 ; i < iRule->getNumberOfIntegrationPoints() ; i++) {
  gp  = iRule->getIntegrationPoint(i) ;
  volume += this -> computeVolumeAround(gp) ;
 } 
 volume /= this -> giveCrossSection() -> give('t');

 return sqrt(volume);
}


int
Quad1PlaneStrain::EIPrimaryUnknownMI_computePrimaryUnknownVectorAt (ValueModeType mode,
                                 TimeStep* stepN, const FloatArray& coords,
                                 FloatArray& answer)
{
 FloatArray lcoords, u;
 FloatMatrix n;
 double ksi, eta;
 int result;

 result = this->computeLocalCoordinates (lcoords, coords);

 ksi = lcoords.at(1);
 eta = lcoords.at(2);

 n.resize (2,8);
 n.zero();
  
 n.at(1,1) = n.at(2,2) = 0.25*(1.0+ksi)*(1.0+eta);
 n.at(1,3) = n.at(2,4) = 0.25*(1.0-ksi)*(1.0+eta);
 n.at(1,5) = n.at(2,6) = 0.25*(1.0-ksi)*(1.0-eta);
 n.at(1,7) = n.at(2,8) = 0.25*(1.0+ksi)*(1.0-eta);
 
 this->computeVectorOf (EID_MomentumBalance, mode, stepN, u);
 answer.beProductOf (n,u);

 return result;
}

void 
Quad1PlaneStrain::EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID (IntArray& answer)
{
 giveDofManDofIDMask  (1, EID_MomentumBalance, answer);
}
