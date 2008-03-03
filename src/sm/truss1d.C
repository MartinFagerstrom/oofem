/* $Header: /home/cvs/bp/oofem/sm/src/truss1d.C,v 1.6 2003/04/06 14:08:32 bp Exp $ */
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

//   file Truss1d.C

#include "truss1d.h"
#include "domain.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "cltypes.h"
#include "engngm.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <math.h>
#endif

#ifdef __OOFEG 
#include "oofeggraphiccontext.h"
#endif

Truss1d :: Truss1d (int n, Domain* aDomain) 
: StructuralElement (n,aDomain), 
 ZZNodalRecoveryModelInterface(), NodalAveragingRecoveryModelInterface(),
 SpatialLocalizerInterface()
, DirectErrorIndicatorRCInterface(),
 EIPrimaryUnknownMapperInterface(), ZZErrorEstimatorInterface(), ZZRemeshingCriteriaInterface(),
 MMAShapeFunctProjectionInterface(), HuertaErrorEstimatorInterface(), HuertaRemeshingCriteriaInterface()
// Constructor.
{
   numberOfDofMans     = 2 ;
   length              = 0. ;
}


void
Truss1d :: computeBmatrixAt (GaussPoint* aGaussPoint, FloatMatrix& answer, int li, int ui)
   // 
   // Returns linear part of geometrical equations of the receiver at gp.
   // Returns the linear part of the B matrix
   //
{
 double l, x1, x2;
 // FloatMatrix* answer;
 l = this->giveLength();

 x1 = this->giveNode(1)->giveCoordinate(1);
 x2 = this->giveNode(2)->giveCoordinate(1);

 answer.resize(1,2);
 
 answer.at(1,1) = (x1-x2)/l/l;
 answer.at(1,2) = (x2-x1)/l/l;

 return  ;
}

void  Truss1d :: computeGaussPoints ()
   // Sets up the array of Gauss Points of the receiver.
{

  numberOfIntegrationRules = 1 ;
  integrationRulesArray = new IntegrationRule*[1];
  integrationRulesArray[0] = new GaussIntegrationRule (1,this, 1, 2);
  integrationRulesArray[0]->setUpIntegrationPoints (_Line, 1, _1dMat);

}



void
Truss1d :: computeLumpedMassMatrix (FloatMatrix& answer, TimeStep* tStep)
   // Returns the lumped mass matrix of the receiver. This expression is
   // valid in both local and global axes.
{
   Material* mat ;
   double    halfMass ;

   mat        = this -> giveMaterial() ;
   halfMass   = mat->give('d') * this->giveCrossSection()->give('A') * this->giveLength() / 2.;
   answer.resize (2,2) ; answer.zero();
   answer . at(1,1) = halfMass ;
   answer . at(2,2) = halfMass ;

  //if (this->updateRotationMatrix()) answer.rotatedWith(*this->rotationMatrix) ;
   return  ;
}


void
Truss1d :: computeNmatrixAt (GaussPoint* aGaussPoint, FloatMatrix& answer) 
   // Returns the displacement interpolation matrix {N} of the receiver, eva-
   // luated at aGaussPoint.
{
   double       ksi,n1,n2 ;
   //FloatMatrix* answer ;

   ksi = aGaussPoint -> giveCoordinate(1) ;
   n1  = (1. - ksi) * 0.5 ;
   n2  = (1. + ksi) * 0.5 ;
   //answer = new FloatMatrix(2,4) ;
  answer.resize (1,2);
  answer.zero();

   answer.at(1,1) = n1 ;
   answer.at(1,2) = n2 ;

   return  ;
}

int
Truss1d :: computeGlobalCoordinates (FloatArray& answer, const FloatArray& lcoords) 
{
 double       ksi,n1,n2 ;

 ksi = lcoords.at(1) ;
 n1  = (1. - ksi) * 0.5 ;
 n2  = (1. + ksi) * 0.5 ;
 
 answer.resize (1);
 answer.at(1) = n1*this->giveNode(1)->giveCoordinate(1)+n2*this->giveNode(2)->giveCoordinate(1);
 
 return 1;
}


double  
Truss1d :: computeVolumeAround (GaussPoint* aGaussPoint)
   // Returns the length of the receiver. This method is valid only if 1
   // Gauss point is used.
{
 double weight  = aGaussPoint -> giveWeight() ;
   return 0.5 * this->giveLength() * weight * this->giveCrossSection()->give('A');
}


double
Truss1d :: giveLength ()
   // Returns the length of the receiver.
{
   double dx ;
   Node   *nodeA,*nodeB ;

   if (length == 0.) {
   nodeA   = this->giveNode(1) ;
   nodeB   = this->giveNode(2) ;
   dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1) ;
   length  = fabs(dx);
  }
   return length ;
 }


IRResultType
Truss1d :: initializeFrom (InputRecord* ir)
{
 this->StructuralElement :: initializeFrom (ir);
 this -> computeGaussPoints();
 return IRRT_OK;
}

 
void
Truss1d ::   giveDofManDofIDMask  (int inode, EquationID, IntArray& answer) const {
 // returns DofId mask array for inode element node.
 // DofId mask array determines the dof ordering requsted from node.
 // DofId mask array contains the DofID constants (defined in cltypes.h)
 // describing physical meaning of particular DOFs.
 //IntArray* answer = new IntArray (2);
 answer.resize (1);
 
 answer.at(1) = D_u;
 return ;
}




void
Truss1d::HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId, 
                             IntArray &localNodeIdArray, IntArray &globalNodeIdArray, 
                             HuertaErrorEstimatorInterface::SetupMode sMode, TimeStep* tStep,
                             int &localNodeId, int &localElemId, int &localBcId,
                             IntArray &controlNode, IntArray &controlDof,
                             HuertaErrorEstimator::AnalysisMode aMode)
{
 Element *element = this->HuertaErrorEstimatorI_giveElement();
 int inode, nodes = 2;
 FloatArray *corner[2], midNode, cor[2];
 double x = 0.0;

 if(sMode == HuertaErrorEstimatorInterface::NodeMode || 
   (sMode == HuertaErrorEstimatorInterface::BCMode && aMode == HuertaErrorEstimator::HEE_linear)){
  for(inode = 0; inode < nodes; inode++){
   corner[inode] = element -> giveNode(inode + 1) -> giveCoordinates();
   if(corner[inode] -> giveSize() != 3){
    cor[inode].resize(3);
    cor[inode].at(1) = corner[inode] -> at(1);
    cor[inode].at(2) = 0.0;
    cor[inode].at(3) = 0.0;

    corner[inode] = &(cor[inode]);
   }

   x += corner[inode] -> at(1);
  }
  
  midNode.resize(3);

  midNode.at(1) = x / nodes;
  midNode.at(2) = 0.0;
  midNode.at(3) = 0.0;
 }

 this -> setupRefinedElementProblem1D (element, refinedElement, level, nodeId, localNodeIdArray, globalNodeIdArray, 
                    sMode, tStep, nodes, corner, midNode, 
                    localNodeId, localElemId, localBcId,
                    controlNode, controlDof, aMode, "Truss1d");
}




#ifdef __OOFEG
void Truss1d :: drawRawGeometry (oofegGraphicContext& gc)
{
  GraphicObj *go;
//  if (!go) { // create new one
 WCRec p[2];   /* poin */
 if (!gc.testElementGraphicActivity(this)) return; 

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getElementColor());
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p[0].x = (FPNum) this->giveNode(1)->giveCoordinate(1);
    p[0].y = 0.;
    p[0].z = 0.0;
    p[1].x = (FPNum) this->giveNode(2)->giveCoordinate(1);
    p[1].y = 0.;
    p[1].z = 0.0;
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, (EObjectP) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

 
void Truss1d :: drawDeformedGeometry (oofegGraphicContext& gc, UnknownType type)
{
  GraphicObj *go;
  TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
 double defScale = gc.getDefScale();
  //  if (!go) { // create new one
  WCRec p[2];   /* poin */
 if (!gc.testElementGraphicActivity(this)) return; 

  EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
  EASValsSetColor(gc.getDeformedElementColor());
  EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
  p[0].x = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
  p[0].y = 0.;
  p[0].z = 0.;

  p[1].x = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
  p[1].y = 0.;
  p[1].z = 0.;
  go = CreateLine3D(p);
  EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
  EMAddGraphicsToModel(ESIModel(), go);
} 


void Truss1d :: drawScalar   (oofegGraphicContext& context)
{
 int i, indx, result = 0;
 WCRec p[2];
  GraphicObj *tr;
 TimeStep* tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
 FloatArray v1,v2;
 double s[2], defScale;
 IntArray map;

 if (!context.testElementGraphicActivity(this)) return; 
 if (context.giveIntVarMode() == ISM_recovered) {
  result+= this->giveInternalStateAtNode (v1, context.giveIntVarType(), context.giveIntVarMode(), 1, tStep);
  result+= this->giveInternalStateAtNode (v2, context.giveIntVarType(), context.giveIntVarMode(), 2, tStep);
 } else if (context.giveIntVarMode() == ISM_local) {
  GaussPoint* gp = integrationRulesArray[0]-> getIntegrationPoint(0);
  result+= giveIPValue (v1, gp, context.giveIntVarType(), tStep);
  v2 = v1; 
  result *= 2;
 }
 if (result != 2) return;

 this->giveIntVarCompFullIndx (map, context.giveIntVarType());

 if ((indx = map.at(context.giveIntVarIndx())) == 0) return;

 s[0] = v1.at(indx);
 s[1] = v2.at(indx);

 EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

 if ((context.getScalarAlgo() == SA_ISO_SURF)|| (context.getScalarAlgo() == SA_ISO_LINE)) {
  
  for (i=0; i< 2; i++) {
   if (context.getInternalVarsDefGeoFlag()) {
    // use deformed geometry
    defScale = context.getDefScale();
    p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
    p[i].y = 0.;
    p[i].z = 0.;
    
   } else {
    p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
    p[i].y = 0.;
    p[i].z = 0.;
   }
  }
  
  //EASValsSetColor(gc.getYieldPlotColor(ratio));
  tr =  CreateLine3D(p);
  EGWithMaskChangeAttributes(LAYER_MASK, tr);
  EMAddGraphicsToModel(ESIModel(), tr);

 } else if ((context.getScalarAlgo() == SA_ZPROFILE)||(context.getScalarAlgo() == SA_COLORZPROFILE)) {
  double landScale= context.getLandScale();
 
  for (i=0; i< 2; i++) {
   if (context.getInternalVarsDefGeoFlag()) {
    // use deformed geometry
    defScale = context.getDefScale();
    p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
    p[i].y = 0.0;
    p[i].z = s[i]*landScale;
     
   } else {
    p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
    p[i].y = 0.0;
    p[i].z = s[i]*landScale;
   }
  }
   
  if (context.getScalarAlgo() == SA_ZPROFILE) {
   /*
   EASValsSetColor(context.getDeformedElementColor());
   EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
   tr =  CreateLine3D(p);
   EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
   */
   WCRec pp[4];
   pp[0].x = p[0].x; pp[0].y = 0.0; pp[0].z = 0.0;
   pp[1].x = p[0].x; pp[1].y = 0.0; pp[1].z = p[0].z;
   pp[2].x = p[1].x; pp[2].y = 0.0; pp[2].z = p[1].z;
   pp[3].x = p[1].x; pp[3].y = 0.0; pp[3].z = 0.0;
   tr = CreateQuad3D(pp);
   EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
   EASValsSetColor(context.getDeformedElementColor());
   //EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
   EASValsSetFillStyle (FILL_HOLLOW);
   EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | LAYER_MASK, tr);
   EMAddGraphicsToModel(ESIModel(), tr);
  
  } else {
   //tr =  CreateTriangleWD3D(p, s[0], s[1], s[2]);
   EASValsSetColor(context.getDeformedElementColor());
   tr =  CreateLine3D(p);
   EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
   EMAddGraphicsToModel(ESIModel(), tr);
  }
 }
}



#endif


Interface* 
Truss1d ::giveInterface (InterfaceType interface)
{
 if (interface == ZZNodalRecoveryModelInterfaceType) return (ZZNodalRecoveryModelInterface*) this;
 else if (interface == NodalAveragingRecoveryModelInterfaceType) return (NodalAveragingRecoveryModelInterface*) this;
 else if (interface == SpatialLocalizerInterfaceType) return (SpatialLocalizerInterface*) this;
 else if (interface == DirectErrorIndicatorRCInterfaceType) return (DirectErrorIndicatorRCInterface*) this;
 else if (interface == EIPrimaryUnknownMapperInterfaceType) return (EIPrimaryUnknownMapperInterface*) this;
 else if (interface == ZZErrorEstimatorInterfaceType) return (ZZErrorEstimatorInterface*)this;
 else if (interface == ZZRemeshingCriteriaInterfaceType) return (ZZRemeshingCriteriaInterface*) this;
 else if (interface == MMAShapeFunctProjectionInterfaceType) return (MMAShapeFunctProjectionInterface*) this;
 else if (interface == HuertaErrorEstimatorInterfaceType) return (HuertaErrorEstimatorInterface*)this;
 else if (interface == HuertaRemeshingCriteriaInterfaceType) return (HuertaRemeshingCriteriaInterface*) this;
 return NULL;
}
 
int 
Truss1d :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
 if ((type == IST_StressTensor)||(type ==IST_StrainTensor)) return 1;

 GaussPoint *gp = integrationRulesArray[0]-> getIntegrationPoint(0) ;
 return this->giveIPValueSize (type, gp);
}


void
Truss1d :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx  
(FloatMatrix& answer, GaussPoint* aGaussPoint, InternalStateType type)
{
  // evaluates N matrix (interpolation estimated stress matrix)
  // according to Zienkiewicz & Zhu paper
  // N(nsigma, nsigma*nnodes)
  // Definition : sigmaVector = N * nodalSigmaVector
   double ksi, n1,n2;

   ksi = aGaussPoint -> giveCoordinate(1) ;
   n1  = (1. - ksi) * 0.5 ;
   n2  = (1. + ksi) * 0.5 ;

  if (this->giveIPValueSize(type, aGaussPoint)) answer.resize(1,2) ;
  else return;
  
  answer.at(1,1) = n1 ;
   answer.at(1,2) = n2 ;
  
   return ;
}

void 
Truss1d::NodalAveragingRecoveryMI_computeNodalValue (FloatArray& answer, int node,
                               InternalStateType type, TimeStep* tStep)
{
 GaussPoint* gp;
 gp = integrationRulesArray[0]-> getIntegrationPoint(0) ;
 this->giveIPValue (answer, gp, type, tStep);
/*
 if (type == IST_StressTensor) {
  gp = integrationRulesArray[0]-> getIntegrationPoint(0) ;
  answer = ((StructuralMaterialStatus*) this->giveMaterial()->giveStatus(gp)) -> giveStressVector();
 } else if (type == IST_StrainTensor) {
  gp = integrationRulesArray[0]-> getIntegrationPoint(0) ;
  answer = ((StructuralMaterialStatus*) this->giveMaterial()->giveStatus(gp)) -> giveStrainVector();
 }else answer.resize(0);
*/
}

void 
Truss1d::NodalAveragingRecoveryMI_computeSideValue (FloatArray& answer, int side,
                              InternalStateType type, TimeStep* tStep)
{
 answer.resize(0);
}


#define POINT_TOL 1.e-3

int
Truss1d::computeLocalCoordinates (FloatArray& answer, const FloatArray& gcoords) 
{
   Node    *node1,*node2;
   double  ksi, x1,x2;

  answer.resize(1);

   node1 = this -> giveNode(1) ;
   node2 = this -> giveNode(2) ;

   x1 = node1 -> giveCoordinate(1) ;
   x2 = node2 -> giveCoordinate(1) ;

  answer.at(1) = ksi = (2.0*gcoords.at(1)-(x1+x2))/(x2-x1);
  
  if (ksi<(-1.-POINT_TOL)) return 0;
  if (ksi>(1.+POINT_TOL)) return 0;
  return 1;
}

int 
Truss1d::SpatialLocalizerI_containsPoint (const FloatArray& coords) 
{
 FloatArray lc; 
 return this->computeLocalCoordinates (lc, coords);
}

double 
Truss1d::SpatialLocalizerI_giveDistanceFromParametricCenter (const FloatArray& coords)
{
 FloatArray lcoords(1), gcoords;
 double dist;
 int size, gsize;

 lcoords.at(1) = 0.0;
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
Truss1d::DirectErrorIndicatorRCI_giveCharacteristicSize () 
{
 return this->giveLength();
}

int
Truss1d::EIPrimaryUnknownMI_computePrimaryUnknownVectorAt (ValueModeType mode,
                              TimeStep* stepN, const FloatArray& coords,
                              FloatArray& answer)
{
 FloatArray u, ksi;
 FloatMatrix n;
 int result;

 result = this->computeLocalCoordinates (ksi, coords);

 n.resize (1,2);
 n.zero();
 
 n.at(1,1) = (1-ksi.at(1))*0.5;
 n.at(1,2) = (1+ksi.at(1))*0.5 ;

 this->computeVectorOf (EID_MomentumBalance, mode, stepN, u);
 answer.beProductOf (n,u);

 return result;
}

void 
Truss1d::EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID (IntArray& answer)
{
 giveDofManDofIDMask  (1, EID_MomentumBalance, answer);
}


void 
Truss1d::MMAShapeFunctProjectionInterface_interpolateIntVarAt (FloatArray& answer, FloatArray& coords,
                                coordType ct, nodalValContainerType& list,
                                InternalStateType type, TimeStep* tStep)
{
 int i, n;
 double n1, n2;
 FloatArray ksi(1);
 if (ct == MMAShapeFunctProjectionInterface::coordType_local)
  ksi.at(1) = coords.at(1);
 else 
  computeLocalCoordinates (ksi, coords);
 
 n1  = (1. - ksi.at(1)) * 0.5 ;
 n2  = (1. + ksi.at(1)) * 0.5 ;

 n = list.at(1)->giveSize();
 answer.resize(n);

 for (i=1; i<= n; i++) 
  answer.at(i) = n1*list.at(1)->at(i)+n2*list.at(2)->at(i);
 
 return ;
}
