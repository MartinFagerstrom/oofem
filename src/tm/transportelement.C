/* $Header: /home/cvs/bp/oofem/tm/src/transportelement.C,v 1.3.4.1 2004/04/05 15:19:53 bp Exp $ */
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

//   file transportelement.C

#include "transportelement.h"
#include "domain.h"
#include "timestep.h"
#include "node.h"
#include "dof.h"
#include "transportmaterial.h"
#include "load.h"
#include "boundaryload.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "debug.h"
#include "verbose.h"
#include "cltypes.h"
#include "elementside.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <stdio.h>
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "conTable.h"
#endif

TransportElement :: TransportElement (int n, Domain* aDomain, ElementMode em)
: Element (n, aDomain)
   // Constructor. Creates an element with number n, belonging to aDomain.
{
  emode = em;
}


TransportElement :: ~TransportElement ()
   // Destructor.
{}


void
TransportElement ::   giveElementDofIDMask  (EquationID, IntArray& answer) const {
// returns DofId mask array for inode element node.
// DofId mask array determines the dof ordering requsted from node.
// DofId mask array contains the DofID constants (defined in cltypes.h)
// describing physical meaning of particular DOFs.
 if (emode == HeatTransferEM) {
  answer.resize (1);
  answer.at(1) = T_f;
 } else if (emode == HeatMass1TransferEM) {
  answer.resize (2);
  answer.at(1) = T_f;
  answer.at(2) = C_1;
 } else {
  _error ("Unknown ElementMode");
 }
}


void
TransportElement ::  giveCharacteristicMatrix (FloatMatrix& answer, 
                        CharType mtrx, TimeStep *tStep) 
// 
// returns characteristics matrix of receiver accordind to mtrx
//
{
  if (mtrx == ConductivityMatrix) 
    this -> computeConductivityMatrix(answer, Conductivity, tStep); 
  else if (mtrx == CapacityMatrix) 
    this -> computeCapacityMatrix(answer, tStep);
  else if (mtrx == LHSBCMatrix)
    this -> computeBCMtrxAt (answer, tStep, VM_Total);
  else if (mtrx == IntSourceLHSMatrix)
    this -> computeIntSourceLHSMatrix (answer, tStep);
  else _error2("giveCharacteristicMatrix: Unknown Type of characteristic mtrx (%s)", __CharTypeToString(mtrx));
  
  return ;
}

void
TransportElement ::  giveCharacteristicVector (FloatArray& answer, CharType mtrx, ValueModeType mode,
                        TimeStep *tStep) 
// 
// returns characteristics vector of receiver according to requested type
//
{
  if (mtrx == ElementBCTransportVector) this->computeBCVectorAt (answer, tStep, mode); 
  else if (mtrx == ElementInternalSourceVector) this->computeInternalSourceRhsVectorAt (answer, tStep, mode);
  else _error2("giveCharacteristicVector: Unknown Type of characteristic mtrx (%s)",
	       __CharTypeToString(mtrx));
  
  return ;
}

int
TransportElement :: checkConsistency ()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
  int result =1;
  if (!this->giveMaterial()->testMaterialExtension(Material_TransportCapability)) {
    _warning("checkConsistency : material without support for transport problems");
    result =0;
  }
  /*
     if (!this->giveCrossSection()->testCrossSectionExtension(CS_TransportCapability)) {
     this->warning("checkConsistency : cross-section without support for transport problems", 1);
     result =0;
     }
     */
  return result;
}

void
TransportElement :: printOutputAt (FILE * file, TimeStep* stepN)
   // Performs end-of-step operations.
{
  Element::printOutputAt(file, stepN);
}

void
TransportElement :: computeCapacityMatrix (FloatMatrix& answer, TimeStep* tStep)
{
 answer.resize(computeNumberOfDofs(EID_ConservationEquation), computeNumberOfDofs(EID_ConservationEquation)); 
 answer.zero();

 if (emode == HeatTransferEM) {
   this->computeCapacitySubMatrix (answer, Capacity, 0, tStep);
 } else if (emode == HeatMass1TransferEM) {
   FloatMatrix subAnswer;
   int i;
   MatResponseMode rmode[2]={Capacity_hh, Capacity_ww};
   double coeff = 1.0; //this->giveMaterial()->give('d');
   
   for (i=1; i<=2; i++) {
     this->computeCapacitySubMatrix (subAnswer, rmode[i-1], 0, tStep);
     this->assembleLocalContribution (answer, subAnswer, 2, i, i, coeff);
   }
 } else {
   _error ("Unknown ElementMode");
 }
}

void
TransportElement :: computeConductivityMatrix (FloatMatrix& answer, MatResponseMode rMode, TimeStep* tStep)
{
  answer.resize(computeNumberOfDofs(EID_ConservationEquation), computeNumberOfDofs(EID_ConservationEquation)); 
  answer.zero();
  if (emode == HeatTransferEM) {
    this->computeConductivitySubMatrix (answer, 2, 0, Conductivity_hh, tStep);
  } else if (emode == HeatMass1TransferEM) {
    FloatMatrix subAnswer;
    MatResponseMode rmode[2][2]={{Conductivity_hh, Conductivity_hw},{Conductivity_wh, Conductivity_ww}};
    int i, j;
    
    for (i=1; i<=2; i++) 
      for (j=1; j<=2; j++) {
        this->computeConductivitySubMatrix (subAnswer, 2, 0, rmode[i-1][j-1], tStep);
        this->assembleLocalContribution (answer, subAnswer, 2, i, j, 1.0);
      }    
  } else {
    _error ("Unknown ElementMode");
  }
}



void
TransportElement::computeCapacitySubMatrix (FloatMatrix& answer, MatResponseMode rmode, int iri, TimeStep* tStep)
{
  int         i ;
  double      dV, c ;
  FloatMatrix n ;
  GaussPoint  *gp ;
  IntegrationRule* iRule = integrationRulesArray[iri];
  
  answer.resize (0,0);
  answer.zero();
  for (i=0 ; i<iRule->getNumberOfIntegrationPoints() ; i++) {
    gp      = iRule->getIntegrationPoint(i) ;
    this -> computeNSubMatrixAt(n, gp->giveCoordinates()) ;
    // ask for capacity coefficient
    c = ((TransportMaterial*)this->giveMaterial())->giveCharacteristicValue(rmode,gp,tStep);
    dV      = this -> computeVolumeAround(gp) ;
    answer.plusProductSymmUpper(n, n, dV * c) ;
  }
  
  answer.symmetrized() ;
}

void
TransportElement::computeConductivitySubMatrix (FloatMatrix& answer, int nsd, int iri, MatResponseMode rmode, TimeStep* tStep)
{
  int         i;
  double      dV ;
  FloatMatrix b, d, db;
  GaussPoint  *gp ;
  IntegrationRule* iRule = integrationRulesArray[iri];
  
  answer.resize (this->giveNumberOfDofManagers(),this->giveNumberOfDofManagers());
  answer.zero();
  for (i=0 ; i < iRule->getNumberOfIntegrationPoints() ; i++) {
    gp = iRule-> getIntegrationPoint(i) ;
    this->computeConstitutiveMatrixAt (d,rmode,gp,tStep);
    this -> computeGradientMatrixAt(b, gp) ;
    dV = this -> computeVolumeAround(gp) ;
    
    db.beProductOf(d,b);
    answer.plusProductSymmUpper(b, db, dV) ;
    //answer.plusProductUnsym(b,db,dV) ;
  }
  answer.symmetrized() ;
  return  ;
}

void
TransportElement::computeInternalSourceRhsSubVectorAt (FloatArray& answer, TimeStep* atTime, ValueModeType mode, int indx)
{
 // Computes numerically the generator Rhs vector of the receiver due to the generator
 //  at stepN.
 // // load is firrst transformed to local cs.
 // // load vector is then transformed to coordinate system in each node.
 // // (should be global coordinate system, but there may be defined 
 // //  different coordinate system in each node)
  int         i,j,igp,n,nLoads ;
  double dV;
  bcGeomType  ltype;
  Load*       load;
  IntegrationRule* iRule = integrationRulesArray[giveDefaultIntegrationRule()];
  TransportMaterial* mat = ((TransportMaterial*)this->giveMaterial());
  GaussPoint* gp;
  
  
  FloatArray  val, helpLoadVector, globalIPcoords ;
  FloatMatrix nm;
  answer.resize (0) ;
  
  nLoads    = this -> giveBodyLoadArray()->giveSize() ;
  for (i=1 ; i<=nLoads ; i++) {
    n     = bodyLoadArray.at(i) ;
    load  = (Load*) domain->giveLoad(n) ;
    ltype = load->giveBCGeoType();
    if (ltype == BodyLoadBGT) {
      
      for (igp=0 ; igp<iRule->getNumberOfIntegrationPoints() ; igp++) {
        gp  = iRule->getIntegrationPoint(igp) ;
        this -> computeNSubMatrixAt(nm, gp->giveCoordinates()) ;
        dV  = this -> computeVolumeAround(gp) ;
        this->computeGlobalCoordinates (globalIPcoords, *gp->giveCoordinates());
        load->computeValueAt (val, atTime, globalIPcoords, mode);
        
        nm.times(val.at(indx)*dV);
        if (helpLoadVector.isEmpty()) helpLoadVector.resize(nm.giveNumberOfColumns());
        for (j=1; j<=nm.giveNumberOfColumns(); j++) helpLoadVector.at(j)+=nm.at(1,j) ;
      }
      answer.add(helpLoadVector);
    }
  }

  // add internal source produced by material (if any)
  if (mat->hasInternalSource()) {
    for (igp=0 ; igp<iRule->getNumberOfIntegrationPoints() ; igp++) {
      gp  = iRule->getIntegrationPoint(igp) ;
      this -> computeNSubMatrixAt(nm, gp->giveCoordinates()) ;
      dV  = this -> computeVolumeAround(gp) ;
      mat->computeInternalSourceVector(val, gp, atTime, mode);
      nm.times(val.at(indx)*dV);
      if (helpLoadVector.isEmpty()) helpLoadVector.resize(nm.giveNumberOfColumns());
      for (j=1; j<=nm.giveNumberOfColumns(); j++) helpLoadVector.at(j)+=nm.at(1,j) ;
    }
    answer.add(helpLoadVector);
  }

 return;
}


 
void
TransportElement :: computeIntSourceLHSMatrix (FloatMatrix& answer, TimeStep* tStep)
{
  TransportMaterial *mat = (TransportMaterial*) this->giveMaterial();
  if (mat -> hasInternalSource()) {
    answer.resize(computeNumberOfDofs(EID_ConservationEquation), computeNumberOfDofs(EID_ConservationEquation));
    answer.zero();
    
    if (emode == HeatTransferEM) {
      this->computeIntSourceLHSSubMatrix (answer, IntSource, 0, tStep);
    } else if (emode == HeatMass1TransferEM) {
      FloatMatrix subAnswer;
      int i;
      MatResponseMode rmode[2]={IntSource_hh, IntSource_ww};
      double coeff = 1.0; //this->giveMaterial()->give('d');
      
      for (i=1; i<=2; i++) {
        this->computeIntSourceLHSSubMatrix (subAnswer, rmode[i-1], 0, tStep);
        this->assembleLocalContribution (answer, subAnswer, 2, i, i, coeff);
      }
    } else {
      _error ("Unknown ElementMode");
    }
  } else {
    answer.resize(0,0);
  }
}

void
TransportElement::computeIntSourceLHSSubMatrix (FloatMatrix& answer, MatResponseMode rmode, int iri, TimeStep* tStep)
// computes LHS matrix due to material internal source (dHeat/dT, dWaterSource/dw)
// IntSource - Heat transfer
// IntSource_hh - HeMo heat source
// IntSource_ww - HeMo water content source
// hw, wh is not used now
{
  int         i ;
  double      dV, c ;
  FloatMatrix n ;
  GaussPoint  *gp ;
  IntegrationRule* iRule = integrationRulesArray[iri];

  answer.resize (0,0);
  answer.zero();
  for (i=0 ; i<iRule->getNumberOfIntegrationPoints() ; i++) {
    gp      = iRule->getIntegrationPoint(i) ;
    this -> computeNSubMatrixAt(n, gp->giveCoordinates()) ;
    // ask for coefficient from material
    c = ((TransportMaterial*)this->giveMaterial())->giveCharacteristicValue(rmode,gp,tStep);
    dV      = this -> computeVolumeAround(gp) ;
    answer.plusProductSymmUpper(n, n, dV * c) ;
  }

  answer.symmetrized() ;
}


 
void
TransportElement :: computeConstitutiveMatrixAt (FloatMatrix& answer,
                         MatResponseMode rMode, GaussPoint* gp,
                         TimeStep* tStep)
{
  ((TransportMaterial*)this->giveMaterial())->giveCharacteristicMatrix(answer,FullForm,rMode,gp,tStep);
}

/*
void
TransportElement :: computeDirichletBcRhsVectorAt (FloatArray& answer, TimeStep* stepN, CharTypeMode mode)
   // Computes the load vector due to the Dirichlet boundary conditions acting on the
   // receiver's nodes, at stepN. 
{
 FloatArray  d, dp;
 FloatMatrix s;
 UnknownTypeMode umode;

 if (mode == TotalMode) umode = UnknownMode_Total;
 else if (mode == IncrementalMode) umode = UnknownMode_Incremental;
 else _error ("computeBcRhsVectorAt: unknown mode encountered");

 this -> computeVectorOfPrescribed(umode,stepN, d) ;

 if (d.containsOnlyZeroes())
   answer.resize (0);
 else {
   this -> computeConductivityMatrix(s, Conductivity, stepN);
   answer.beProductOf (s, d); answer.negated() ;
 }
 
 return  ;
}
*/


void 
TransportElement :: computeBCVectorAt (FloatArray& answer, TimeStep* tStep, ValueModeType mode)
{
 answer.resize(computeNumberOfDofs(EID_ConservationEquation));
 answer.zero();

 if (emode == HeatTransferEM) {
   this->computeBCSubVectorAt (answer, tStep, mode, 1);
 } else if (emode == HeatMass1TransferEM) {
   FloatArray subAnswer;
   int i;
   
   for (i=1; i<=2; i++) {
     this->computeBCSubVectorAt (subAnswer, tStep, mode, i);
     this->assembleLocalContribution (answer, subAnswer, 2, i, 1.0);
   }    
 } else {
   _error ("Unknown ElementMode");
 }

}

void 
TransportElement :: computeBCMtrxAt (FloatMatrix& answer, TimeStep* tStep, ValueModeType mode)
{
  int ndofs = computeNumberOfDofs(EID_ConservationEquation);
  answer.resize (ndofs, ndofs);
  answer.zero();

 if (emode == HeatTransferEM) {
   this->computeBCSubMtrxAt (answer, tStep, mode, 1);
 } else if (emode == HeatMass1TransferEM) {
   FloatMatrix subAnswer;
   
   for (int i=1; i<=2; i++) {
     this->computeBCSubMtrxAt (subAnswer, tStep, mode, i);
     if (subAnswer.isNotEmpty()) {
       this->assembleLocalContribution (answer, subAnswer, 2, i, i, 1.0);
     }
   }
 } else {
   _error ("Unknown ElementMode");
 }
}


void 
TransportElement :: computeBCSubVectorAt (FloatArray& answer, TimeStep* tStep, ValueModeType mode, int indx)
{
  int n, id;
  GeneralBoundaryCondition*       load; 
  bcGeomType  ltype;
  FloatArray vec;

  answer.resize (this->giveNumberOfDofManagers());
  answer.zero();

  // loop over boundary load array
  int nLoads    = this -> giveBoundaryLoadArray() -> giveSize() / 2;
  for (int i=1 ; i<=nLoads ; i++) {
    n     = boundaryLoadArray.at(1+(i-1)*2) ;
    id    = boundaryLoadArray.at(i*2) ;
    load  = (GeneralBoundaryCondition*) domain->giveLoad(n) ;
    ltype = load->giveBCGeoType();
    if (ltype == EdgeLoadBGT) {
      this->computeEdgeBCSubVectorAt (vec, (Load*)load, id, tStep, mode, indx);
    } else if (ltype == SurfaceLoadBGT) {
      this->computeSurfaceBCSubVectorAt (vec, (Load*)load, id, tStep, mode, indx);
    } else {
      _error("computeBCSubVectorAt : unsupported bc type encountered");
    }
    answer.add(vec);
  } // end loop over applied bc
}

void 
TransportElement::computeEdgeBCSubVectorAt (FloatArray& answer, Load* load, int iEdge, 
                                            TimeStep* tStep, ValueModeType mode, int indx)
{
  int i, j, approxOrder, numberOfEdgeIPs;
  
  answer.resize (this->giveNumberOfDofManagers());
  answer.zero();

  if (((load->giveType() == TransmissionBC) || (load->giveType() == ConvectionBC))) {
    BoundaryLoad* edgeLoad = static_cast<BoundaryLoad*>(load);
    if (edgeLoad->isDofExcluded(indx)) return;

    
    approxOrder = edgeLoad->giveApproxOrder()+this->giveApproxOrder(indx); 
    numberOfEdgeIPs = (int) ceil ((approxOrder+1.)/2.);
    GaussIntegrationRule iRule(1, this, 1, 1);
    iRule.setUpIntegrationPoints (_Line, numberOfEdgeIPs, _Unknown);
    GaussPoint* gp ;
    FloatArray reducedAnswer, val, ntf;
    IntArray mask;
    FloatMatrix n;
    double dV, coeff = 1.0;
    
    if (load->giveType() == TransmissionBC) coeff = -1.0;
    else coeff = edgeLoad->giveProperty ('a');
    
    for (i=0 ; i < iRule.getNumberOfIntegrationPoints() ; i++) {
      gp  = iRule.getIntegrationPoint(i) ;
      this -> computeEgdeNMatrixAt(n, gp) ;
      dV  = this -> computeEdgeVolumeAround(gp, iEdge) ;
      //nt.beTranspositionOf (n);
      
      if (edgeLoad->giveFormulationType() == BoundaryLoad::BL_EntityFormulation)
        edgeLoad->computeValueAt (val, tStep, *(gp->giveCoordinates()), mode);
      else {
        FloatArray globalIPcoords;
        this->computeEdgeIpGlobalCoords (globalIPcoords, gp, iEdge);
        edgeLoad->computeValueAt (val, tStep, globalIPcoords, mode);
      }
      
      n.times(val.at(indx)*coeff*dV);
      if (reducedAnswer.isEmpty()) reducedAnswer.resize(n.giveNumberOfColumns());
      for (j=1; j<=n.giveNumberOfColumns(); j++) reducedAnswer.at(j)+=n.at(1,j) ;
    }
    this -> giveEdgeDofMapping (mask, iEdge);
    answer.assemble(reducedAnswer, mask);
  } else {
    _error("computeBCSubVectorAt : unsupported bc type encountered");
  }
}

void
TransportElement::computeSurfaceBCSubVectorAt (FloatArray& answer, Load* load, 
                                               int iSurf, TimeStep* tStep, ValueModeType mode, int indx)
{
  int i, j, approxOrder;
  double dV, coeff = 1.0 ;

  if (!this->testElementExtension(Element_SurfaceLoadSupport))
  _error("computeSurfaceBCSubVectorAt : no surface load support");

  BoundaryLoad* surfLoad = dynamic_cast<BoundaryLoad*>(load);
  if (surfLoad) {
    IntegrationRule* iRule;
    GaussPoint* gp ;
    FloatArray reducedAnswer, val, globalIPcoords;
    IntArray mask;
    FloatMatrix  n;
    
    answer.resize (this->giveNumberOfDofManagers());
    answer.zero();
    
    if (surfLoad->isDofExcluded(indx)) return;

    if (load->giveType() == TransmissionBC) coeff = -1.0;
    else coeff = surfLoad->giveProperty ('a');

    approxOrder = surfLoad->giveApproxOrder()+this->giveApproxOrder(indx);
  
    iRule = this -> GetSurfaceIntegrationRule (approxOrder);
    for (i=0 ; i < iRule->getNumberOfIntegrationPoints() ; i++) {
      gp  = iRule-> getIntegrationPoint(i) ;
      this -> computeSurfaceNMatrixAt(n, gp) ;
      dV  = this -> computeSurfaceVolumeAround(gp, iSurf) ;
      //nt.beTranspositionOf ( n );
      
      if (surfLoad->giveFormulationType() == BoundaryLoad::BL_EntityFormulation)
        surfLoad->computeValueAt (val, tStep, *(gp->giveCoordinates()), mode);
      else {
        this->computeSurfIpGlobalCoords (globalIPcoords, gp, iSurf);
        surfLoad->computeValueAt (val, tStep, globalIPcoords, mode);
      }

      n.times(val.at(indx)*coeff*dV);
      if (reducedAnswer.isEmpty()) reducedAnswer.resize(n.giveNumberOfColumns());
      for (j=1; j<=n.giveNumberOfColumns(); j++) reducedAnswer.at(j)+=n.at(1,j) ;
    }
    this -> giveSurfaceDofMapping (mask, iSurf);
    answer.assemble(reducedAnswer, mask);
  
    delete iRule;
    
  } else {
    _error("computeSurfaceBCSubVectorAt : unsupported bc type encountered");
  }
}




void 
TransportElement :: computeBCSubMtrxAt (FloatMatrix& answer, TimeStep* tStep, ValueModeType mode, int indx)
{
 int i, igp, n, id, defined = 0;
 GeneralBoundaryCondition*       load;
 double dV;
 bcGeomType  ltype;

 answer.resize (this->giveNumberOfDofManagers(), this->giveNumberOfDofManagers());
 answer.zero();

  // loop over boundary load array
  int nLoads    = this -> giveBoundaryLoadArray() -> giveSize() / 2;
  for (i=1 ; i<=nLoads ; i++) {
    n     = boundaryLoadArray.at(1+(i-1)*2) ;
    id    = boundaryLoadArray.at(i*2) ;
    load  = (Load*) domain->giveLoad(n) ;
    if ((load->giveType() == ConvectionBC)) {
      ltype = load->giveBCGeoType();
      if (ltype == EdgeLoadBGT) {
        BoundaryLoad* edgeLoad = static_cast<BoundaryLoad*>(load);
        if (edgeLoad->isDofExcluded(indx)) continue;

        defined = 1;
        int approxOrder = 2*this->giveApproxOrder(indx); 
        int numberOfEdgeIPs = (int) ceil ((approxOrder+1.)/2.);
        GaussIntegrationRule iRule(1, this, 1, 1);
        iRule.setUpIntegrationPoints (_Line, numberOfEdgeIPs, _Unknown);
        GaussPoint* gp ;
        FloatArray val;
        IntArray mask;
        FloatMatrix subAnswer, n;
        
        for (igp=0 ; igp < iRule.getNumberOfIntegrationPoints() ; igp++) {
          gp  = iRule.getIntegrationPoint(igp) ;
          this -> computeEgdeNMatrixAt(n, gp) ;
          dV  = this -> computeEdgeVolumeAround(gp, id) ;
          subAnswer.plusProductSymmUpper(n, n, dV*edgeLoad->giveProperty ('a')) ;
        }
        subAnswer.symmetrized();
        this -> giveEdgeDofMapping (mask, id);
        answer.assemble(subAnswer, mask);
      } else if (ltype == SurfaceLoadBGT) {
        IntegrationRule* iRule;
        GaussPoint* gp ;
        FloatArray val;
        IntArray mask;
        FloatMatrix subAnswer, n;
        
        BoundaryLoad* surfLoad = static_cast<BoundaryLoad*>(load);
        if (surfLoad->isDofExcluded(indx)) continue;

        defined = 1;
        int approxOrder = 2*this->giveApproxOrder(indx);
        iRule = this -> GetSurfaceIntegrationRule (approxOrder);
        
        for (igp=0 ; igp < iRule->getNumberOfIntegrationPoints() ; igp++) {
          gp  = iRule->getIntegrationPoint(igp) ;
          this -> computeSurfaceNMatrixAt(n, gp) ;
          dV  = this -> computeSurfaceVolumeAround(gp, id) ;
          subAnswer.plusProductSymmUpper(n, n, dV*surfLoad->giveProperty ('a')) ;
        }
        delete iRule;
        subAnswer.symmetrized();
        this -> giveSurfaceDofMapping (mask, id);
        answer.assemble(subAnswer, mask);
        
      } else _error("computeBCSubMtrxAt : unsupported bc type encountered");
    }
  } // end loop over applied bc
 if (!defined) answer.resize(0,0);
}


void
TransportElement :: assembleLocalContribution (FloatMatrix& answer, FloatMatrix& src, 
                        int ndofs, int rdof, int cdof, double coeff)
{
 int i,j, ti, tj;
 int nnodes = this->giveNumberOfDofManagers();

 for (i=1; i<=nnodes; i++) {
  ti = (i-1)*ndofs+rdof;
  for (j=1; j<=nnodes; j++) {
   tj = (j-1)*ndofs+cdof;
   answer.at(ti,tj)+= src.at(i,j)*coeff;
  }
 }
}


void
TransportElement :: assembleLocalContribution (FloatArray& answer, FloatArray& src, 
                        int ndofs, int rdof, double coeff)
{
 int i, ti;
 int nnodes = this->giveNumberOfDofManagers();

 for (i=1; i<=nnodes; i++) {
  ti = (i-1)*ndofs+rdof;
  answer.at(ti)+= src.at(i)*coeff;
 }
}


void  
TransportElement :: updateInternalState (TimeStep* stepN)
   // Updates the receiver at end of step.
{
  int i,j ;
  IntegrationRule* iRule;
  FloatArray f,r;
  FloatMatrix n;
  TransportMaterial* mat = ((TransportMaterial*)this->giveMaterial());
  GaussPoint* gp;
  
  // force updating ip values
  for (i=0 ; i < numberOfIntegrationRules ; i++) {
    iRule = integrationRulesArray[i];
    for (j=0; j < iRule->getNumberOfIntegrationPoints(); j++) {
      gp = iRule->getIntegrationPoint(j);
      this -> computeNmatrixAt(n, gp->giveCoordinates()) ;
      this -> computeVectorOf (EID_ConservationEquation, VM_Total, stepN, r);
      f.beProductOf(n,r);
      mat->updateInternalState (f, gp, stepN);
    }
  }
}

void 
TransportElement::EIPrimaryFieldI_evaluateFieldVectorAt (FloatArray& answer, PrimaryField& pf, 
                             FloatArray& coords, IntArray& dofId, ValueModeType mode, 
                             TimeStep* atTime)
{
 int i,j,indx;
 double sum;
 FloatArray elemvector, f, lc;
 FloatMatrix n;
 IntArray elemdofs;
 // determine element dof ids
 this->giveElementDofIDMask  (pf.giveEquationID(), elemdofs) ;
 // first evaluate element unknown vector
 this->computeVectorOf (pf, mode, atTime, elemvector);
 // determine corresponding local coordinates
 if (this->computeLocalCoordinates (lc, coords)) {
  // compute interpolation matrix
  this->computeNmatrixAt(n, &lc);
  // compute answer
  answer.resize(dofId.giveSize());
  for (i = 1; i <= dofId.giveSize(); i++) {
   if ((indx = elemdofs.findFirstIndexOf(dofId.at(i)))) {
    for (j=1, sum = 0.0; j<=elemvector.giveSize(); j++)
     sum += n.at(indx,j)*elemvector.at(j);
    answer.at(i) = sum;
   } else _error ("EIPrimaryFieldI_evaluateFieldVectorAt: unknown dof id encountered");
  }  
 } else _error ("EIPrimaryFieldI_evaluateFieldVectorAt: target point not in receiver volume");
}


#ifdef __OOFEG
int 
TransportElement::giveInternalStateAtNode (FloatArray& answer, InternalStateType type, InternalStateMode mode, 
                                           int node, TimeStep* atTime)
{

  Node* n = this->giveNode(node);
  if (type == IST_Temperature) {
    int dofindx;
    if ((dofindx = n->findDofWithDofId (T_f))) {
      answer.resize(1);
      answer.at(1) = n->giveDof(dofindx)->giveUnknown(EID_ConservationEquation, VM_Total, atTime);
      return 1;
    } else return 0;
  } else if (type == IST_MassConcentration_1) {
    int dofindx;
    if ((dofindx = n->findDofWithDofId (C_1))) {
      answer.resize(1);
      answer.at(1) = n->giveDof(dofindx)->giveUnknown(EID_ConservationEquation, VM_Total, atTime);
      return 1;
    } else return 0;
  } else return Element::giveInternalStateAtNode (answer, type, mode, node, atTime);
}

#endif

int 
TransportElement::giveIntVarCompFullIndx (IntArray& answer, InternalStateType type)
{
  if ((type == IST_Temperature) || (type == IST_MassConcentration_1)) {
    answer.resize(1); answer.at(1) = 1; return 1;
  } else {
    return Element::giveIntVarCompFullIndx (answer, type);
  }
}