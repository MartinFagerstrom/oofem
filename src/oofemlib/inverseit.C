/* $Header: /home/cvs/bp/oofem/oofemlib/src/subspaceit.C,v 1.3 2003/04/06 14:08:26 bp Exp $ */
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

//
// file subspaceit.cc
//

#define DETAILED_REPORT

#include "inverseit.h"
#include "engngm.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <math.h>
#endif
#include "cltypes.h"
#include "verbose.h"
#include "flotmtrx.h"
#include "skyline.h"
#include "flotarry.h"
#include "gjacobi.h"
#include "mathfem.h"

InverseIteration ::  InverseIteration (int i, Domain* d,EngngModel* m) : 
SparseGeneralEigenValueSystemNM (i,d,m) {

  nc     = 0 ;
  rtol   = 10.E-6 ;   // convergence tolerance
  n = 0 ;
  nitem = 100;  // max number of iterations
  solved = 0 ;
}


InverseIteration ::  ~InverseIteration () {}

NM_Status 
InverseIteration::solve (SparseMtrx* a, SparseMtrx* b, FloatArray* _eigv, FloatMatrix* _r, double rtol, int nroot)
//void  SubspaceIteration :: solveYourselfAt (TimeStep* tNow)
//
//
//
{
  int nn,nc,i,j,k,ii, ac;
  double c;

  //
  // check matrix size
  // 
  nc = min(2*nroot,nroot+8) ;
  if ((!a) || (!b)) _error ("SubspaceIteration :: solveYourselfAt : matrices are not defined\n");
  if (a->giveNumberOfColumns() != b->giveNumberOfColumns())
    _error("SubspaceIteration :: solveYourselfAt : matrices size mismatch\n");
  if (!a->canBeFactorized()) _error ("SubspaceIteration :: a matrix not support factorization");
  //
  // check for wery small problem
  //
  nn = a->giveNumberOfColumns() ;
  if (nc > nn ) nc = nn;
  
  FloatArray w(nc), ww(nc), t(nn), tt(nn), *ptr, *ptr2;
  FloatArray z[nc], zz[nc], x[nc];
  
  for (j=0; j<nc; j++) {
    z[j].resize(nn);
    zz[j].resize(nn);
    x[j].resize(nn);
  }

  /*  initial setting  */
  for (i=1;i<=nc;i++){
    ww.at(i)=1.0;
  }
  
  for (j=0; j<nc;j++) {
    ptr = &z[j];
    for (i=1;i<=nn;i++)
      ptr->at(i)=1.0;
  }
  /*  stiffness factorization  */
  a-> factorized () ;
  
  for (i=0;i<nitem;i++) {

    /*  copy zz=z  */
    for (j=0; j<nc; j++) zz[j]=z[j];

    /*  solve matrix equation K.X = M.X  */
    for (j=0;j<nc;j++) {
      x[j] = z[j];
      a->backSubstitutionWith(x[j]) ;
    }

    /*  evaluation of Rayleigh quotients  */
    for (j=0;j<nc;j++){
      w.at(j+1) = dotProduct (zz[j],x[j],nn);
    }

    //mmc_sky (m,x,z,adrm,n,niv);
    for (j=0;j<nc;j++) {
      b->times(x[j], z[j]) ;
    }
    
    for (j=0;j<nc;j++){
      c= dotProduct(z[j], x[j], nn);
      w.at(j+1) /= c;
    }
    
    /*  check convergence  */
    ac=0;
    for (j=1;j<=nc;j++){
      if (fabs((ww.at(j)-w.at(j))/w.at(j))< rtol)  ac++;
      ww.at(j)=w.at(j);
    }
    
    //printf ("\n iteration  %d   %d",i,ac);
    //w.printYourself();
    
    /*  Gramm-Schmidt ortogonalization   */
    for (j=0;j<nc;j++) {
      if (j!=0) {
        b->times(x[j], t);
      }
      for (ii=0;ii<j;ii++) {
        c = dotProduct (x[ii], t, nn);
        ptr = &x[j];
        ptr2= &x[ii];
        for (k = 1; k<= nn; k++) ptr->at(k) -= ptr2->at(k)*c;
      }
      b->times(x[j],t);
      c =   dotProduct (x[j], t, nn);
      x[j].times(1.0/sqrt(c));
    }
    
    if (ac>nroot){
      break;
    }
    
    /*  compute new approximation of Z  */
    for (j=0;j<nc;j++) {
      b->times(x[j],z[j]);
    }
  }

  // copy results
  for (i = 1; i<= nroot; i++) {
    _eigv->at(i) = w.at(i);
    ptr = &x[i-1];
    for (j=1; j<= nn; j++) _r->at(j,i) = ptr->at(j);
  }
  
  if (i < nitem) 
    OOFEM_LOG_INFO ("InverseIt info: convergence reached in %d iterations\n",i);
  else
    OOFEM_LOG_WARNING2 ("InverseIt info: convergence not reached after %d iterations\n",i);
 
  solved = 1;
  return NM_Success;;
}

IRResultType
InverseIteration :: initializeFrom (InputRecord* ir)
//
// 
//
{
  return IRRT_OK;
}


