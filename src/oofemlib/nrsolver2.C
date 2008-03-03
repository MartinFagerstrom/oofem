/* $Header: /home/cvs/bp/oofem/oofemlib/src/nrsolver2.C,v 1.8.4.1 2004/04/05 15:19:43 bp Exp $ */
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
// file nrsolver.C
//

#include "nrsolver2.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <math.h>
#endif
#include "cltypes.h"
#include "verbose.h"
#include "ldltfact.h"
#include "imlsolver.h"
#include "timestep.h"
#include "flotmtrx.h"
//#include "nlinearstatic.h"
#include "mathfem.h"
#include "usrdefsub.h"


#define nrsolver_SMALL_NUM 1.e-20
#define NRSOLVER_MAX_REL_ERROR_BOUND 1.e10
#define NRSOLVER_MAX_RESTARTS 4
#define NRSOLVER_RESET_STEP_REDUCE 0.25
#define NRSOLVER_DEFAULT_NRM_TICKS 10


NRSolver2 ::NRSolver2 (int i, Domain* d,EngngModel* m, EquationID ut) : 
SparseNonLinearSystemNM (i,d,m, ut) {
// 
// constructor
//
nsmax  = 60 ;        // default maximum number of sweeps allowed
rtol   = 10.E-3  ;   // convergence tolerance
//Psi    = 0.1;       // displacement control on
solved = 0 ;
NR_Mode = NR_OldMode = nrsolverModifiedNRM;
NR_ModeTick = -1; // do not swith to calm_NR_OldMode
MANRMSteps = 0;

linSolver = NULL;
linesearchSolver = NULL;
lsFlag = 0;    // no line-search
}

NRSolver2 ::  ~NRSolver2 () {
//
// destructor
//
 if (linSolver) delete linSolver;
 if (linesearchSolver) delete linesearchSolver;
}


NM_Status 
NRSolver2::solve (SparseMtrx* k, FloatArray* R, FloatArray* R0,
         FloatArray* Rr, FloatArray* r, FloatArray* DeltaR, FloatArray* F,
         double& l, double rtol, referenceLoadInputModeType rlm,
         int& nite, TimeStep* tNow)
//
// this function solve the problem of the unbalanced equilibrium 
// using NR scheme
//
//
{
 FloatArray rhs, deltaR, RT;
 FloatArray rInitial;
 //FloatArray F;
 double RRT, forceErr, dispErr = 0.;
 double drr;
 int neq = R->giveSize() ;
 int irest = 0;
 NM_Status status;

 OOFEM_LOG_INFO("Time       Iteration       ForceError      DisplError\n__________________________________________________________\n");

 l = 1.0;
 rInitial = *r;


 status = NM_None;
 this->giveLinearSolver();

 // compute total load R = R+R0
 RT = *R;
 if (R0) RT.add(R0);

 restart:
 DeltaR -> zero();

 //linSolver -> setSparseMtrxAsComponent (LinearEquationLhs,k);

 deltaL = tNow->giveTimeIncrement();

 deltaR.resize(neq);
// if (tNow ->giveNumber() == 1) {
 rhs =  *R;
 // if (R0) rhs.add(*R0);

 //engngModel->updateComponent (tNow, NonLinearRhs_Total);

 RRT = dotProduct(RT.givePointer(),RT.givePointer(),neq);

 //if (R0) RR0 = dotProduct(R0->givePointer(),R0->givePointer(),neq);
 //else RR0 = 0.0;

 nite = 0;

 do {
  nite ++;
  
  if (nite > 1) {
   if ((NR_Mode == nrsolverFullNRM) || ((NR_Mode == nrsolverAccelNRM) && (nite%MANRMSteps == 0))) {
    engngModel->updateComponent (tNow, NonLinearLhs, domain);
    //linSolver -> setSparseMtrxAsComponent (LinearEquationLhs,k);
   }
  }
/*
  linSolver -> setFloatArrayAsComponent (LinearEquationRhs,&rhs);
  linSolver -> setFloatArrayAsComponent (LinearEquationSolution,&deltaR);
  linSolver -> solveYourselfAt (tNow);
  linSolver -> updateYourselfExceptLhs ();
*/
  if ((nite == 1) && (Rr->giveSize())) {
   rhs.add (*Rr);
  }

  linSolver -> solve (k, &rhs, &deltaR);
  //
  // update solution
  //
  if (this->lsFlag && (nite != 1)) {
   // linesearch 
   LineSearchNM::LS_status status;
   IntArray prescribedEqs(0);
   double eta;

   this->giveLineSearchSolver()->solve (r, &deltaR, F, R, R0, prescribedEqs, 1.0, eta, status, tNow);
   DeltaR -> add (deltaR);

  } else {
   r -> add(deltaR);
   DeltaR -> add (deltaR);
   tNow->incrementStateCounter();              // update solution state counter
   //
   // convergency check
   //
   //((NonLinearStatic *)engngModel) -> giveInternalForces(F, *DeltaR, tNow);
   engngModel->updateComponent (tNow, InternalRhs, domain);
   //F->negated();
  }

  rhs = RT;
  //if (R0) rhs.add(*R0);
  rhs.substract(F);
  
  //
  // compute forceError
  //
    // err is relative error of unbalanced forces
  forceErr = dotProduct (rhs.givePointer(),rhs.givePointer(),neq);
  // we compute a relative error norm 
  if ((RRT) > nrsolver_SMALL_NUM) forceErr = sqrt (forceErr / (RRT));
  else forceErr = sqrt (forceErr); // absolute norm
  //
  // compute displacement error
  // 
  // err is relative displacement change
  drr = dotProduct (r->givePointer(), r->givePointer(),neq);
  if (drr < nrsolver_SMALL_NUM) dispErr = 1.;
  else {
   dispErr = dotProduct (deltaR.givePointer(),deltaR.givePointer(),neq) / drr;
   dispErr = sqrt(dispErr);
  }
  // 
  // Restart if nite >= nsmax of if force or displacement error is bigger 
  // than allowed limit (rtol * CALM_MAX_REL_ERROR_BOUND)
  //
  if ((nite >= nsmax) ||
    (fabs(forceErr) > rtol * NRSOLVER_MAX_REL_ERROR_BOUND) ||
    (fabs(dispErr)  > rtol * NRSOLVER_MAX_REL_ERROR_BOUND))
   {

    irest ++;
    if (irest <= NRSOLVER_MAX_RESTARTS) {
     // convergence problems
     // there must be step restart followed by decrease of step length
     // status |= NM_ForceRestart;    
     // reduce step length

     /*
       double time;
       time = tNow->giveTime() - tNow->giveTimeIncrement()*(1.0-NRSOLVER_RESET_STEP_REDUCE) ;
       deltaL =  deltaL * NRSOLVER_RESET_STEP_REDUCE ;
       if (deltaL < minStepLength)  deltaL = minStepLength;
       
       tNow -> setTime(time);
       tNow -> setTimeIncrement(tNow->giveTimeIncrement()*NRSOLVER_RESET_STEP_REDUCE);
       tNow->incrementStateCounter();              // update solution state counter
       */
     
     // restore previous total displacement vector
     r -> times(0.);
     r -> add (rInitial);
     // reset all changes fro previous equilibrium state
     engngModel -> initStepIncrements();
     DeltaR -> zero();
     // restore initial stiffness
     engngModel->updateComponent (tNow, NonLinearLhs,domain);
     // recalculate new Load Vector R
     engngModel->updateComponent (tNow, NonLinearRhs_Incremental, domain);
     //delete F; F = NULL;
#ifdef VERBOSE
     OOFEM_LOG_INFO("NRSolver2 iteration Reset ...\n");
#endif
     NR_OldMode  = NR_Mode;
     NR_Mode     = nrsolverFullNRM;
     NR_ModeTick = NRSOLVER_DEFAULT_NRM_TICKS;
     goto restart;
    } else {
     status = NM_NoSuccess;
     _warning2 ("NRSolver2 - convergence not reached after %d iterations", nsmax);
     // exit(1);
     break;
    }
   }
  OOFEM_LOG_INFO("%-10d %-15d %-15e %-15e\n",(int)tNow->giveTime(),nite,forceErr,dispErr);

 } while ((fabs(forceErr) > rtol) || (fabs(dispErr) > rtol));
 
 //delete F;
 //
 // end of iteration
 //
 // ls ->letSolutionBe(deltar);
 // Lambda += DeltaLambda ;      // *
 //
 // update dofs,nodes,Elemms and print result
 //
#ifdef VERBOSE
 // printf ("\nCALM - step iteration finished") ;
#endif

 status |= NM_Success;
 solved = 1;
 return status;
}    

IRResultType
NRSolver2 :: initializeFrom (InputRecord* ir)
//
// 
//
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 nsmax = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, nsmax, IFT_NRSolver_maxiter, "maxiter"); // Macro
  if (nsmax < 30) nsmax = 30;

  minStepLength = 0.0;
 IR_GIVE_OPTIONAL_FIELD (ir, minStepLength, IFT_NRSolver_minsteplength, "minsteplength"); // Macro

  // read if MANRM method is used
 MANRMSteps = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, MANRMSteps, IFT_NRSolver_manrmsteps, "manrmsteps"); // Macro
  if (MANRMSteps > 0) {
   NR_Mode = NR_OldMode = nrsolverAccelNRM;
  } else {
   NR_Mode = nrsolverModifiedNRM;
  }

 int _val = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, _val, IFT_NRSolver_lstype, "lstype"); // Macro
 solverType = (LinSystSolverType) _val;
 this->giveLinearSolver ()-> initializeFrom (ir);

 this->lsFlag = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, lsFlag, IFT_NRSolver_linesearch, "linesearch"); // Macro

 if (this->lsFlag) this->giveLineSearchSolver() -> initializeFrom (ir);

  return IRRT_OK;
}

contextIOResultType    
NRSolver2 :: saveContext (DataStream* stream, ContextMode mode, void *obj) {
 return CIO_OK;
}

contextIOResultType    
NRSolver2 :: restoreContext(DataStream* stream, ContextMode mode, void *obj) {
 return CIO_OK;
}


SparseLinearSystemNM*
NRSolver2 :: giveLinearSolver() {
  
  if (linSolver) {
    if (linSolver->giveLinSystSolverType()==solverType) {
      return linSolver;
    } else {
      delete linSolver;
    }
  }
  
  linSolver = :: CreateUsrDefSparseLinSolver (solverType, 1, domain, engngModel);
  if (linSolver==NULL) _error ("giveLinearSolver: linear solver creation failed");
  return linSolver;
}

LineSearchNM* 
NRSolver2 :: giveLineSearchSolver()
{
 if (linesearchSolver == NULL) 
  linesearchSolver = new LineSearchNM (1, this->giveDomain(), engngModel);

 return linesearchSolver;
}