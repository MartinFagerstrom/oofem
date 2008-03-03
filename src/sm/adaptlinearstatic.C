/* $Header: /home/cvs/bp/oofem/sm/src/adaptlinearstatic.C,v 1.6 2003/04/06 14:08:30 bp Exp $ */
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
// file adaptivelinearstatic.cc
//

#include "verbose.h"
#include "adaptlinearstatic.h"
#include "timestep.h"
#include "nummet.h"

#include "remeshingcrit.h"
#include "t3dinterface.h"
#include "targe2interface.h"
#include "freeminterface.h"

#include "eleminterpunknownmapper.h"
#include "errorestimator.h"
#include "usrdefsub.h"

void
AdaptiveLinearStatic::solveYourselfAt (TimeStep* tStep)
{
 LinearStatic::solveYourselfAt (tStep);

 // perform error evaluation
 // evaluate error of the reached solution
 this->ee->estimateError (temporaryEM, this->giveCurrentStep());
// this->ee->estimateError (equilibratedEM, this->giveCurrentStep());
 RemeshingStrategy strategy = this->ee->giveRemeshingCrit () -> giveRemeshingStrategy (this->giveCurrentStep());
 
 if (strategy == NoRemeshing_RS) return;
 else {
  // do remeshing
  MesherInterface* mesher;
  if (this->meshPackage == MPT_TARGE2) mesher = new Targe2Interface();
  else if (this->meshPackage == MPT_FREEM) mesher = new FreemInterface();
  else mesher = new T3DInterface();

   mesher->createMesh (this->giveDomain(1), this->giveCurrentStep(),1,this->giveDomain(1)->giveSerialNumber()+1);
  // terminate step
  this->terminate(this->giveCurrentStep());
  this->terminateAnalysis();
  exit (1);
 }
}


int  
AdaptiveLinearStatic::initializeAdaptive (int stepNumber)
{
 /*
   Due to linear character of the problem,
   the whole analysis is restarted from beginning.
   The solution steps represent the adaptive steps and for each adaptive step 
   new domain with corresponding domainSerNum is generated.
   */
 int result = 1;
/*
 this -> initStepIncrements();

 int sernum = stepNumber + 1;
 printf ("\nrestoring domain %d.%d\n", 1, sernum);
 Domain* dNew = new Domain (1, sernum, this);
 FILE* domainInputFile;
 this->giveDomainFile (&domainInputFile, 1, sernum, contextMode_read);
 if (!dNew -> instanciateYourself(domainInputFile)) _error ("initializeAdaptive: domain Instanciation failed");
 fclose (domainInputFile);

 printf ("\ndeleting old domain\n");
 delete domainList->at(1);
 domainList->put(1, dNew);

 // init equation numbering
 this->forceEquationNumbering();

 // set time step
 this->giveCurrentStep()->setTime (stepNumber+1);

 // init equation numbering
 // this->forceEquationNumbering();
 this->giveNumericalMethod(giveCurrentStep())->setDomain (dNew);
 this->ee->setDomain (dNew);
*/
 return result;
}
 

contextIOResultType AdaptiveLinearStatic:: restoreContext (DataStream* stream, ContextMode mode, void* obj) {
  contextIOResultType iores;
  if ((iores = LinearStatic::restoreContext (stream, mode, obj)) != CIO_OK) THROW_CIOERR (iores);
  return CIO_OK;
}

IRResultType
AdaptiveLinearStatic :: initializeFrom (InputRecord* ir)
// input from inputString
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro
 ErrorEstimatorType eeType;
 
 LinearStatic::initializeFrom (ir);

 int eeTypeId = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, eeTypeId, IFT_AdaptiveLinearStatic_eetype, "eetype"); // Macro
 eeType = (ErrorEstimatorType) eeTypeId;
 this-> ee = ::CreateUsrDefErrorEstimator (eeType, 1, this->giveDomain(1)); 

 ee->initializeFrom (ir);

 int meshPackageId = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, meshPackageId, IFT_AdaptiveLinearStatic_meshpackage, "meshpackage"); // Macro

 if (meshPackageId == 1) meshPackage = MPT_TARGE2; 
 else if (meshPackageId == 2) meshPackage = MPT_FREEM; 
 else meshPackage = MPT_T3D;

 return IRRT_OK;
}


void
 AdaptiveLinearStatic::updateDomainLinks()
{
 LinearStatic::updateDomainLinks();
 // associate ee to possibly newly restored mesh
 this->ee->setDomain(this->giveDomain(1));
}