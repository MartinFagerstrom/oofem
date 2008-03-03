/* $Header: /home/cvs/bp/oofem/sm/src/isodamagemodel.C,v 1.4.4.1 2004/04/05 15:19:47 bp Exp $ */
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

// file: isodamagemodel.C


#include "isodamagemodel.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "structuralcrosssection.h"
#include "mathfem.h"
#include "datastream.h"

IsotropicDamageMaterial :: IsotropicDamageMaterial (int n, Domain *d) : StructuralMaterial (n,d)
//
// constructor
//
{
 linearElasticMaterial = NULL;
 llcriteria = idm_strainLevelCR;
}


IsotropicDamageMaterial :: ~IsotropicDamageMaterial ()
//
// destructor
//
{
 delete linearElasticMaterial;
}

int
IsotropicDamageMaterial :: hasMaterialModeCapability (MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
 if ((mode == _3dMat) || (mode == _PlaneStress) || 
   (mode == _PlaneStrain) || (mode == _1dMat)) return 1;
 return 0;
}


void
IsotropicDamageMaterial :: give3dMaterialStiffnessMatrix (FloatMatrix& answer, 
                             MatResponseForm form,
                             MatResponseMode mode,
                             GaussPoint* gp,
                             TimeStep* atTime)
//
// computes full constitutive matrix for case of gp stress-strain state.
//
{
 IsotropicDamageMaterialStatus *status = (IsotropicDamageMaterialStatus*) this -> giveStatus (gp);
 double om;
 if (mode == ElasticStiffness) {
  om = 0.0;
 } else {
  om = status->giveTempDamage();
  om = min (om, 0.999999);
 }

 this->giveLinearElasticMaterial()->give3dMaterialStiffnessMatrix (answer, form, mode, gp, atTime);
 answer.times (1.0-om);

 return ;
 
}


void
IsotropicDamageMaterial :: giveRealStressVector (FloatArray& answer, MatResponseForm form, GaussPoint* gp, 
                         const FloatArray& totalStrain, 
                         TimeStep* atTime)
//
// returns real stress vector in 3d stress space of receiver according to 
// previous level of stress and current
// strain increment, the only way, how to correctly update gp records
//
{
 IsotropicDamageMaterialStatus *status = (IsotropicDamageMaterialStatus*) this -> giveStatus (gp);
 //StructuralCrossSection *crossSection = (StructuralCrossSection*) gp -> giveElement()->giveCrossSection();
 LinearElasticMaterial* lmat = this->giveLinearElasticMaterial();
 FloatArray strainVector, reducedTotalStrainVector;
 FloatMatrix de;
 double f, equivStrain, tempKappa=0.0, omega=0.0; 

 this->initGpForNewStep(gp);

 // substract stress independent part
 // note: eigenStrains (tepmerature) is not contained in mechanical strain stored in gp
 // therefore it is necessary to substract always the total eigen strain value
 this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain,atTime, VM_Total);

 //crossSection->giveFullCharacteristicVector(totalStrainVector, gp, reducedTotalStrainVector);
 
 // compute equivalent strain
 this->computeEquivalentStrain (equivStrain, reducedTotalStrainVector, gp, atTime);

 if (llcriteria == idm_strainLevelCR) {
  // compute value of loading function if strainLevel crit apply
  f = equivStrain - status->giveKappa();

  if (f <= 0.0) {
   // damage do not grow
   tempKappa = status->giveKappa();
   omega     = status->giveDamage();
  } else {
   // damage grow
   tempKappa = equivStrain;
   this->initDamaged (tempKappa, reducedTotalStrainVector, gp);
   // evaluate damage parameter
   this->computeDamageParam (omega, tempKappa, reducedTotalStrainVector, gp);
  }
 } else if (llcriteria == idm_damageLevelCR) {

  // evaluate damage parameter first
  tempKappa = equivStrain;
  this->initDamaged (tempKappa, reducedTotalStrainVector, gp);
  this->computeDamageParam (omega, tempKappa, reducedTotalStrainVector, gp);
  if (omega < status->giveDamage()) {
   // unloading takes place
   omega = status->giveDamage();
   //printf (".");
  }
 } else _error ("giveRealStressVector: unsupported loading/uloading criteria");


 lmat->giveCharacteristicMatrix(de, ReducedForm, SecantStiffness, gp, atTime);
 de.times(1.0-omega);
 answer.beProductOf (de, reducedTotalStrainVector);

 // update gp
 status-> letTempStrainVectorBe (totalStrain);
 status-> letTempStressVectorBe (answer);
 status-> setTempKappa (tempKappa);
 status-> setTempDamage(omega);
 return ;
}


void IsotropicDamageMaterial::givePlaneStressStiffMtrx (FloatMatrix& answer, MatResponseForm form, MatResponseMode mode,
                            GaussPoint* gp, TimeStep* atTime)
{
 IsotropicDamageMaterialStatus *status = (IsotropicDamageMaterialStatus*) this -> giveStatus (gp);
 double om;
 if (mode == ElasticStiffness) {
  om = 0.0;
 } else {
  om = status->giveTempDamage();
  om = min (om, 0.999999);
 }

 this->giveLinearElasticMaterial()->giveCharacteristicMatrix (answer, form, mode, gp, atTime);
 answer.times (1.0-om);

 return ;
 
}


void IsotropicDamageMaterial::givePlaneStrainStiffMtrx (FloatMatrix& answer, MatResponseForm form,MatResponseMode mode,
                            GaussPoint* gp, TimeStep* atTime)
{
 IsotropicDamageMaterialStatus *status = (IsotropicDamageMaterialStatus*) this -> giveStatus (gp);
 double om;
 if (mode == ElasticStiffness) {
  om = 0.0;
 } else {
  om = status->giveTempDamage();
  om = min (om, 0.999999);
 }

 this->giveLinearElasticMaterial()->giveCharacteristicMatrix (answer, form, mode, gp, atTime);
 answer.times (1.0-om);

 return ;

}
void IsotropicDamageMaterial::give1dStressStiffMtrx (FloatMatrix& answer, MatResponseForm form,MatResponseMode mode,
                           GaussPoint* gp, TimeStep* atTime)
{
 IsotropicDamageMaterialStatus *status = (IsotropicDamageMaterialStatus*) this -> giveStatus (gp);
 double om;
 if (mode == ElasticStiffness) {
  om = 0.0;
 } else {
  om = status->giveTempDamage();
  om = min (om, 0.999999);
 }

 this->giveLinearElasticMaterial()->giveCharacteristicMatrix (answer, form, mode, gp, atTime);
 answer.times (1.0-om);

 return ;

}
 
#ifdef __OOFEG
#endif



int 
IsotropicDamageMaterial::giveIPValue (FloatArray& answer, GaussPoint* aGaussPoint, InternalStateType type, TimeStep* atTime)
{
 IsotropicDamageMaterialStatus* status = (IsotropicDamageMaterialStatus*) this -> giveStatus (aGaussPoint);
 if ((type == IST_DamageTensor)||(type == IST_PrincipalDamageTensor)) {
   answer.resize(1);
  answer.at(1) = status->giveDamage();
  return 1;
 } else if ((type == IST_DamageTensorTemp)||(type == IST_PrincipalDamageTempTensor)) {
   answer.resize(1);
  answer.at(1) = status->giveTempDamage();
  return 1;
 } else if (type == IST_MaxEquivalentStrainLevel) {
  answer.resize(1);
  answer.at(1) = status->giveKappa ();
  return 1;
 }else return StructuralMaterial::giveIPValue (answer, aGaussPoint, type, atTime);
}
  



InternalStateValueType 
IsotropicDamageMaterial::giveIPValueType (InternalStateType type)
{
 if ((type == IST_DamageTensor)||(type ==IST_DamageTensorTemp) ||
   (type == IST_PrincipalDamageTensor) || (type == IST_PrincipalDamageTempTensor)) return ISVT_TENSOR_S3;
 else if (type == IST_MaxEquivalentStrainLevel) return ISVT_SCALAR;
 else return StructuralMaterial::giveIPValueType (type);
}


int 
IsotropicDamageMaterial::giveIntVarCompFullIndx (IntArray& answer, InternalStateType type, MaterialMode mmode)
{
 if ((type == IST_DamageTensor)||(type ==IST_DamageTensorTemp) ||
   (type == IST_PrincipalDamageTensor) || (type == IST_PrincipalDamageTempTensor)) {
  answer.resize (9);
  answer.at(1) = 1;
  return 1;
 } else if (type == IST_MaxEquivalentStrainLevel) {
  answer.resize (1);
  answer.at(1) = 1;
  return 1;
 } else 
  return StructuralMaterial::giveIntVarCompFullIndx (answer, type, mmode);
}


int
IsotropicDamageMaterial::giveIPValueSize (InternalStateType type, GaussPoint* aGaussPoint)
{
 if ((type == IST_DamageTensor)||(type ==IST_DamageTensorTemp) ||
   (type == IST_PrincipalDamageTensor) || (type == IST_PrincipalDamageTempTensor)||
   (type == IST_MaxEquivalentStrainLevel)) return 1;
 else return StructuralMaterial::giveIPValueSize (type, aGaussPoint);
}


void
IsotropicDamageMaterial :: giveThermalDilatationVector (FloatArray& answer, 
                          GaussPoint * gp,  TimeStep* tStep)
   //
   // returns a FloatArray(6) of initial strain vector
   // eps_0 = {exx_0, eyy_0, ezz_0, gyz_0, gxz_0, gxy_0}^T
   // caused by unit temperature in direction of 
   // gp (element) local axes
   //
{
 //FloatArray *result = new FloatArray (6);
 answer.resize (6); answer.zero();
 answer.at(1) = this->tempDillatCoeff;
 answer.at(2) = this->tempDillatCoeff;
 answer.at(3) = this->tempDillatCoeff;
 
 return ;
}


IRResultType
IsotropicDamageMaterial :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 IR_GIVE_FIELD (ir, tempDillatCoeff, IFT_IsotropicDamageMaterial_talpha, "talpha"); // Macro
 return StructuralMaterial::initializeFrom (ir);
} 


int
IsotropicDamageMaterial::giveInputRecordString(std::string &str, bool keyword)
{
 char buff[1024];

 StructuralMaterial::giveInputRecordString(str, keyword);
 sprintf(buff, " talpha %e", this -> tempDillatCoeff);
 str += buff;

 return 1;
}



IsotropicDamageMaterialStatus::IsotropicDamageMaterialStatus (int n, Domain*d, GaussPoint* g):StructuralMaterialStatus(n,d,g)
{
 kappa = tempKappa = 0.0;
 damage = tempDamage = 0.0;
}


IsotropicDamageMaterialStatus::~IsotropicDamageMaterialStatus () 
{}


void 
IsotropicDamageMaterialStatus :: printOutputAt  (FILE *file, TimeStep* tStep)
{
 
 StructuralMaterialStatus :: printOutputAt (file, tStep);
 fprintf (file,"status { ");
 if (this->damage > 0.0) {
  fprintf (file,"kappa %f, damage %f ",this->kappa, this->damage);
 }
 
 fprintf (file,"}\n");
}
 

void 
IsotropicDamageMaterialStatus::initTempStatus ()
{
 StructuralMaterialStatus :: initTempStatus();
 this->tempKappa = this->kappa;
 this->tempDamage= this->damage;
}

void 
IsotropicDamageMaterialStatus::updateYourself(TimeStep* atTime)
{
 StructuralMaterialStatus::updateYourself(atTime);
 this->kappa = this->tempKappa;
 this->damage= this->tempDamage;
}


contextIOResultType
IsotropicDamageMaterialStatus::saveContext (DataStream* stream, ContextMode mode, void *obj)
{
 contextIOResultType iores;

 // save parent class status
 if ((iores = StructuralMaterialStatus :: saveContext (stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);

 // write a raw data
 if (!stream->write(&kappa,1)) THROW_CIOERR(CIO_IOERR);
 if (!stream->write(&damage,1)) THROW_CIOERR(CIO_IOERR);

 return CIO_OK;
}

contextIOResultType
IsotropicDamageMaterialStatus::restoreContext(DataStream* stream, ContextMode mode, void *obj)
{
 contextIOResultType iores;

 // read parent class status
 if ((iores = StructuralMaterialStatus :: restoreContext (stream,mode,obj)) != CIO_OK) THROW_CIOERR(iores);

 // read raw data 
 if (!stream->read (&kappa,1)) THROW_CIOERR(CIO_IOERR);
 if (!stream->read (&damage,1)) THROW_CIOERR(CIO_IOERR);

 return CIO_OK;
}

