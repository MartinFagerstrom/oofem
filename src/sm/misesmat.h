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

#ifndef misesmat_h
#define misesmat_h

#include "structuralmaterial.h"
#include "structuralms.h"
#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
class GaussPoint;
class Domain;

/**
 * This class implements an isotropic elastoplastic material
 * with Mises yield condition, associated flow rule
 * and linear isotropic hardening.
 *
 * It differs from other similar materials (such as J2Mat)
 * by implementation - here we use the radial return, which
 * is the most efficient algorithm for this specific model.
 * Also, an extension to large strain will be available.
 */
class MisesMat : public StructuralMaterial
{
protected:
    /// Reference to the basic elastic material.
    LinearElasticMaterial *linearElasticMaterial;

    /// Elastic shear modulus.
    double G;

    /// Elastic bulk modulus.
    double K;

    /// Hardening modulus.
    double H;

    /// Initial (uniaxial) yield stress.
    double sig0;

    double omega_crit;
    double a;

public:
    MisesMat(int n, Domain *d);
    virtual ~MisesMat();

    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, MaterialMode mode);
    double computeDamage(GaussPoint *gp, TimeStep *atTime);
    double computeDamageParam(double tempKappa);
    double computeDamageParamPrime(double tempKappa);
    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *atTime);

    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual const char *giveClassName() const { return "MisesMat"; }
    virtual classType giveClassID() const { return MisesMatClass; }

    /// Returns a reference to the basic elastic material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    //virtual int giveSizeOfFullHardeningVarsVector();
    //virtual int giveSizeOfReducedHardeningVarsVector(GaussPoint *gp);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                       MatResponseForm form, MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep);

protected:
    /// Evaluates the stress from Green-Lagrange strain E.
    void giveRealStressVectorComputedFromStrain(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                                const FloatArray &E, TimeStep *tStep);

    /// evaluates the stress from deformation gradient F.
    void giveRealStressVectorComputedFromDefGrad(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                                 const FloatArray & F, TimeStep *tStep);

    /// Converts the deformation gradient F into the Green-Lagrange strain E
    void convertDefGradToGLStrain(const FloatMatrix &F, FloatMatrix &E);
    void computeGLPlasticStrain(const FloatMatrix &F, FloatMatrix &Ep, FloatMatrix b, double J);

    void give3dSSMaterialStiffnessMatrix(FloatMatrix &answer,
                                         MatResponseForm form, MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *tStep);
    virtual void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void give3dLSMaterialStiffnessMatrix(FloatMatrix &answer,
                                         MatResponseForm form, MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);

    virtual InternalStateValueType giveIPValueType(InternalStateType type);

    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
};

//=============================================================================


class MisesMatStatus : public StructuralMaterialStatus
{
protected:
    /// Plastic strain (initial).
    FloatArray plasticStrain;

    /// Plastic strain (final).
    FloatArray tempPlasticStrain;

    /// Deviatoric trial stress - needed for tangent stiffness.
    FloatArray trialStressD;

    /**************************************************/
    double trialStressV;
    /**************************************************/

    FloatArray effStress;
    FloatArray tempEffStress;

    /// Cumulative plastic strain (initial).
    double kappa;

    /// Cumulative plastic strain (final).
    double tempKappa;

    /// Deformation gradient(final).
    FloatMatrix tempDefGrad, defGrad;

    /************************/
    double tempDamage, damage;
    /******************************/

    /// Left Cauchy-Green deformation gradient (initial).
    FloatMatrix leftCauchyGreen;
    /// Left Cauchy-Green deformation gradient (final).
    FloatMatrix tempLeftCauchyGreen;

public:
    MisesMatStatus(int n, Domain *d, GaussPoint *g);
    virtual ~MisesMatStatus();

    void givePlasticStrain(FloatArray &answer) { answer = plasticStrain; }

    void giveTrialStressDev(FloatArray &answer) { answer = trialStressD; }

    /*******************************************/
    void giveTrialStressVol(double &answer) { answer = trialStressV; }
    /*******************************************/
    double giveDamage() { return damage; }
    double giveTempDamage() { return tempDamage; }

    double giveCumulativePlasticStrain() { return kappa; }
    double giveTempCumulativePlasticStrain() { return tempKappa; }

    void giveTempDefGrad(FloatMatrix &answer) { answer = tempDefGrad; }
    void giveDefGrad(FloatMatrix &answer) { answer = defGrad; }

    void giveTempLeftCauchyGreen(FloatMatrix &answer) { answer = tempLeftCauchyGreen; }
    void giveLeftCauchyGreen(FloatMatrix &answer) { answer = leftCauchyGreen; }

    void giveTempEffectiveStress(FloatArray &answer) { answer = tempEffStress; }
    void giveEffectiveStress(FloatArray &answer) { answer = effStress; }

    void letTempPlasticStrainBe(FloatArray &values) { tempPlasticStrain = values; }

    void letTrialStressDevBe(FloatArray &values) { trialStressD = values; }

    void letEffectiveStressBe(FloatArray &values) { effStress = values; }

    void letTempEffectiveStressBe(FloatArray &values) { tempEffStress = values; }


    void setTrialStressVol(double value) { trialStressV = value; }

    void setTempCumulativePlasticStrain(double value) { tempKappa = value; }
    /****************************************/
    void setTempDamage(double value) { tempDamage = value; }
    /************************************************/
    void letDefGradBe(FloatMatrix values) { defGrad = values; }
    void letTempDefGradBe(FloatMatrix values) { tempDefGrad = values; }

    void letTempLeftCauchyGreenBe(FloatMatrix values) { tempLeftCauchyGreen = values; }
    void letLeftCauchyGreenBe(FloatMatrix values) { leftCauchyGreen = values; }

    const FloatArray *givePlasDef() { return & plasticStrain; }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "MisesMatStatus"; }
    virtual classType giveClassID() const { return MisesMatStatusClass; }
};
} // end namespace oofem
#endif // misesmat_h
