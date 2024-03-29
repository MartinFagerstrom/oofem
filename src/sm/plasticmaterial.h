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

#ifndef plasticmaterial_h
#define plasticmaterial_h

#include "structuralmaterial.h"
#include "linearelasticmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "structuralms.h"

namespace oofem {
class GaussPoint;

enum state_flag_values { PM_Elastic, PM_Yielding, PM_Unloading };
/**
 * This class implements associated Material Status to PlasticMaterial.
 */
class PlasticMaterialStatus : public StructuralMaterialStatus
{
protected:
    /// Plastic strain vector (reduced form).
    FloatArray plasticStrainVector;
    FloatArray tempPlasticStrainVector;

    /// Strain space hardening variables.
    FloatArray strainSpaceHardeningVarsVector;
    FloatArray tempStrainSpaceHardeningVarsVector;

    /// Yield function status indicator.
    int state_flag;
    int temp_state_flag;

    /// Plastic consistency parameter.
    double gamma, temp_gamma;

public:
    PlasticMaterialStatus(int n, Domain *d, GaussPoint *g);
    virtual ~PlasticMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    void givePlasticStrainVector(FloatArray &answer) const { answer = plasticStrainVector; }
    void giveTempPlasticStrainVector(FloatArray &answer) const { answer = tempPlasticStrainVector; }
    void giveStrainSpaceHardeningVars(FloatArray &answer) const
    { answer = strainSpaceHardeningVarsVector; }
    void givetempStrainSpaceHardeningVarsVector(FloatArray &answer) const
    { answer = tempStrainSpaceHardeningVarsVector; }

    void letPlasticStrainVectorBe(const FloatArray &v) { plasticStrainVector = v; }
    void letTempPlasticStrainVectorBe(const FloatArray &v) { tempPlasticStrainVector = v; }
    void letTempStrainSpaceHardeningVarsVectorBe(const FloatArray &v) { tempStrainSpaceHardeningVarsVector = v; }
    void letStrainSpaceHardeningVarsVectorBe(const FloatArray &v) { strainSpaceHardeningVarsVector = v; }

    int giveStateFlag() { return state_flag; }
    int giveTempStateFlag() { return temp_state_flag; }
    double givePlasticConsistencyPrameter() { return gamma; }
    double giveTempPlasticConsistencyPrameter() { return temp_gamma; }
    void letTempPlasticConsistencyPrameterBe(double v) { gamma = v; }

    void letTempStateFlagBe(int v) { temp_state_flag = v; }

    // definition
    virtual const char *giveClassName() const { return "PlasticMaterialStatus"; }
    virtual classType giveClassID() const { return PerfectlyPlasticMaterialStatusClass; }
};

/**
 * This class implements a general plastic material.
 * It is assumed to be a base class for many material
 * models based on different yield conditions and hardening laws.
 * Stress return mapping algorithm  is based on general
 * return mapping algorithm, with following assumptions
 * - associative flow rule
 * - general hardening law is supported
 * - Kunt-Tucker conditions apply.
 */
class PlasticMaterial : public StructuralMaterial
{
protected:
    /// Reference to bulk (undamaged) material
    LinearElasticMaterial *linearElasticMaterial;

public:
    PlasticMaterial(int n, Domain *d);
    virtual ~PlasticMaterial();

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour() { return 1; }
    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual const char *giveClassName() const { return "PlasticMaterial"; }
    virtual classType giveClassID() const { return PerfectlyPlasticMaterialClass; }

    /// Returns reference to undamaged (bulk) material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);


    virtual void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);

    virtual InternalStateValueType giveIPValueType(InternalStateType type);

    virtual int giveIPValueSize(InternalStateType type, GaussPoint *gp);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    // add here some auxiliary functions if needed
    FloatArray *ComputeGradientVector(GaussPoint *gp, FloatArray *fullStressVector,
                                      FloatArray *fullStressSpaceHardeningVars);
    FloatArray *ComputeResidualVector(GaussPoint *gp, double Gamma,
                                      FloatArray *plasticStrainVectorR,
                                      FloatArray *strainSpaceHardeningVariables,
                                      FloatArray *gradientVectorR);
    virtual void giveConsistentStiffnessMatrix(FloatMatrix & answer,
                                               MatResponseForm,
                                               MatResponseMode,
                                               GaussPoint * gp,
                                               TimeStep * tStep);

    void  computeConsistentModuli(FloatMatrix &answer,
                                  GaussPoint *gp, FloatMatrix &elasticModuliInverse,
                                  FloatMatrix &hardeningModuliInverse,
                                  double Gamma, const FloatArray &fullStressVector,
                                  const FloatArray &fullStressSpaceHardeningVars);
    void  computeDiagModuli(FloatMatrix &answer,
                            GaussPoint *gp, FloatMatrix &elasticModuliInverse,
                            FloatMatrix &hardeningModuliInverse);

    virtual FloatArray *ComputeStressSpaceHardeningVars(GaussPoint *gp,
                                                        FloatArray *strainSpaceHardeningVariables)
    { return NULL; }
    virtual double computeYieldValueAt(GaussPoint *gp, FloatArray *stressVector,
                                       FloatArray *stressSpaceHardeningVars) { return 0.; }
    virtual void   computeHardeningReducedModuli(FloatMatrix &answer,
                                                 GaussPoint *gp,
                                                 FloatArray *strainSpaceHardeningVariables,
                                                 TimeStep *tStep) = 0;
    virtual FloatArray *ComputeStressGradient(GaussPoint *gp, FloatArray *stressVector,
                                              FloatArray *stressSpaceHardeningVars) { return NULL; }
    virtual FloatArray *ComputeStressSpaceHardeningVarsReducedGradient(GaussPoint *gp,
                                                                       FloatArray *stressVector,
                                                                       FloatArray *stressSpaceHardeningVars)
    { return NULL; }
    virtual int hasHardening() { return 0; }
    virtual void  computeReducedGradientMatrix(FloatMatrix &answer,
                                               GaussPoint *gp,
                                               const FloatArray &stressVector,
                                               const FloatArray &stressSpaceHardeningVars) = 0;

    virtual void computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                             const FloatArray &strainIncrement, TimeStep *tStep);
    virtual void computeReducedElasticModuli(FloatMatrix &answer, GaussPoint *gp,
                                             TimeStep *tStep);
    virtual void compute3dElasticModuli(FloatMatrix &answer, GaussPoint *gp,
                                        TimeStep *tStep) = 0;

    // auxiliary functions
    virtual int         giveSizeOfFullHardeningVarsVector()  { return 0; }
    virtual int         giveSizeOfReducedHardeningVarsVector(GaussPoint *)  { return 0; }

    friend class PlasticMaterialStatus;



    // next functions overloaded rom structural material level
    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer,
                                          MatResponseForm form, MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                          MatResponseForm form, MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void give1dStressStiffMtrx(FloatMatrix &answer,
                                       MatResponseForm form, MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *tStep);

    virtual void give2dBeamLayerStiffMtrx(FloatMatrix &answer,
                                          MatResponseForm form, MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void give2dPlateLayerStiffMtrx(FloatMatrix &answer,
                                           MatResponseForm form, MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep);

    virtual void give1dFiberStiffMtrx(FloatMatrix &answer,
                                      MatResponseForm form, MatResponseMode mode,
                                      GaussPoint *gp,
                                      TimeStep *atTime);

    virtual void give3dShellLayerStiffMtrx(FloatMatrix &answer,
                                           MatResponseForm form, MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep);
};
} // end namespace oofem
#endif // plasticmaterial_h
