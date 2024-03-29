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
#ifndef druckerpragerplasticitysm_h
#define druckerpragerplasticitysm_h

#include "flotarry.h"
#include "flotmtrx.h"

#include "structuralms.h"
#include "strainvector.h"
#include "stressvector.h"
#include "structuralmaterial.h"
#include "isolinearelasticmaterial.h"

namespace oofem {
/**
 * This class implements the material status associated to DruckerPragerPlasticitySM.
 * Tracks volumetric and deviatoric plastic strain and hardening.
 * @author Simon Rolshoven
 */
class DruckerPragerPlasticitySMStatus : public StructuralMaterialStatus
{
public:
    /// Values of history variable state_flag.
    enum state_flag_values { DP_Elastic, DP_Unloading, DP_Yielding, DP_Vertex };

protected:
    /// Volumetric plastic strain.
    double volumetricPlasticStrain;
    double tempVolumetricPlasticStrain;

    /// Deviatoric of plastic strain.
    StrainVector plasticStrainDeviator;
    StrainVector tempPlasticStrainDeviator;

    /// Hardening variable.
    double kappa;
    double tempKappa;

    /// Indicates the state (i.e. elastic, yielding, vertex, unloading) of the Gauss point
    int state_flag;
    int temp_state_flag;

public:
    /// Constructor
    DruckerPragerPlasticitySMStatus(int n, Domain *d, GaussPoint *gp);

    /// Destructor
    virtual ~DruckerPragerPlasticitySMStatus();

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "DruckerPragerPlasticitySMStatus"; }
    virtual classType giveClassID() const { return DruckerPragerPlasticitySMStatusClass; }

    /**
     * Get the full plastic strain vector from the material status.
     * @param answer Plastic strain vector.
     */
    void  giveFullPlasticStrainVector(StrainVector &answer) const
    {
        StrainVector plasticStrainVector(_3dMat);
        plasticStrainDeviator.computeDeviatoricVolumetricSum(plasticStrainVector,
                                                             volumetricPlasticStrain);
        plasticStrainVector.convertToFullForm(answer);
    }
    /**
     * Get the plastic strain deviator from the material status.
     * @param answer Plastic strain deviator.
     */
    void givePlasticStrainDeviator(StrainVector &answer) const { answer = plasticStrainDeviator; }
    /**
     * Get the volumetric plastic strain from the material status.
     * @return Volumetric plastic strain.
     */
    double giveVolumetricPlasticStrain() const { return volumetricPlasticStrain; }
    /**
     * Get the hardening variable from the material status.
     * @return hardening variable kappa
     */
    double giveKappa() const { return kappa; }
    /**
     * Get the state flag from the material status.
     * @return State flag (i.e. elastic, unloading, yielding, vertex case yielding)
     */
    int giveStateFlag() const { return state_flag; }

    /**
     * Get the temp value of the full plastic strain vector from the material status.
     * @param answer Temp value of plastic strain vector.
     */
    void giveTempPlasticStrainVector(StrainVector &answer) const
    {
        tempPlasticStrainDeviator.computeDeviatoricVolumetricSum(answer,
                                                                 tempVolumetricPlasticStrain);
    }
    /**
     * Get the temp value of the plastic strain deviator from the material status.
     * @param answer Temp value of plastic strain deviator.
     */
    void giveTempPlasticStrainDeviator(StrainVector &answer) const { answer = tempPlasticStrainDeviator; }
    /**
     * Get the temp value of the volumetric strain deviator from the material status.
     * @return Temp value of volumetric plastic strain
     */
    double giveTempVolumetricPlasticStrain() const { return tempVolumetricPlasticStrain; }
    /**
     * Get the temp value of the hardening variable from the material status.
     * @return Temp value of hardening variable kappa.
     */
    double giveTempKappa() const { return tempKappa; }
    /**
     * Get the temp value of the state flag from the material status.
     * @return Temp value of the state flag (i.e. elastic, unloading, yielding, vertex case yielding)
     */
    int giveTempStateFlag() const { return temp_state_flag; }

    /**
     * Assign the temp value of deviatoric plastic strain.
     * @param v New temp value of deviatoric plastic strain.
     */
    void letTempPlasticStrainDeviatorBe(const StrainVector &v) { tempPlasticStrainDeviator = v; }
    /**
     * Assign the temp value of volumetric plastic strain.
     * @param v New temp value of volumetric plastic strain.
     */
    void letTempVolumetricPlasticStrainBe(double v) { tempVolumetricPlasticStrain = v; }
    /**
     * Assign the temp value of the hardening variable.
     * @param v New temp value of the hardening variable.
     */
    void letTempKappaBe(double v) { tempKappa = v; }
    /**
     * Assign the temp value of the state flag.
     * @param v New temp value of the state flag (i.e. elastic, unloading, yielding, vertex case yielding).
     */
    void letTempStateFlagBe(int v) { temp_state_flag = v; }
};

/**
 * This class implements a (local) nonassociated plasticity model based on the Drucker-Prager yield criterion with hardening and softening.
 * @author Simon Rolshoven
 */
class DruckerPragerPlasticitySM : public StructuralMaterial
{
public:
    /// Values of history variable state_flag.
    enum state_flag_values { DP_Elastic, DP_Unloading, DP_Yielding, DP_Vertex };

protected:
    /**
     * Controls the hardening function in the yield stress:
     * 1: linear hardening/softening with cutoff at zero stress.
     * 2: exponential hardening/softening to limitYieldStress.
     */
    int hardeningType;
    /// Parameter of the exponential laws.
    double kappaC;
    /// Hardening modulus normalized with the elastic modulus, parameter of the linear hardening/softening law.
    double hardeningModulus;
    /// Parameter of the exponential hardening law.
    double limitYieldStress;
    /// Parameter of all three laws, this is the initial value of the yield stress in pure shear.
    double initialYieldStress;
    /// Friction coefficient, parameter of the yield criterion.
    double alpha;
    /// Dilatancy coefficient, parameter of the flow rule.
    double alphaPsi;
    /// Scalar factor between rate of plastic multiplier and rate of hardening variable.
    double kFactor;

    /// Associated linear elastic material.
    IsotropicLinearElasticMaterial *LEMaterial;

    /// Hardening variable
    double kappa;
    double tempKappa;

    /// Volumetric stress, i.e. 1/3 of the trace of sigma.
    double volumetricStress;
    /// Volumetric part of elastic trial strain.
    double volumetricElasticTrialStrain;
    /// J2-invariant of elastic trial stress.
    double trialStressJTwo;
    /// Deviatoric part of stress.
    StressVector stressDeviator;
    /// Yield tolerance.
    double yieldTol;
    /// Maximum number of iterations for stress return.
    int newtonIter;

public:
    /// Constructor
    DruckerPragerPlasticitySM(int n, Domain *d);
    /// Destructor
    virtual ~DruckerPragerPlasticitySM();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int hasMaterialModeCapability(MaterialMode mMode);

    virtual const char *giveClassName() const { return "DruckerPragerPlasticitySM"; }
    virtual classType giveClassID() const { return DruckerPragerPlasticitySMClass; }

    virtual void giveRealStressVector(FloatArray &answer,
                              MatResponseForm form,
                              GaussPoint *gp,
                              const FloatArray &strainVector,
                              TimeStep *atTime);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form,
                                       MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep);

    /**
     * Perform a standard local stress return using the function computeYieldValue at the specified Gauss point.
     * This function computes the history variables, i.e. the temp variable of plastic strain, hardening variable, state flag, and the temp stress, and stores them in the temp-status.
     * @param gp Gauss point.
     * @param strain Strain vector of this Gauss point.
     */
    void performLocalStressReturn(GaussPoint *gp,
                                  const StrainVector &strain);
    /**
     * Check if the trial stress state falls within the vertex region.
     * @param eM Elasticity modulus.
     * @param gM Shear modulus.
     * @param kM Bulk modulus.
     * @return True for vertex case and false if regular stress return has to be used.
     */
    bool checkForVertexCase(double eM, double gM, double kM);
    /**
     * Perform stress return for regular case, i.e. if the trial stress state does not lie within the vertex region.
     * @param eM Elasticity modulus.
     * @param gM Shear modulus.
     * @param kM Bulk modulus.
     */
    void performRegularReturn(double eM, double gM, double kM);
    /**
     * Perform stress return for vertex case, i.e. if the trial stress state lies within the vertex region.
     * @param eM Elasticity modulus.
     * @param gM Shear modulus.
     * @param kM Bulk modulus.
     */
    void performVertexReturn(double eM, double gM, double kM);
    /**
     * Compute the yield value based on stress and hardening variable.
     * @param meanStress 1/3 of trace of sigma.
     * @param JTwo Second deviatoric invariant.
     * @param kappa Hardening variable.
     * @param eM Elasticity modulus.
     * @return Yield value.
     */
    double computeYieldValue(double meanStress,
                             double JTwo,
                             double kappa,
                             double eM) const;
    /**
     * Compute the current yield stress in pure shear of the Drucker-Prager model according to the used hardening law. The yield stress is tauY in f(sigma, kappa) = F(sigma) - tauY(kappa).
     * @param kappa Hardening variable.
     * @param eM Elasticity modulus.
     * @returns Yield stress in pure shear.
     */
    virtual double computeYieldStressInShear(double kappa, double eM) const;

    /**
     * Compute derivative of yield stress with respect to the hardening variable kappa.
     * @param kappa Hardening variable.
     * @param eM Elasticity modulus.
     * @return Derivative of yield stress with respect to kappa.
     */
    virtual double computeYieldStressPrime(double kappa, double eM) const;

    /**
     * Compute and give back algorithmic stiffness matrix for the regular case (no vertex).
     * @param answer Consistent stiffness matrix.
     * @param form Material response form.
     * @param mode Material reponse mode.
     * @param gp Gauss point.
     * @param tStep Time step.
     */
    void giveRegAlgorithmicStiffMatrix(FloatMatrix &answer,
                                       MatResponseForm form,
                                       MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *tStep);
    /**
     * Compute consistent stiffness matrix for the vertex case.
     * @param answer Consistent stiffness matrix.
     * @param form Material response form.
     * @param mode Material reponse mode.
     * @param gp Gauss point.
     * @param tStep Time step.
     */
    void giveVertexAlgorithmicStiffMatrix(FloatMatrix &answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer,
                    GaussPoint *gp,
                    InternalStateType type,
                    TimeStep *atTime);

    virtual int giveIPValueSize(InternalStateType type,
                        GaussPoint *gp);

    virtual int giveIntVarCompFullIndx(IntArray &answer,
                               InternalStateType type,
                               MaterialMode mmode);

    virtual InternalStateValueType giveIPValueType(InternalStateType type);

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *atTime)
    {
        LEMaterial->giveThermalDilatationVector(answer, gp, atTime);
    }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

#ifdef __PARALLEL_MODE
    virtual double predictRelativeComputationalCost(GaussPoint *gp);
    virtual double predictRelativeRedistributionCost(GaussPoint *gp) { return 1.0; }
#endif
};
} // end namespace oofem
#endif // druckerpragerplasticitysm_h
