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

#ifndef perfectlyplasticmaterial_h
#define perfectlyplasticmaterial_h

#include "structuralmaterial.h"
#include "linearelasticmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "structuralms.h"

namespace oofem {

/**
 * This class implements associated Material Status to PerfectlyPlasticMaterial.
 *
 * This class contains only yield_flag to inform, whether
 * gp is yielding or is in elastic state.
 * temp_yield_flag is used to indicate yield_flag during a process when we try
 * to find a equilibrium state. After This State is reached we update the yield_Flag
 * accordingly.
 * When lastly reached state is equilibrium state we call updateYourself(), so
 * yield_flag containing last reached status of yielding is
 * copied into yield_flag indicating yield status at previously
 * reached equilibrium state.
 *
 * Tasks:
 * - Storing the yield_flag
 * - Storing the plastic strain + increment.
 */
class PerfectlyPlasticMaterialStatus : public StructuralMaterialStatus
{
protected:
    FloatArray plasticStrainVector; // reduced form to save memory
    FloatArray plasticStrainIncrementVector;
    int yield_flag;
    int temp_yield_flag;

public:
    PerfectlyPlasticMaterialStatus(int n, Domain *d, GaussPoint *g);
    virtual ~PerfectlyPlasticMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    int setTempYieldFlag(int i) { return temp_yield_flag = i; }
    int giveTempYieldFlag() { return temp_yield_flag; }
    int giveYieldFlag() { return yield_flag; }

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    void givePlasticStrainVector(FloatArray &answer) const { answer =  plasticStrainVector; }
    void givePlasticStrainIncrementVector(FloatArray &answer) const { answer = plasticStrainIncrementVector; }
    void letPlasticStrainVectorBe(FloatArray &v) { plasticStrainVector = v; }
    void letPlasticStrainIncrementVectorBe(FloatArray &v) { plasticStrainIncrementVector = v; }

    // definition
    virtual const char *giveClassName() const { return "PerfectlyPlasticMaterialStatus"; }
    virtual classType giveClassID() const { return PerfectlyPlasticMaterialStatusClass; }
};

/**
 * This class implements a perfectly plastic material in a finite element problem. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 *
 * This class is mainly abstract, no new methods or data added. This class is
 * mainly intended to be a base class for perfectly plastic materials.
 * (for perfectly plastic material a linear elastic behaviour is assumed
 * before yielding surface is being reached, and also unloading is assumed to
 * be linear elastic).
 * This class takes care about possible failure and fracture and hardening,
 * even if it is not required
 * for perfectly-plastic material, but this class is assumed to be a parent
 * class for more sophisticated materials where these phenomena should appear.
 * Tasks:
 * - Returning standard material stiffness and flexibility matrices for 3d-case.
 *   according to current state determined by using data stored
 *   in Gausspoint.
 * - Returning a material property (method 'give'). Only for non-standard elements.
 * - Returning real stress state vector(tensor) at Gauss point for 3d-case.
 */
class PerfectlyPlasticMaterial : public StructuralMaterial
{
protected:
    int yieldCriteria;
    int loadingCriteria;
    LinearElasticMaterial *linearElasticMaterial;

public:

    PerfectlyPlasticMaterial(int n, Domain *d) : StructuralMaterial(n, d)
    {
        yieldCriteria = 0;
        loadingCriteria = 0;
        linearElasticMaterial = NULL;
    }

    virtual ~PerfectlyPlasticMaterial() { delete linearElasticMaterial; }

    virtual void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void updateYourself(GaussPoint *gp, TimeStep *tStep);

    // can be moved to parent classes
    virtual void updateIfFailure(GaussPoint *gp,
                                 FloatArray *stressVector3d,
                                 FloatArray *PlasticStrainVector3d) { }

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour() { return 1; }
    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual const char *giveClassName() const { return "PerfectlyPlasticMaterial"; }
    virtual classType giveClassID() const { return PerfectlyPlasticMaterialClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual double give(int aProperty, GaussPoint *gp);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    // two functions used to initialize and updating temporary variables in
    // gps status. These variables are used to control process, when
    // we try to find equilibrium state

    virtual void giveMaterialStiffnessMatrix(FloatMatrix &answer,
                                             MatResponseMode mode, GaussPoint *gp,
                                             TimeStep *tStep);
    virtual void giveEffectiveMaterialStiffnessMatrix(FloatMatrix &answer,
                                                      MatResponseMode mode, GaussPoint *gp,
                                                      TimeStep *tStep);

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
    virtual void give3dShellLayerStiffMtrx(FloatMatrix &answer,
                                           MatResponseForm form, MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep);

    void computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                     const FloatArray &strainIncrement, TimeStep *tStep);
    void computePlasticStiffnessAt(FloatMatrix &answer,
                                   GaussPoint *gp,
                                   FloatArray *currentStressVector,
                                   FloatArray *currentPlasticStrainVector,
                                   FloatArray *strainIncrement3d,
                                   TimeStep *tStep,
                                   double &lambda);

    FloatArray *GiveStressCorrectionBackToYieldSurface(GaussPoint *gp,
                                                       FloatArray *stressVector3d,
                                                       FloatArray *plasticVector3d);
    //
    // yield(YC-like functions) and loading(LC-like functions) criteria specific section
    //

    virtual double computeYCValueAt(GaussPoint *gp, FloatArray *, FloatArray *) { return -1.0; }
    virtual FloatArray *GiveYCStressGradient(GaussPoint *gp, FloatArray *, FloatArray *) { return NULL; }
    virtual FloatArray *GiveLCStressGradient(GaussPoint *gp, FloatArray *, FloatArray *) { return NULL; }
    virtual FloatArray *GiveYCPlasticStrainGradient(GaussPoint *gp, FloatArray *, FloatArray *) { return NULL; }
    virtual FloatArray *GiveLCPlasticStrainGradient(GaussPoint *gp, FloatArray *, FloatArray *) { return NULL; }
    //virtual FloatArray* GiveHardeningGradient (GaussPoint *gp,FloatArray *,FloatArray *) {return NULL;}
    virtual void updateTempYC(GaussPoint *gp, FloatArray *, FloatArray *) { }
    virtual void updateTempLC(GaussPoint *gp, FloatArray *, FloatArray *) { }
    //virtual int updateYieldStatus(GaussPoint* gp, FloatArray* strainIncrementIn3d);
};
} // end namespace oofem
#endif // perfectlyplasticmaterial_h
