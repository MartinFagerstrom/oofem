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

#ifndef kelvinchsol_h
#define kelvinchsol_h

#include "rheoChM.h"

namespace oofem {
/**
 * This class implements associated Material Status to KelvinChainSolidMaterial,
 * which corresponds to a solidifying Kelvin chain model (framework for creep with aging).
 */
class KelvinChainSolidMaterialStatus : public RheoChainMaterialStatus
{
public:
    KelvinChainSolidMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits);
    virtual ~KelvinChainSolidMaterialStatus() {}
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // definition
    virtual const char *giveClassName() const { return "KelvinChainSolidMaterialStatus"; }
    virtual classType giveClassID() const { return KelvinChainSolidMaterialStatusClass; }
};


/**
 * This class implements a solidifying Kelvin chain model
 * describing a viscoelastic material.
 */
class KelvinChainSolidMaterial : public RheoChainMaterial
{
public:
    KelvinChainSolidMaterial(int n, Domain *d);
    virtual ~KelvinChainSolidMaterial() {}

    virtual void updateYourself(GaussPoint *gp, TimeStep *tStep);

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour() { return 0; }
    virtual const char *giveClassName() const { return "KelvinChainSolidMaterial"; }
    virtual classType giveClassID() const { return KelvinChainSolidMaterialClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void  giveShrinkageStrainVector(FloatArray &answer,
                                            MatResponseForm form,
                                            GaussPoint *gp,
                                            TimeStep *tStep,
                                            ValueModeType mode)
    { answer.resize(0); }

    virtual void  giveEigenStrainVector(FloatArray &answer, MatResponseForm form,
                                        GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;


protected:
    virtual int hasIncrementalShrinkageFormulation() { return 0; }

    /// Evaluation of the creep compliance function - function useless here
    virtual double computeCreepFunction(GaussPoint *gp, double ofAge, double atTime);

    virtual double giveEModulus(GaussPoint *gp, TimeStep *tStep);

    /// Evaluation of the relative volume of the solidified material
    virtual double computeSolidifiedVolume(GaussPoint *gp, TimeStep *tStep) = 0;

    /// factors for exponential algorithm
    virtual double computeBetaMu(GaussPoint *gp, TimeStep *tStep, int Mu);
    virtual double computeLambdaMu(GaussPoint *gp, TimeStep *tStep, int Mu);

    LinearElasticMaterial *giveLinearElasticMaterial();
};
} // end namespace oofem
#endif // kelvinchsol_h
