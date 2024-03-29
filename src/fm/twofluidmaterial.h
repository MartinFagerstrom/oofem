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

#ifndef twofluidmaterial_h
#define twofluidmaterial_h

#include "fluiddynamicmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "matconst.h"
#include "structuralelement.h"
#include "matstatus.h"

namespace oofem {
class GaussPoint;

/**
 * Material coupling the behavior of two particular materials based on
 * rule of mixture. The weighting factor is VOF fraction.
 */
class TwoFluidMaterial : public FluidDynamicMaterial
{
protected:
    int masterMat;
    int slaveMaterial [ 2 ];

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    TwoFluidMaterial(int n, Domain *d) : FluidDynamicMaterial(n, d) { masterMat = 0; }
    /// Destructor.
    virtual ~TwoFluidMaterial() { }

    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime) { }

    virtual double  giveCharacteristicValue(MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *atTime);

    virtual void computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep);
    virtual void giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual double give(int aProperty, GaussPoint *gp);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual const char *giveClassName() const { return "TwoFluidMaterial"; }
    virtual classType giveClassID() const { return TwoFluidMaterialClass; }
    virtual int checkConsistency();
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    FluidDynamicMaterial *giveMaterial(int i) const { return static_cast< FluidDynamicMaterial * >( domain->giveMaterial(slaveMaterial [ i ]) ); }
    double giveTempVOF(GaussPoint *gp);
};
} // end namespace oofem
#endif // twofluidmaterial_h
