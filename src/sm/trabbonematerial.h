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

#ifndef trabbonematerial_h

#include "structuralmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "matconst.h"
#include "matstatus.h"
#include "strainvector.h"

#include "linearelasticmaterial.h"
#include "dictionr.h"

#include "structuralms.h"
#include "cltypes.h"

namespace oofem {

/**
 * This class implements associated Material Status to TrabBoneMaterial.
 */
class TrabBoneMaterialStatus : public StructuralMaterialStatus
{
protected:
    double tempAlpha, alpha;
    double tempDam, dam;
    double smtrx, slope;
    double sigC, matConstC;
    FloatArray tempEpsp, epsp, tempDepsp;

public:
    TrabBoneMaterialStatus(int n, Domain *d, GaussPoint *g);
    virtual ~TrabBoneMaterialStatus();

    void printOutputAt(FILE *file, TimeStep *tStep);

    double giveAlpha();
    double giveTempAlpha();
    double giveDam();
    double giveTempDam();
    double giveSmtrx();
    double giveSlope();
    double giveSigC();
    double giveMatConstC();
    FloatArray givePlasStrainVector();
    FloatArray giveTempPlasStrainVector();
    FloatArray giveTempIncPlasStrainVector();

    void setTempAlpha(double al) { tempAlpha = al; }
    void setTempDam(double da) { tempDam = da; }
    void setSmtrx(double smt) { smtrx = smt; }
    void setSlope(double slp) { slope = slp; }
    void setSigC(double sc) { sigC = sc; }
    void setMatConstC(double mcc) { matConstC = mcc; }
    void setTempEpsp(double epsip) { tempEpsp.at(1) = epsip; }
    void setTempDepsp(double depsip) { tempDepsp.at(1) = depsip; }


    // definition
    virtual const char *giveClassName() const { return "TrabBoneMaterialStatus"; }
    virtual classType giveClassID() const { return TrabBoneMaterialStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};



/**
 * Trabecular bone material model.
 */
class TrabBoneMaterial : public StructuralMaterial
{
protected:
    double E0, Eil, Eie, kie, Ek, Cc, Cc2, EpsC, SigYp, SigYn, adam;

public:
    TrabBoneMaterial(int n, Domain *d);

    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain);

    void computeDensification(GaussPoint *gp, const FloatArray &totalStrain);

    double computeDamageParam(double alpha, GaussPoint *gp);

    double computeDamage(GaussPoint *gp, TimeStep *tStep);

    virtual void computeCumPlastStrain(double &alpha, GaussPoint *gp, TimeStep *atTime);

    virtual void give1dStressStiffMtrx(FloatMatrix &answer,
                                       MatResponseForm form, MatResponseMode mode, GaussPoint *gp,
                                       TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual int hasMaterialModeCapability(MaterialMode);

    virtual const char *giveClassName() const { return "TrabBoneMaterial"; }
    virtual classType giveClassID() const { return TrabBoneMaterialClass; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
};
} // end namespace oofem
#define trabbonematerial_h
#endif
