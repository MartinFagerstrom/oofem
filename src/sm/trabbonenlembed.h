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

#ifndef trabbonenlembed_h
#define trabbonenlembed_h

#include "trabboneembed.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "cltypes.h"

#include "sparsemtrx.h"
#include "dynalist.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
class GaussPoint;


/**
 * Trabecular bone nonlocal material status.
 */
class TrabBoneNLEmbedStatus : public TrabBoneEmbedStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    /// Equivalent strain for averaging.
    double localCumPlastStrainForAverage;

public:
    TrabBoneNLEmbedStatus(int n, Domain *d, GaussPoint *g);
    virtual ~TrabBoneNLEmbedStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    /// Gives the local cumulative plastic strain.
    double giveLocalCumPlastStrainForAverage() { return localCumPlastStrainForAverage; }
    /// Sets the local cumulative plastic strain.
    void setLocalCumPlastStrainForAverage(double ls) { localCumPlastStrainForAverage = ls; }

    // definition
    virtual const char *giveClassName() const { return "TrabBoneNLEmbedStatus"; }
    virtual classType giveClassID() const { return TrabBoneEmbedStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual Interface *giveInterface(InterfaceType it);
};

/**
 * Trabecular bone nonlocal material.
 */
class TrabBoneNLEmbed : public TrabBoneEmbed, public StructuralNonlocalMaterialExtensionInterface
{
protected:
    double R;
    double mParam;

public:
    TrabBoneNLEmbed(int n, Domain *d);
    virtual ~TrabBoneNLEmbed();

    virtual const char *giveClassName() const { return "TrabBoneNLEmbed"; }
    virtual classType giveClassID() const { return TrabBoneNLEmbedClass; }
    virtual const char *giveInputRecordName() const { return "trabbonenlembed"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    virtual Interface *giveInterface(InterfaceType);

    virtual void computeCumPlastStrain(double &alpha, GaussPoint *gp, TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer,  MatResponseForm form, GaussPoint *gp, const FloatArray &strainVector, TimeStep *tStep);

    void computeLocalCumPlastStrain(double &alpha, const StrainVector &strain, GaussPoint *gp, TimeStep *tStep)
    {
        TrabBoneEmbed :: computeCumPlastStrain(alpha, gp, tStep);
    }

    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep);
    virtual double computeWeightFunction(const FloatArray &src, const FloatArray &coord);

    virtual int hasBoundedSupport() { return 1; }

    /// Determines the width (radius) of limited support of weighting function.
    virtual void giveSupportRadius(double &radius) { radius = this->R; }

protected:
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new TrabBoneNLEmbedStatus(1, TrabBoneEmbed :: domain, gp); }
};
} // end namespace oofem
#endif
