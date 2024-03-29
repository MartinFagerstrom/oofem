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


#ifndef trabbonenl3d_h
#define trabbonenl3d_h

#include "trabbone3d.h"
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

/**
 * Trabecular bone nonlocal material status.
 */
class TrabBoneNL3DStatus : public TrabBone3DStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    /// Equivalent strain for averaging.
    double localCumPlastStrainForAverage;

public:
    TrabBoneNL3DStatus(int n, Domain *d, GaussPoint *g);
    virtual ~TrabBoneNL3DStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    double giveLocalCumPlastStrainForAverage() { return localCumPlastStrainForAverage; }
    const FloatArray *giveLTangentContrib();
    void setLocalCumPlastStrainForAverage(double ls) { localCumPlastStrainForAverage = ls; }

    // definition
    virtual const char *giveClassName() const { return "TrabBoneNL3DStatus"; }
    virtual classType giveClassID() const { return TrabBone3DStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
    virtual Interface *giveInterface(InterfaceType it);
};


/**
 * Trabecular bone nonlocal material model.
 */
class TrabBoneNL3D : public TrabBone3D,
    public StructuralNonlocalMaterialExtensionInterface,
    public NonlocalMaterialStiffnessInterface
{
protected:
    double R;
    double mParam;

public:
    TrabBoneNL3D(int n, Domain *d);
    virtual ~TrabBoneNL3D();

    virtual const char *giveClassName() const { return "TrabBoneNL3D"; }
    virtual classType giveClassID() const { return TrabBoneNL3DClass; }
    virtual const char *giveInputRecordName() const { return "trabbonenl3d"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    virtual Interface *giveInterface(InterfaceType it);


    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);

    /**
     * Computes the local cumulated plastic strain from given strain vector (full form).
     * @param[out] kappa Return parameter, containing the corresponding cumulated plastic strain.
     * @param strain Total strain vector in full form.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    void computeLocalCumPlastStrain(double &kappa, const StrainVector &strain, GaussPoint *gp, TimeStep *tStep)
    {
        TrabBone3D :: computeCumPlastStrain(kappa, gp, tStep);
    }

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,  MatResponseForm form, MatResponseMode mode, GaussPoint *gp,  TimeStep *tStep);

#ifdef __OOFEG
    // Plots the sparse structure of stiffness contribution.
    //virtual void NonlocalMaterialStiffnessInterface_showSparseMtrxStructure(GaussPoint *gp, oofegGraphicContext &gc, TimeStep *tStep);
#endif

    virtual void NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s,
                                                                      GaussPoint *gp, TimeStep *atTime);

    virtual dynaList< localIntegrationRecord > *NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp);

    /**
     * Computes the "local" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp Source integration point.
     * @param loc Local code numbers.
     * @param lcontrib "Local" contribution.
     * @param s  Numbering scheme for unknowns.
     * @param tStep Time step.
     * @return Nonzero if local point contributes (loading) or zero if not (unloading in elastic range, elastic).
     */
    int giveLocalNonlocalStiffnessContribution(GaussPoint *gp, IntArray &loc, const UnknownNumberingScheme &s,
                                               FloatArray &lcontrib, TimeStep *tStep);

    /**
     * Computes the "remote" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp Remote integration point.
     * @param rloc Remote element code numbers.
     * @param rcontrib "Remote" contribution.
     * @param s  Numbering scheme for unknowns.
     * @param tStep Time step.
     */
    void giveRemoteNonlocalStiffnessContribution(GaussPoint *gp, IntArray &rloc, const UnknownNumberingScheme &s,
                                                 FloatArray &rcontrib, TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer,  MatResponseForm form, GaussPoint *gp, const FloatArray &strainVector, TimeStep *tStep);

    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep);

    virtual double computeWeightFunction(const FloatArray &src, const FloatArray &coord);

    virtual int hasBoundedSupport() { return 1; }

    /// Determines the width (radius) of limited support of weighting function.
    virtual void giveSupportRadius(double &radius) { radius = this->R; }

#ifdef __PARALLEL_MODE
    virtual int packUnknowns(CommunicationBuffer &buff, TimeStep *tStep, GaussPoint *ip);
    virtual int unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *tStep, GaussPoint *ip);
    virtual int estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip);
#endif

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new TrabBoneNL3DStatus(1, TrabBone3D :: domain, gp); }
};
} // end namespace oofem
#endif
