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

#ifndef misesmatnl_h

#include "misesmat.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"
#include "cltypes.h"

#include "dynalist.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {

/**
 * Mises Nonlocal material status.
 * @author Milan
 */
class MisesMatNlStatus : public MisesMatStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    // STATE VARIABLE DECLARATION
    // Equivalent strain for avaraging
    double localCumPlasticStrainForAverage;

public:
    MisesMatNlStatus(int n, Domain *d, GaussPoint *g);
    virtual ~MisesMatNlStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // STATE VARIABLE
    // declare state variable access and modification methods
    double giveLocalCumPlasticStrainForAverage() { return localCumPlasticStrainForAverage; }
    const FloatArray *giveLTangentContrib();
    void setLocalCumPlasticStrainForAverage(double ls) { localCumPlasticStrainForAverage = ls; }

    // DEFINITION
    virtual const char *giveClassName() const { return "MisesMatNlStatus"; }
    virtual classType giveClassID() const { return MisesMatClass; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual Interface *giveInterface(InterfaceType);
};


/**
 * Mises nonlocal material.
 * @author Milan
 */
class MisesMatNl : public MisesMat, public StructuralNonlocalMaterialExtensionInterface,
    public NonlocalMaterialStiffnessInterface
{
protected:
    double Rf;
    double exponent;
    int averType;

public:
    MisesMatNl(int n, Domain *d);
    virtual ~MisesMatNl();

    virtual const char *giveClassName() const { return "MisesMatNl"; }
    virtual classType giveClassID() const { return MisesMatClass; }
    virtual const char *giveInputRecordName() const { return "MisesMatNl"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    virtual Interface *giveInterface(InterfaceType);

    /**
     * Computes the nonlocal cumulated plastic strain from its local form.
     * @param[out] kappa Return param, containing the corresponding cumulated plastic strain.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void computeCumPlasticStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);
    double computeDamage(GaussPoint *gp, TimeStep *atTime);
    void modifyNonlocalWeightFunctionAround(GaussPoint *gp);
    double computeDistanceModifier(double damage);
    void computeLocalCumPlasticStrain(double &kappa, GaussPoint *gp, TimeStep *atTime)
    {
        MisesMat :: computeCumPlastStrain(kappa, gp, atTime);
    }

    virtual void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep);
    //virtual void givePlaneStrainStiffMtrx(FloatMatrix& answer,MatResponseForm, MatResponseMode, GaussPoint *gp,TimeStep *atTime);
    //virtual void give3dMaterialStiffnessMatrix(FloatMatrix& answer, MatResponseForm, MatResponseMode, GaussPoint *gp, TimeStep *atTime);

#ifdef __OOFEG
#endif

    virtual void NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s,
                                                                      GaussPoint *gp, TimeStep *atTime);

    virtual dynaList< localIntegrationRecord > *NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp);

    /**
     * Computes the "local" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp Source integration point
     * @param loc Local code numbers
     * @param lcontrib "Local" contribution
     * @return Nonzero if local point contributes (loading) or zero if not (unloading in elastic range, elastic)
     */
    int giveLocalNonlocalStiffnessContribution(GaussPoint *gp, IntArray &loc, const UnknownNumberingScheme &s,
                                               FloatArray &lcontrib, TimeStep *atTime);

    /**
     * Computes the "remote" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp Remote integration point.
     * @param loc Remote element code numbers.
     * @param rcontrib "Remote" contribution.
     */
    void giveRemoteNonlocalStiffnessContribution(GaussPoint *gp, IntArray &rloc, const UnknownNumberingScheme &s,
                                                 FloatArray &rcontrib, TimeStep *atTime);

    virtual void giveRealStressVector(FloatArray &answer,  MatResponseForm form, GaussPoint *gp, const FloatArray &strainVector, TimeStep *atTime);

    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *atTime);

    virtual int hasBoundedSupport() { return 1; }

#ifdef __PARALLEL_MODE
    int packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip);
    int unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip);
    int estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip);
#endif

protected:
    // Creates the corresponding material status
    MaterialStatus *CreateStatus(GaussPoint *gp) const { return new MisesMatNlStatus(1, MisesMat :: domain, gp); }
};
} // end namespace oofem
#define misesmatnl_h
#endif
