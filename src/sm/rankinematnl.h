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


#ifndef rankinematnl_h

#include "rankinemat.h"
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
 * Rankine nonlocal material status.
 */
class RankineMatNlStatus : public RankineMatStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    /// Equivalent strain for averaging.
    double localCumPlasticStrainForAverage;

    /// For printing only
    double kappa_nl;
    double kappa_hat;

public:
    RankineMatNlStatus(int n, Domain *d, GaussPoint *g);
    virtual ~RankineMatNlStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    double giveLocalCumPlasticStrainForAverage() { return localCumPlasticStrainForAverage; }
    const FloatArray *giveLTangentContrib();
    void setLocalCumPlasticStrainForAverage(double ls) { localCumPlasticStrainForAverage = ls; }

    virtual const char *giveClassName() const { return "RankineMatNlStatus"; }
    virtual classType giveClassID() const { return RankineMatClass; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);

    void setKappa_nl(double kap) { kappa_nl = kap; }
    void setKappa_hat(double kap) { kappa_hat = kap; }
    double giveKappa_nl() { return kappa_nl; }
    double giveKappa_hat() { return kappa_hat; }

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual Interface *giveInterface(InterfaceType);
};


/**
 * Rankine nonlocal material
 */
class RankineMatNl : public RankineMat, public StructuralNonlocalMaterialExtensionInterface,
    public NonlocalMaterialStiffnessInterface
{
public:
    RankineMatNl(int n, Domain *d);
    virtual ~RankineMatNl() {; }

    virtual const char *giveClassName() const { return "RankineMatNl"; }
    virtual classType giveClassID() const { return RankineMatNlClass; }
    virtual const char *giveInputRecordName() const { return "RankineMatNl"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    virtual Interface *giveInterface(InterfaceType);

    /**
     * Computes the nonlocal cumulated plastic strain from its local form.
     * @param kappa return param, containing the corresponding cumulated plastic strain.
     * @param gp integration point.
     * @param atTime time step.
     */
    virtual void computeCumPlasticStrain(double &kappa, GaussPoint *gp, TimeStep *atTime);
    double computeDamage(GaussPoint *gp, TimeStep *atTime);
    void modifyNonlocalWeightFunctionAround(GaussPoint *gp);
    double computeDistanceModifier(double damage);
    void computeLocalCumPlasticStrain(double &kappa, GaussPoint *gp, TimeStep *atTime)
    {
        RankineMat :: computeCumPlastStrain(kappa, gp, atTime);
    }


    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep);
    //virtual void givePlaneStrainStiffMtrx(FloatMatrix& answer,MatResponseForm, MatResponseMode,GaussPoint * gp,TimeStep * atTime);
    // virtual void give3dMaterialStiffnessMatrix (FloatMatrix& answer,  MatResponseForm, MatResponseMode,GaussPoint* gp,  TimeStep* atTime);

#ifdef __OOFEG
    // Plots the sparse structure of stiffness contribution.
    //  virtual void NonlocalMaterialStiffnessInterface_showSparseMtrxStructure(GaussPoint* gp, oofegGraphicContext& gc, TimeStep* atTime);
#endif

    virtual void NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s,
                                                                      GaussPoint *gp, TimeStep *atTime);

    virtual dynaList< localIntegrationRecord > *NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp);

    /**
     * Computes the "local" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp Source integration point.
     * @param loc Local code numbers.
     * @param lcontrib "Local" contribution.
     * @return Nonzero if local point contributes (loading) or zero if not (unloading in elastic range, elastic).
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

    // Computes elastic stiffness for normal stress components
    // @param answer result of size (3,3)
    // @param mode determines the MatResponseMode
    // @param gp integration point
    // @param atTime time step
    //  void giveNormalElasticStiffnessMatrix (FloatMatrix& answer, MatResponseMode rMode, GaussPoint*gp, TimeStep* atTime) ;

    virtual void giveRealStressVector(FloatArray &answer,  MatResponseForm form, GaussPoint *gp, const FloatArray &strainVector, TimeStep *atTime);

    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *atTime);

    virtual int hasBoundedSupport() { return 1; }

    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *gp);

#ifdef __PARALLEL_MODE
    virtual int packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip);
    virtual int unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip);
    virtual int estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip);
#endif

protected:
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new RankineMatNlStatus(1, RankineMat :: domain, gp); }
};
} // end namespace oofem
#define rankinematnl_h
#endif
