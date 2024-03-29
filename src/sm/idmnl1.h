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

#ifndef idmnl1_h
#define idmnl1_h

#include "idm1.h"
#include "structuralnonlocalmaterialext.h"
#include "nonlocmatstiffinterface.h"

namespace oofem {
class GaussPoint;

/**
 * This class implements associated Material Status to IDNLMaterial (Nonlocal isotropic damage).
 * Stores local equivalent strain for averaging.
 */
class IDNLMaterialStatus : public IsotropicDamageMaterial1Status, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    /// Equivalent strain for averaging.
    double localEquivalentStrainForAverage;

    /* // Variables used to track loading/reloading
     * public:
     * enum LastStateType {LST_elastic, LST_loading, LST_unloading};
     * LastStateType lst;
     */

public:
    /// Constructor.
    IDNLMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor.
    virtual ~IDNLMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    /// Returns the local  equivalent strain to be averaged.
    double giveLocalEquivalentStrainForAverage() { return localEquivalentStrainForAverage; }
    /// Sets the localEquivalentStrainForAverage to given value.
    void setLocalEquivalentStrainForAverage(double ls) { localEquivalentStrainForAverage = ls; }

    // definition
    virtual const char *giveClassName() const { return "IDNLMaterialStatus"; }
    virtual classType giveClassID() const { return IsotropicDamageMaterialStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /**
     * Interface requesting service.
     * In the case of nonlocal constitutive models,
     * the use of multiple inheritance is assumed. Typically, the class representing nonlocal
     * constitutive model status is derived both from class representing local status and from class
     * NonlocalMaterialStatusExtensionInterface or from one of its derived classes
     * (which declare services and variables corresponding to specific analysis type).
     * @return In both cases, this function returns pointer to this object, obtained by
     * returning address of component or using pointer conversion from receiver to base class
     * NonlocalMaterialStatusExtensionInterface.
     */
    virtual Interface *giveInterface(InterfaceType);
};


/**
 * This class implements a Nonlocal Isotropic Damage Model for Concrete in Tension
 * Model based on nonlocal averaging of equivalent strain.
 */
class IDNLMaterial : public IsotropicDamageMaterial1, public StructuralNonlocalMaterialExtensionInterface,
    public NonlocalMaterialStiffnessInterface
{
protected:
    /// Final value of interaction radius, for a model with evolving characteristic length.
    double Rf;
    /// Parameter used as an exponent by models with evolving characteristic length.
    double exponent;
    /// Parameter specifying how the weight function should be adjusted due to damage.
    int averType;

public:
    /// Constructor
    IDNLMaterial(int n, Domain *d);
    /// Destructor
    virtual ~IDNLMaterial();

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "IDNLMaterial"; }
    virtual classType giveClassID() const { return IDNLMaterialClass;}
    virtual const char *giveInputRecordName() const { return "idmnl1"; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    virtual Interface *giveInterface(InterfaceType it);

    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    void computeLocalEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
    { IsotropicDamageMaterial1 :: computeEquivalentStrain(kappa, strain, gp, tStep); }

    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep);

    void modifyNonlocalWeightFunctionAround(GaussPoint *gp);
    double computeDistanceModifier(double damage);

#ifdef __OOFEG
    /// Plots the sparse structure of stiffness contribution.
    virtual void NonlocalMaterialStiffnessInterface_showSparseMtrxStructure(GaussPoint *gp, oofegGraphicContext &gc, TimeStep *atTime);
#endif

    virtual void computeDamageParam(double &omega, double kappa, const FloatArray &strain, GaussPoint *gp);

    /**@name Services required by NonlocalMaterialStiffnessInterface and related ones to support Nonlocal Stiffness*/
    //@{
    virtual void NonlocalMaterialStiffnessInterface_addIPContribution(SparseMtrx &dest, const UnknownNumberingScheme &s,
                                                                      GaussPoint *gp, TimeStep *tStep);
    /**
     * Returns integration list of receiver. Contains localIntegrationRecord structures, containing
     * references to integration points and their weights that influence to nonlocal average in
     * receiver's associated integration point.
     */
    virtual dynaList< localIntegrationRecord > *NonlocalMaterialStiffnessInterface_giveIntegrationDomainList(GaussPoint *gp);
    /**
     * Computes the "local" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp Source integration point.
     * @param loc Local code numbers.
     * @param s Determines the equation numbering scheme.
     * @param lcontrib "Local" contribution.
     * @param tStep Time step.
     * @return Nonzero if local point contributes (loading) or zero if not (unloading in elastic range, elastic)
     */
    int giveLocalNonlocalStiffnessContribution(GaussPoint *gp, IntArray &loc, const UnknownNumberingScheme &s,
                                               FloatArray &lcontrib, TimeStep *tStep);
    /**
     * Computes the "remote" part of nonlocal stiffness contribution assembled for given integration point.
     * @param gp Remote integration point.
     * @param rloc Remote element code numbers.
     * @param s Determines the equation numbering scheme.
     * @param rcontrib "Remote" contribution.
     * @param tStep Time step.
     */
    void giveRemoteNonlocalStiffnessContribution(GaussPoint *gp, IntArray &rloc, const UnknownNumberingScheme &s,
                                                 FloatArray &rcontrib, TimeStep *tStep);
    /**
     * Computes elastic stiffness for normal stress components.
     * @param answer Result of size (3,3).
     * @param mode Determines the MatResponseMode.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    void giveNormalElasticStiffnessMatrix(FloatMatrix &answer,
                                          MatResponseMode mode,
                                          GaussPoint *gp, TimeStep *tStep);
    //@}

#ifdef __PARALLEL_MODE
    int packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip);
    int unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip);
    int estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip);
    virtual double predictRelativeComputationalCost(GaussPoint *gp);
    virtual double predictRelativeRedistributionCost(GaussPoint *gp) { return 1.0; }
#endif

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new IDNLMaterialStatus(1, IsotropicDamageMaterial1 :: domain, gp); }

protected:
    virtual void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp) { }
};
} // end namespace oofem
#endif // idmnl1_h
