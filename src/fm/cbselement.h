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

#ifndef cbselement_h
#define cbselement_h

#include "fmelement.h"

namespace oofem {
class TimeStep;
class GaussPoint;
class FloatMatrix;
class FloatArray;

/**
 * This abstract class represent a general base element class for
 * fluid dynamic problems solved using CBS algorithm.
 */
class CBSElement : public FMElement
{
public:
    CBSElement(int n, Domain *aDomain);
    virtual ~CBSElement();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep);
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual double giveCharacteristicValue(CharType type, TimeStep *tStep);

    /** Calculates consistent mass matrix. */
    virtual void computeConsistentMassMtrx(FloatMatrix &answer, TimeStep *tStep) = 0;
    /** Calculates diagonal mass matrix. */
    virtual void computeDiagonalMassMtrx(FloatArray &answer, TimeStep *tStep) = 0;
    /** Calculates rhs due to prescribed (*) velocities for (*) velocities. */
    virtual void computePrescribedTermsI(FloatArray &answer, ValueModeType mode, TimeStep *tStep);
    /* Calculates rhs due to prescribed density/pressure. */
    //virtual void computePrescribedTermsII(FloatArray& answer, ValueModeType, TimeStep*);
    /** Calculates convection component for (*) velocities. */
    virtual void computeConvectionTermsI(FloatArray &answer, TimeStep *tStep) = 0;
    /** Calculates contribution from diffusion terms for (*) velocities. */
    virtual void computeDiffusionTermsI(FloatArray &answer, TimeStep *tStep) = 0;
    /// Computes velocity terms on RHS for density equation.
    virtual void computeDensityRhsVelocityTerms(FloatArray &answer, TimeStep *tStep) = 0;
    /// Computes pressure terms on RHS for density equation.
    virtual void computeDensityRhsPressureTerms(FloatArray &answer, TimeStep *tStep) = 0;
    /// Computes prescribed pressure due to applied tractions.
    virtual void computePrescribedTractionPressure(FloatArray &answer, TimeStep *tStep) = 0;
    /// Computes number of edges/sides with prescribed traction contributing to node with prescribed pressure.
    virtual void computeNumberOfNodalPrescribedTractionPressureContributions(FloatArray &answer, TimeStep *tStep) = 0;
    /// Calculates the pressure LHS.
    virtual void computePressureLhs(FloatMatrix &answer, TimeStep *tStep) = 0;
    /// Calculates the RHS of velocity correction step.
    virtual void computeCorrectionRhs(FloatArray &answer, TimeStep *tStep) = 0;
    /// Calculates critical time step.
    virtual double computeCriticalTimeStep(TimeStep *tStep) = 0;

    // time step termination
    virtual void updateInternalState(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual int checkConsistency();

    // definition
    virtual const char *giveClassName() const { return "CBSElement"; }
    virtual classType giveClassID() const { return CBSElementClass; }

    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *atTime);
    // Graphics output
    //void drawYourself(oofegGraphicContext&);
    //virtual void drawRawGeometry(oofegGraphicContext&) {}
    //virtual void drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

protected:
    /// Computes deviatoric stress vector in given integration point and solution step from given total strain vector
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) = 0;
};
} // end namespace oofem
#endif // cbselement_h
