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

#ifndef lumpedmasselement_h
#define lumpedmasselement_h

#include "structuralelement.h"

namespace oofem {

/**
 * This class implements a simple lumped mass element. Its purpose is to introduce
 * an additional mass (mass components or rotary inertias) into a node.
 * The mass element is defined by a single node.
 * At present, mass is defined in the nodal coordinate system.
 * The same element can be used to add an additional stiffness if needed (Not yet implemented).
 */
class LumpedMassElement : public StructuralElement
{
protected:
    ///Mass and moments of inertia corresponding to nodal DOFs
    FloatArray components;

public:
    LumpedMassElement(int n, Domain *d);
    virtual ~LumpedMassElement() { }

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
    { answer.resize(0, 0); }
    virtual void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
    { answer.resize(0, 0); }
    virtual void computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType)
    { answer.resize(0); }
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0)
    { answer.resize(0); }

    virtual int computeNumberOfDofs(EquationID ut);
    virtual void giveDofManDofIDMask(int inode, EquationID eid, IntArray &answer) const;

    virtual void updateInternalState(TimeStep *tStep) {}
    virtual void updateYourself(TimeStep *tStep) {}
    virtual int checkConsistency();

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
    void drawDeformedGeometry(oofegGraphicContext &, UnknownType);
    void drawScalar(oofegGraphicContext &context);
#endif

    // definition & identification
    virtual const char *giveClassName() const { return "LumpedMassElement"; }
    virtual classType giveClassID() const { return LumpedMassElementClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_point; }

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                          int lowerIndx = 1, int upperIndx = ALL_STRAINS)
    {}
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer) {}
};
} // end namespace oofem
#endif // lumpedmasselement_h
