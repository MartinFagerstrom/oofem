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

#ifndef igaelements_h
#define igaelements_h

#include "iga.h"
#include "planestresselementevaluator.h"
#include "space3delementevaluator.h"
#include "feibspline.h"
#include "feinurbs.h"
#include "feitspline.h"

#include "flotarry.h"
#include "flotmtrx.h"
#include "matresponsemode.h"
#include "mathfem.h"

namespace oofem {
class BsplinePlaneStressElement : public IGAElement, public PlaneStressStructuralElementEvaluator
{
protected:
    BSplineInterpolation interpolation;

public:
    BsplinePlaneStressElement(int n, Domain *aDomain);
    IRResultType initializeFrom(InputRecord *ir);

    void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicVector(answer, type, mode, t);
    }

    FEInterpolation *giveInterpolation() { return & this->interpolation; }
    virtual Element *giveElement() { return this; }
    void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const {
        PlaneStressStructuralElementEvaluator :: giveDofManDofIDMask(inode, u, answer);
    }
    virtual int computeNumberOfDofs(EquationID ut) { return numberOfDofMans * 2; }
    void updateInternalState(TimeStep *stepN) { PlaneStressStructuralElementEvaluator :: updateInternalState(stepN); }
    // definition & identification
    virtual const char *giveClassName() const { return "BsplinePlaneStressElement"; }
    virtual classType giveClassID() const { return BsplinePlaneStressElementClass; }

#ifdef __OOFEG
    //
    // Graphics output
    //
    virtual void  drawScalar(oofegGraphicContext &context);
    virtual void drawDeformedGeometry(oofegGraphicContext &mode, UnknownType ut) {
        drawIGAPatchDeformedGeometry(this, this, mode, ut);
    }
#endif

protected:
    virtual int giveNsd() { return 2; }
};


class NURBSPlaneStressElement : public IGAElement, public PlaneStressStructuralElementEvaluator
{
protected:
    NURBSInterpolation interpolation;

public:
    NURBSPlaneStressElement(int n, Domain *aDomain);
    IRResultType initializeFrom(InputRecord *ir);

    void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicVector(answer, type, mode, t);
    }

    FEInterpolation *giveInterpolation() { return & this->interpolation; }
    virtual Element *giveElement() { return this; }
    void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const {
        PlaneStressStructuralElementEvaluator :: giveDofManDofIDMask(inode, u, answer);
    }
    virtual int computeNumberOfDofs(EquationID ut) { return numberOfDofMans * 2; }
    void updateInternalState(TimeStep *stepN) { PlaneStressStructuralElementEvaluator :: updateInternalState(stepN); }
    // definition & identification
    virtual const char *giveClassName() const { return "NURBSPlaneStressElement"; }
    virtual classType giveClassID() const { return NURBSPlaneStressElementClass; }
#ifdef __OOFEG
    //
    // Graphics output
    //
    virtual void  drawScalar(oofegGraphicContext &context);
    virtual void drawDeformedGeometry(oofegGraphicContext &mode, UnknownType ut) {
        drawIGAPatchDeformedGeometry(this, this, mode, ut);
    }

#endif

protected:
    virtual int giveNsd() { return 2; }
};



class TSplinePlaneStressElement : public IGATSplineElement, public PlaneStressStructuralElementEvaluator
{
protected:
    TSplineInterpolation interpolation;

public:
    TSplinePlaneStressElement(int n, Domain *aDomain);
    IRResultType initializeFrom(InputRecord *ir) {
        IGATSplineElement :: initializeFrom(ir);
        //PlaneStressStructuralElementEvaluator::initializeFrom(ir);
        return IRRT_OK;
    }

    void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicVector(answer, type, mode, t);
    }

    FEInterpolation *giveInterpolation() { return & this->interpolation; }
    virtual Element *giveElement() { return this; }
    void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const {
        PlaneStressStructuralElementEvaluator :: giveDofManDofIDMask(inode, u, answer);
    }
    virtual int computeNumberOfDofs(EquationID ut) { return numberOfDofMans * 2; }
    void updateInternalState(TimeStep *stepN) { PlaneStressStructuralElementEvaluator :: updateInternalState(stepN); }
    // definition & identification
    virtual const char *giveClassName() const { return "TSplinePlaneStressElement"; }
    virtual classType giveClassID() const { return TSplinePlaneStressElementClass; }
#ifdef __OOFEG
    // Graphics output
    virtual void  drawScalar(oofegGraphicContext &context);
#endif

protected:
    virtual int giveNsd() { return 2; }
};


class NURBSSpace3dElement : public IGAElement, public Space3dStructuralElementEvaluator
{
protected:
    NURBSInterpolation interpolation;

public:
    NURBSSpace3dElement(int n, Domain *aDomain);
    IRResultType initializeFrom(InputRecord *ir);

    void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) {
        Space3dStructuralElementEvaluator :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) {
        Space3dStructuralElementEvaluator :: giveCharacteristicVector(answer, type, mode, t);
    }

    FEInterpolation *giveInterpolation() { return & this->interpolation; }
    virtual Element *giveElement() { return this; }
    void giveDofManDofIDMask(int inode, EquationID u, IntArray &answer) const {
        Space3dStructuralElementEvaluator :: giveDofManDofIDMask(inode, u, answer);
    }
    virtual int computeNumberOfDofs(EquationID ut) { return numberOfDofMans * 3; }
    void updateInternalState(TimeStep *stepN) { Space3dStructuralElementEvaluator :: updateInternalState(stepN); }
    // definition & identification
    virtual const char *giveClassName() const { return "NURBSSpace3dElement"; }
    virtual classType giveClassID() const { return NURBSSpace3dElementClass; }
#ifdef __OOFEG
    // Graphics output
    virtual void drawScalar(oofegGraphicContext &context);
    virtual void drawDeformedGeometry(oofegGraphicContext &mode, UnknownType ut) {
        drawIGAPatchDeformedGeometry(this, this, mode, ut);
    }
#endif

protected:
    virtual int giveNsd() { return 3; }
};
} // end namespace oofem
#endif //igaelements_h
