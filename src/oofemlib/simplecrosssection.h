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

#ifndef simplecrosssection_h
#define simplecrosssection_h

#include "structuralcrosssection.h"
#include "structuralmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
class GaussPoint;

/**
 * Class implementing "simple" cross section model in finite element problem.
 * A cross section  is an attribute of a domain. It is usually also attribute of many
 * elements.
 *
 * The simple cross section implementation does not perform any integration over cross-section volume,
 * it represents a cross section model, where the whole cross section model is represented by single integration point.
 * and therefore all requests for characteristic contributions (stiffness) and for real stress computations are simply
 * passed to parent StructuralCrossSection class, which invokes corresponding material mode services.
 * Please note, that it is assumed that material model will support these material modes and provide
 * corresponding services for characteristic components and stress evaluation.
 * For description, how to incorporate more elaborate models of cross section, please read
 * base CrossSection documentation.
 *
 * The overloaded methods giveFullCharacteristicVector and giveFullCharacteristicVector add some additional support
 * for integrated cross section models - _3dShell, _3dBeam, _2dPlate and _2dBeam.
 *
 * This class also reads into its property dictionary necessary geometric cross section characteristics,
 * which are requested by particular material models. These parameters can be requested using get service and
 * include those defined by CrossSectionProperty.
 */
class SimpleCrossSection : public StructuralCrossSection
{
public:
    /** 
     * Constructor.
     * @param n Cross section number.
     * @param d Associated domain.
     */
    SimpleCrossSection(int n, Domain *d) : StructuralCrossSection(n, d) { }

    virtual void giveFullCharacteristicVector(FloatArray &answer,  GaussPoint *gp, const FloatArray &strainVector);
    virtual void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                          const FloatArray &charVector3d);
    virtual void giveRealStresses(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                          const FloatArray &reducedStrainIncrement, TimeStep *tStep);

    virtual void giveCharMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                   MatResponseForm form, MatResponseMode mode,
                                                   GaussPoint *gp, StructuralMaterial *mat,
                                                   TimeStep *tStep);

    virtual void computeStressIndependentStrainVector(FloatArray &answer,
                                                      GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    virtual double give(CrossSectionProperty a);

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "SimpleCrossSection"; }
    virtual classType giveClassID() const { return SimpleCrossSectionClass; }
    virtual const char *giveInputRecordName() const { return "SimpleCS"; }

    /**
     * Initializes receiver acording to object description stored in input record.
     * Calls CrossSection initializeFrom service and reads the values of
     * - 'thick' thickness
     * - 'width' width
     * - 'area' area
     * - 'iy' Moment of inertia around y
     * - 'iz' Moment of inertia around z
     * - 'ik' Torsion moment around x
     * - 'beamshearcoeff' Beam shear coefficient
     * @param ir Record to read off.
     */
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

protected:
    virtual void giveMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                       MatResponseForm form,
                                       MatResponseMode rMode,
                                       GaussPoint *gp,
                                       StructuralMaterial *mat,
                                       TimeStep *tStep);
};
} // end namespace oofem
#endif // simplecrosssection_h
