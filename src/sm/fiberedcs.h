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

#ifndef fiberedcs_h
#define fiberedcs_h


#include "structuralcrosssection.h"
#include "structuralmaterial.h"
#include "element.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "interface.h"

namespace oofem {
class GaussPoint;
class FiberedCrossSectionModelInterface;

/**
 * This class implements a fibered cross section in a finite element problem. A cross
 * section  is an attribute of a domain. It is usually also attribute of many
 * elements.
 *
 * The attribute 'propertyDictionary' contains all the properties of a
 * layered cross section, like thickness and width of each layer.
 * The attribute 'layerMaterials' contains an array of Materials corresponding
 * to each layer.
 *
 * it uses master - slave GaussPoint approach, where master gp has more slaves gp.
 * slave gp represent for each fiber material point. It's coordinate sections
 * contains y,z-coordinates from mid-section. the slaves are manageg completely
 * ( created, saved their context.,,,) from this class. Master gp only deletes
 * slaves in destructor.
 *
 * Tasks:
 * - Returning standard material stiffness marices (like 3d stress-strain, 2d plane ,
 *   plate, 3dbeam, 2d beam ..) according to current state determined by parameter
 *   StressMode by calling gp->material->GiveMaterialStiffnessMatrix (....) and by
 *   possible modifying returned matrix. (for example in layered mode approach
 *   each layer  is asked for 3dMatrialStiffnes and this is integrated for example
 *   over thickness for plate bending problems)
 * - Returning RealStress state in gauss point and for given Stress mode.
 * - Returning a properties of cross section like thickness or area.
 */
class FiberedCrossSection : public StructuralCrossSection
{
protected:
    IntArray fiberMaterials; ///< Material of each fiber.
    FloatArray fiberThicks; ///< Thickness for each fiber.
    FloatArray fiberWidths; ///< Width for each fiber.
    int numberOfFibers;     ///< Number of fibers.
    double thick; ///< Total thickness.
    double width; ///< Total width.
    double area;  ///< Total area.
    FloatArray fiberYcoords, fiberZcoords;
public:

    FiberedCrossSection(int n, Domain *d) : StructuralCrossSection(n, d), fiberMaterials(), fiberThicks(), fiberWidths(),
        fiberYcoords(), fiberZcoords()
    {
        thick = 0.;
        width = 0.;
        area = -1.0;
    }

    virtual ~FiberedCrossSection()  { }

    virtual void giveRealStresses(FloatArray & answer, MatResponseForm, GaussPoint *,
                          const FloatArray &, TimeStep * tStep);

    virtual void giveCharMaterialStiffnessMatrix(FloatMatrix &answer,
                                         MatResponseMode rMode,
                                         GaussPoint *,
                                         TimeStep *tStep);

    virtual void giveCharMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                   MatResponseForm form, MatResponseMode rMode,
                                                   GaussPoint *, StructuralMaterial *,
                                                   TimeStep *tStep);



    virtual void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                            const FloatArray &charVector3d);
    virtual void giveFullCharacteristicVector(FloatArray &answer,
                                      GaussPoint *gp, const FloatArray &charVector);
    virtual FloatArray *imposeStressConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStressVector3d);
    virtual FloatArray *imposeStrainConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStrainVector3d);
    virtual void giveStressStrainMask(IntArray &answer, MatResponseForm form,
                                      MaterialMode mmode, StructuralMaterial *mat) const;

    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode, Material *mat);
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);

    virtual void giveFiberMaterialStiffnessMatrix(FloatMatrix &fiberMatrix, MatResponseForm FullForm,
                                                  MatResponseMode rMode, GaussPoint *layerGp,
                                                  TimeStep *tStep);

    virtual void computeStressIndependentStrainVector(FloatArray &answer,
                                                      GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    virtual double give(CrossSectionProperty a);

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "FiberedCrossSection"; }
    virtual classType giveClassID() const { return FiberedCrossSectionClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void printYourself();
    double computeIntegralThickWidth();
    MaterialMode giveCorrespondingSlaveMaterialMode(MaterialMode);
    GaussPoint *giveSlaveGaussPoint(GaussPoint *gp, int);

    virtual contextIOResultType saveIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);
    virtual contextIOResultType restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);

#ifdef __PARALLEL_MODE
    int packUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
    {
        _error("packUnknowns: not implemented");
        return 0;
    }

    int unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN, GaussPoint *ip)
    {
        _error("unpackAndUpdateUnknowns: not implemented");
        return 0;
    }

    int estimatePackSize(CommunicationBuffer &buff, GaussPoint *ip)
    {
        _error("estimatePackSize: not implemented");
        return 0;
    }
#endif

protected:
    virtual void giveMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                               MatResponseForm form,
                                               MatResponseMode rMode,
                                               GaussPoint *gp,
                                               StructuralMaterial *mat,
                                               TimeStep *tStep);

    void give3dBeamMaterialStiffnessMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode rMode,
                                           GaussPoint *gp,
                                           StructuralMaterial *mat,
                                           TimeStep *tStep);

    FloatArray *GiveIntegrated3dBeamStress(GaussPoint *gp);

    double giveArea();

    friend class Material;
};

/**
 * The element interface required by FiberedCrossSection.
 */
class FiberedCrossSectionInterface : public Interface
{
public:
    FiberedCrossSectionInterface() { }

    /**
     * Computes full 3d strain vector in element fiber. This function is necessary
     * if layered cross section is specified.
     * @param answer Full fiber strain vector.
     * @param masterGp Element integration point.
     * @param slaveGp Slave integration point representing particular fiber.
     * @param tStep Time step.
     */
    virtual void FiberedCrossSectionInterface_computeStrainVectorInFiber(FloatArray &answer, GaussPoint *masterGp,
                                                                         GaussPoint *slaveGp, TimeStep *tStep) = 0;
};
} // end namespace oofem
#endif // fiberedcs_h
