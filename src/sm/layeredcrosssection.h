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

#ifndef layeredcrosssection_h
#define layeredcrosssection_h

#include "structuralcrosssection.h"
#include "structuralmaterial.h"
#include "element.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "interface.h"

namespace oofem {
class GaussPoint;
class LayeredCrossSectionModelInterface;

/**
 * This class implements a layered cross section in a finite element problem. A cross
 * section  is an attribute of a domain. It is usually also attribute of many
 * elements.
 *
 * The attribute 'propertyDictionary' contains all the properties of a
 * layered cross section, like thickness and width of each layer.
 * The attribute 'layerMaterials' contains an array of Materials corresponding
 * to each layer.
 *
 * It uses master - slave GaussPoint approach, where master gp has more slaves gp.
 * slave gp represent for each layer material point. It's coordinate sections
 * contains z-coordinate (-1,1) from mid-section. the slaves are manageg completely
 * ( created, saved their context.,,,) from this class. Master gp only deletes
 * slaves in destructor.
 *
 * Tasks:
 * - Returning standard material stiffness marices (like 3dstress-strain, 2d plane ,
 *   plate, 3dbeam, 2d beam ..) according to current state determined by parametr
 *   StressMode by calling gp->material->GiveMaterialStiffnessMatrix (....) and by
 *   possible modifiing returned matrix. (for example in layerde mode aproach
 *   each layer  is asked for 3dMatrialStiffnes and this is integrated for example
 *   over thickness for plate bending broblems)
 * - Returning RealStress state in gauss point and for given Stress mode.
 * - Returning a properties of cross section like thickness or area.
 */
class LayeredCrossSection : public StructuralCrossSection
{
protected:
    IntArray layerMaterials; ///< Material of each layer.
    FloatArray layerThicks; ///< Thickness for each layer.
    FloatArray layerWidths; ///< Width for each layer.
    int numberOfLayers;
    double midSurfaceZcoordFromBottom, totalThick;
    double area;

public:
    LayeredCrossSection(int n, Domain *d) : StructuralCrossSection(n, d), layerMaterials(), layerThicks(), layerWidths()
    {
        numberOfLayers = 0;
        totalThick = 0.;
        area = -1.0;
    }

    virtual ~LayeredCrossSection() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveRealStresses(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                          const FloatArray &reducedStrainIncrement, TimeStep *tStep);

    virtual void giveCharMaterialStiffnessMatrix(FloatMatrix &answer,
                                         MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *tStep);

    // next function is intended to be used if we would like to obtain
    // char matrix form different material which is not associated with gp and its element.
    // (mainly for obtaining linear elastic matrix)
    // stress-strain mode is taken from gp.
    // NORMALLY - PLEASE USE GiveCharMaterialStiffnessMatrix function
    virtual void giveCharMaterialStiffnessMatrixOf(FloatMatrix &answer,
                                                   MatResponseForm form, MatResponseMode rMode,
                                                   GaussPoint *, StructuralMaterial *,
                                                   TimeStep *tStep);



    virtual void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                         const FloatArray &charVector3d);
    virtual void giveFullCharacteristicVector(FloatArray &answer,
                                      GaussPoint *gp, const FloatArray &strainVector);

    virtual FloatArray *imposeStressConstrainsOnGradient(GaussPoint *gp, FloatArray *);
    virtual FloatArray *imposeStrainConstrainsOnGradient(GaussPoint *gp, FloatArray *);

    virtual void giveStressStrainMask(IntArray &answer, MatResponseForm form,
                                      MaterialMode mmode, StructuralMaterial *mat) const;
    virtual void giveLayerMaterialStiffnessMatrix(FloatMatrix &layerMatrix, MatResponseForm FullForm,
                                                  MatResponseMode rMode, GaussPoint *layerGp,
                                                  TimeStep *tStep);

    virtual void computeStressIndependentStrainVector(FloatArray &answer,
                                                      GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    virtual double give(CrossSectionProperty a);

    /// Returns the total thickness of all layers.
    double computeIntegralThick();

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "LayeredCrossSection"; }
    virtual classType giveClassID() const { return LayeredCrossSectionClass; }
    virtual void printYourself();

    MaterialMode giveCorrespondingSlaveMaterialMode(MaterialMode);
    GaussPoint *giveSlaveGaussPoint(GaussPoint *gp, int slaveIndex);

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
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               StructuralMaterial *mat,
                                               TimeStep *tStep);
    void giveDerivedMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseForm form,
                                            MatResponseMode mode,
                                            GaussPoint *, StructuralMaterial *mat,
                                            TimeStep *tStep);

    void give2dPlateMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseForm form,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            StructuralMaterial *mat,
                                            TimeStep *tStep);
    void give3dShellMaterialStiffness(FloatMatrix &answer,
                                      MatResponseForm form,
                                      MatResponseMode mode,
                                      GaussPoint *gp,
                                      StructuralMaterial *mat,
                                      TimeStep *tStep);
    void give2dBeamMaterialStiffnessMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           StructuralMaterial *mat,
                                           TimeStep *tStep);

    void giveIntegrated3dShellStress(FloatArray &answer, GaussPoint *gp);

    double giveArea();

    friend class Material;
};

/**
 * The element interface required by LayeredCrossSection.
 */
class LayeredCrossSectionInterface : public Interface
{
public:
    LayeredCrossSectionInterface() { }

    /**
     * Computes full 3D strain vector in element layer. This function is necessary
     * if layered cross section is specified. If it is implemented, the testElementExtension
     * service should return nonzero for Element_LayeredSupport parameter. This service is used by
     * layered cross section models.
     * @param answer Full layer strain vector.
     * @param masterGp Element integration point.
     * @param slaveGp Slave integration point representing particular layer.
     * @param tStep Time step.
     */
    virtual void computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                            GaussPoint *slaveGp, TimeStep *tStep) = 0;
};
} // end namespace oofem
#endif // layeredcrosssection_h
