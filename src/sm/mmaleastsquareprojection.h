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

#ifndef mmaleastsquareprojection_h
#define mmaleastsquareprojection_h

#include "alist.h"
#include "materialmappingalgorithm.h"
#include "interface.h"
#include "compiler.h"

#include "dynalist.h"

namespace oofem {
class Domain;
class Element;
class TimeStep;

enum MMALeastSquareProjectionPatchType { MMALSPPatchType_1dq, MMALSPPatchType_2dq };
/**
 * Defines, whether only the necessary number of closest points will be used
 * to fit a polynomial. If not defined, more points from gp neighborhood will
 * be used, based on element connectivity.
 */
//#define MMALSP_ONLY_CLOSEST_POINTS

/**
 * The class implements the transfer of state variables based on
 * Least square fit over old mesh integration points (in the neighborhood of point of interest).
 *
 * The algorithm of projecting internal vars (q) can be summarized as follows:
 * -# The "source" element on old mesh containing point of interest is found
 * -# The element patch is constructed from "source" element neighbourhood
 * -# The least square fit is done
 * -# Value in point of interest is evaluated.
 *
 * It is obvious, that this mapper operates locally and therefore there is no need to declared this
 * mapper as static material member.
 *
 */
class MMALeastSquareProjection : public MaterialMappingAlgorithm
{
protected:
    /// If set, then only IP in the neighbourhood with same state can be used to interpolate the values.
    int stateFilter;
    /// If set, then only IP in the same region are taken into account.
    int regionFilter;
    /// List of Gp participating in patch.
    dynaList< GaussPoint * >patchGPList;
    /// Patch domain.
    Domain *patchDomain;
    /// Type of patch.
    MMALeastSquareProjectionPatchType patchType;
public:
    /// Constructor
    MMALeastSquareProjection();
    /// Destructor
    virtual ~MMALeastSquareProjection();

    virtual void __init(Domain *dold, IntArray &type, FloatArray &coords, int region, TimeStep *tStep);

    virtual void finish(TimeStep *tStep);

    virtual int __mapVariable(FloatArray &answer, FloatArray &coords, InternalStateType type, TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    virtual const char *giveClassName() const { return "MMALeastSquareProjectionPatchType"; }

protected:
    void computePolynomialTerms(FloatArray &P, FloatArray &coords, MMALeastSquareProjectionPatchType type);
    int giveNumberOfUnknownPolynomialCoefficients(MMALeastSquareProjectionPatchType regType);
};
} // end namespace oofem
#endif // mmaleastsquareprojection_h
