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


//   *****************************************
//   *** CLASS ERROR ESTIMATOR (INDICATOR) ***
//   *****************************************

#ifndef errorestimator_h
#define errorestimator_h

#include "femcmpnn.h"
#include "compiler.h"

#include "interface.h"
#include "remeshingcrit.h"
#include "classtype.h"
#include "errorestimatortype.h"
#include "intarray.h"

namespace oofem {
class Domain;
class Element;
class TimeStep;

/** Type characterizing different type of errors. */
enum EE_ValueType { globalNormEEV, globalErrorEEV, globalWeightedErrorEEV };
/** Type characterizing different type of element errors. */
enum EE_ErrorType { unknownET, indicatorET, internalStressET, primaryUnknownET };
/** Type determining whether temporary or equilibrated variables are used for error evaluation. */
enum EE_ErrorMode { equilibratedEM, temporaryEM };


/**
 * The base class for all error estimation or error indicator algorithms.
 * The basic task is to evaluate the error on associated domain. If this task requires
 * the special element algorithms, these should be included using interface concept.
 *
 * This estimator should also provide the compatible Remeshing Criteria class, which
 * based on various error measures will evaluate the required mesh density of a new domain.
 */
class ErrorEstimator : public FEMComponent
{
protected:
    ErrorEstimatorType eeType;
    RemeshingCriteria *rc;
    /**
     * Map indicating regions to skip (region - cross section model).
     * Do not access this variable directly, since this variable is read from input and could have size different
     * from actual number of regions - use always the skipRegion method, since it performs size check.  and handles
     * all regions correctly.
     */
    IntArray regionSkipMap;
    /// Number of skipped elements.
    int skippedNelems;

public:
    /// Constructor
    ErrorEstimator(int n, Domain *d) : FEMComponent(n, d) {
        rc = NULL;
        skippedNelems = 0;
        regionSkipMap.resize(0);
    }
    /// Destructor
    virtual ~ErrorEstimator() { if ( rc ) { delete rc; } }
    /// Sets Domain; should also re-initialize attributes if necessary.
    void setDomain(Domain *d);
    /**
     * Returns the element error. The estimateError service should be called before.
     * @param type Error type.
     * @param elem Element for which error requested.
     * @param tStep Time step.
     * @return Element error.
     */
    virtual double giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep) = 0;
    /**
     * Returns the characteristic value of given type. The estimateError service should be called before.
     * This method is supposed to be used by associated remeshingCriteria to access some characteristic
     * values already computed or known at error estimator level.
     * @param type Error type.
     * @param tStep Time step.
     * @return Error value for given type.
     */
    virtual double giveValue(EE_ValueType type, TimeStep *tStep) = 0;
    /**
     * Returns number of elements skipped in error estimation.
     * @return Number of skipped elements.
     */
    int giveNumberOfSkippedElements() { return skippedNelems; }
    /**
     * Estimates the error on associated domain at given time step. The estimated values can be
     * requested using giveElementError and giveValue methods. The type of errors provided
     * depends on error estimator type implementing the service.
     * @param mode Error mode.
     * @param tStep Time step.
     */
    virtual int estimateError(EE_ErrorMode mode, TimeStep *tStep) = 0;
    /**
     * Returns reference to associated remeshing criteria.
     */
    virtual RemeshingCriteria *giveRemeshingCrit() = 0;

    /**
     * Returns error estimation type of receiver.
     */
    ErrorEstimatorType giveErrorEstimatorType() const { return eeType; }

    /**
     * Returns nonzero if region has been skipped in error estimation (user option).
     * It is strongly recommended to use this function, instead of direct access to
     * regionSkipMap variable by derived classes, since the size check is done here.
     * @param reg Region to check.
     * @return True if region should be skipped.
     */
    bool skipRegion(int reg) { if ( reg <= regionSkipMap.giveSize() ) { return regionSkipMap.at(reg) > 0; } else { return false; } }
    virtual void reinitialize() { this->rc->reinitialize(); }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "ErrorEstimator"; }
    virtual classType giveClassID() const { return ErrorEstimatorClass; }
};
} // end namespace oofem
#endif // errorestimator_h






