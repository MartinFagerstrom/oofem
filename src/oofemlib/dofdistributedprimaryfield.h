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

#ifndef dofdistributedprimaryfield_h
#define dofdistributedprimaryfield_h

#include "primaryfield.h"

namespace oofem {
/**
 * Class representing field of primary variables, which are typically allocated on nodes.
 * The field is determined by DOF values stored in DOF repositories (unknown dictionary).
 * These repositories are maintained and updated by engineering models since the algorithms are very specific to each model.
 * The class can return several variables stored in DOF. The purpose of this class is to provide
 * a shell that allows to access these repositories using field services.
 * The class contains also a solution vector for temporal storage of unknowns. The vector needs to be projected back to DOFs.
 */
class DofDistributedPrimaryField : public PrimaryField
{
public:
    /**
     * Constructor. Creates a field of given type associated to given domain.
     * Not using pointer to domain, because this will prevent the use of PrimaryField as an
     * EngngModel attribute. This is because the domain does not exists when
     * PrimaryField is created (this is when EngngModel is created).
     * @param a Engineering model which field belongs to.
     * @param idomain Index of domain for field.
     * @param ft Type of stored field.
     * @param ut Equation ID for unknowns in field.
     * @param nHist Number of old time steps to store.
     */
    DofDistributedPrimaryField(EngngModel *a, int idomain, FieldType ft, EquationID ut, int nHist);
    ~DofDistributedPrimaryField();

    virtual void initialize(ValueModeType mode, TimeStep *atTime, FloatArray &answer);

    virtual double giveUnknownValue(Dof *dof, ValueModeType mode, TimeStep *atTime);

    virtual void update(ValueModeType mode, TimeStep *atTime, FloatArray &vectorToStore);

    virtual FloatArray *giveSolutionVector(TimeStep *atTime);

    virtual void advanceSolution(TimeStep *atTime);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode);
};
} // end namespace oofem
#endif // dofdistributedprimaryfield_h
