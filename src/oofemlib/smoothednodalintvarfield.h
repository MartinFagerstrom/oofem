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

#ifndef smoothednodalintvarfield_h
#define smoothednodalintvarfield_h

#include "domain.h"
#include "field.h"
#include "zznodalrecoverymodel.h"

namespace oofem {
/**
 * Class representing a field of an internal variable smoothed from integration points
 * into nodes. It is relying on NodalRecoveryModel to smooth the data.
 */
class SmoothedNodalInternalVariableField : public Field
{
protected:
    /// Smoother type
    NodalRecoveryModel::NodalRecoveryModelType stype;
    /// Smoother
    NodalRecoveryModel *smoother;
    /// InternalStateType.
    InternalStateType istType;
    /// Source domain.
    Domain *domain;

public:
    /**
     * Constructor. Creates a internal variable field of given type associated to given domain.
     * @param ist Physical meaning of field.
     * @param b Field type.
     * @param st detrmines the type of nodal recovery model used.
     * @param d Domain which field belongs to.
     */
    SmoothedNodalInternalVariableField(InternalStateType ist, FieldType b, NodalRecoveryModel::NodalRecoveryModelType st, Domain *d);
    ~SmoothedNodalInternalVariableField();

    virtual int evaluateAt(FloatArray &answer, FloatArray &coords, ValueModeType mode, TimeStep *atTime);
    virtual int evaluateAt(FloatArray &answer, DofManager* dman, ValueModeType mode, TimeStep *atTime);

    InternalStateType giveType() { return istType; }
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode);
    virtual const char *giveClassName() const { return "SmoothedNodalInternalVariableField"; }
};
} // end namespace oofem
#endif // smoothednodalintvarfieldintvarfield_h
