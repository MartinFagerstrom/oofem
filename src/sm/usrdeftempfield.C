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

#include "usrdeftempfield.h"
#include "timestep.h"
#include <sstream>

namespace oofem {
void
UserDefinedTemperatureField :: computeValueAt(FloatArray &answer, TimeStep *stepN, FloatArray &coords, ValueModeType mode)
// Returns the value of the receiver at time and given position respecting the mode.
{
    int i, err;
    double result;

    if ( ( mode != VM_Incremental ) && ( mode != VM_Total ) ) {
        _error2( "computeComponentArrayAt: unknown mode (%s)", __ValueModeTypeToString(mode) );
    }

    answer.resize(this->size);
    std::ostringstream buff;
    for ( i = 1; i <= size; i++ ) {
        buff << "x=" << coords.at(1) << ";y=" << coords.at(2) << ";z=" << coords.at(3) <<
            ";t=" << stepN->giveTargetTime() << ";" <<ftExpression[i-1];
        result = myParser.eval(buff.str().c_str(), err);
        if ( err ) {
            _error("computeValueAt: parser syntax error");
        }

        answer.at(i) = result;

        if ( ( mode == VM_Incremental ) && ( !stepN->isTheFirstStep() ) ) {
            buff << "x=" << coords.at(1) << ";y=" << coords.at(2) << ";z=" << coords.at(3) <<
                ";t=" << (stepN->giveTargetTime() - stepN->giveTimeIncrement()) << ";" << ftExpression[i-1];
            result = myParser.eval(buff.str().c_str(), err);
            if ( err ) {
                _error("computeValueAt: parser syntax error");
            }

            answer.at(i) -= result;
        }
    }
}

IRResultType
UserDefinedTemperatureField :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, size, IFT_UserDefinedTemperatureField_size, "size"); // Macro
    if ( size > 3 ) {
        size = 3;
    }

    if ( size > 0 ) {
        IR_GIVE_FIELD(ir, ftExpression [ 0 ], IFT_UserDefinedTemperatureField_t1, "t1(txyz)");
    }

    if ( size > 1 ) {
        IR_GIVE_FIELD(ir, ftExpression [ 1 ], IFT_UserDefinedTemperatureField_t1, "t2(txyz)");
    }

    if ( size > 2 ) {
        IR_GIVE_FIELD(ir, ftExpression [ 2 ], IFT_UserDefinedTemperatureField_t1, "t3(txyz)");
    }

    return IRRT_OK;
}
} // end namespace oofem
