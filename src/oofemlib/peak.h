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

#ifndef peak_h
#define peak_h

#include "loadtime.h"

namespace oofem {
/**
 * This class implements a function that is 0 everywhere, except in a single
 * point.
 */
class PeakFunction : public LoadTimeFunction
{
private:
    /// Specific time when function is nonzero.
    double t;
    /// Value of function at nonzero time.
    double value;

public:
    PeakFunction(int i, Domain *d) : LoadTimeFunction(i, d)
    {
        t = 0.0;
        value = 0.0;
    }
    virtual ~PeakFunction() { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual classType giveClassID() const { return PeakFunctionClass; }
    virtual const char *giveClassName() const { return "PeakFunction"; }

    virtual double  __at(double);
};
} // end namespace oofem
#endif // peak_h
