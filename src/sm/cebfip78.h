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

#ifndef cebfip78_h
#define cebfip78_h

#include "maxwellChM.h"

namespace oofem {

/**
 * This class implements a CEB-FIP 78 rheologic Maxwell chain model in a finite
 * element problem.
 */
class CebFip78Material : public MaxwellChainMaterial
{
protected:
    double E28;    ///< Young modulus at age of 28 days [MPa].
    double fibf;   ///< Basic creep coefficient.
    double kap_a;  ///< Coefficient of hydrometric conditions.
    double kap_c;  ///< Coefficient of type of cement.
    double kap_tt; ///< Coefficient of temperature effects.
    double u;      ///< Surface imposed to environment [mm^2]; temporary here ; should be in crosssection level.

public:
    CebFip78Material(int n, Domain *d) : MaxwellChainMaterial(n, d) { }
    virtual ~CebFip78Material() { }

    virtual const char *giveClassName() const { return "CebFip78Material"; }
    virtual classType giveClassID() const { return CebFip78MaterialClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);

protected:
    virtual double computeCreepFunction(GaussPoint *gp, double atTime, double ofAge);
};
} // end namespace oofem
#endif // cebfip78_h
