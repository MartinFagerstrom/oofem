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

#ifndef m4_h
#define m4_h

#include "microplanematerial_bazant.h"
#include "structuralms.h"

namespace oofem {
/**
 * Related material model status to BazantMBCMaterial class
 * for storing history  variables in particular integration point
 * (microplane). Unique copy for each microplane must exist.
 */
class M4MaterialStatus : public StructuralMaterialStatus
{
public:
    M4MaterialStatus(int n, Domain *d, GaussPoint *g);
    virtual ~M4MaterialStatus();

    // definition
    virtual const char *giveClassName() const { return "M4MaterialStatus"; }
    virtual classType giveClassID() const { return M4MaterialStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};


/**
 * Implementation of microplane material model according to Bazant's boundary curve
 * approach.
 */
class M4Material : public MicroplaneMaterial_Bazant
{
protected:

    double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12;
    double c13, c14, c15, c16, c17, c18, c19, c20; /*c... fixed empirical constants*/
    double k1, k2, k3, k4, k5, E, nu, mu;
    double EV, ED, ET;
    double talpha;

public:
    /**
     * Constructor. Creates  Bazant's Boundary Curve Microplane Material belonging
     * to domain d, with number n.
     * @param n Material number.
     * @param d Domain to which newly created material belongs.
     */
    M4Material(int n, Domain *d);
    /// Destructor.
    virtual ~M4Material() { }

    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual void giveRealMicroplaneStressVector(FloatArray &answer, Microplane *mplane, const FloatArray &strain, TimeStep *tStep);

    double macbra(double x);
    double FVplus(double ev, double k1, double c13, double c14, double c15, double Ev);
    double FVminus(double ev, double k1, double k3, double k4, double E);
    double FDminus(double ed, double k1, double c7, double c8, double c9, double E);
    double FDplus(double ed, double k1, double c5, double c6, double c7, double c20, double E);
    double FN(double en, double sv, double k1, double c1, double c2, double c3, double c4,
              double E, double Ev);
    double FT(double sn, double ev, double k1, double k2, double c10,
              double c11, double c12, double Et);

    virtual void updateVolumetricStressTo(Microplane *mPlane, double sigv);

    virtual int giveSizeOfReducedStressStrainVector(MaterialMode);
    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "M4Material"; }
    virtual classType giveClassID() const { return M4MaterialClass; }

protected:
    MaterialStatus *CreateMicroplaneStatus(GaussPoint *gp) { return new M4MaterialStatus(1, domain, gp); }
};
} // end namespace oofem
#endif // m4_h
