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

#ifndef xfemelementinterface_h
#define xfemelementinterface_h

#include "interface.h"
#include "alist.h"

namespace oofem {
class FloatArray;
class FloatMatrix;
class Triangle;
class Element;
class GaussPoint;

/**
 * Provides Xfem interface for an element.
 */
class XfemElementInterface : public Interface
{
public:
    /// Constructor.
    XfemElementInterface(Element *e) : Interface() { this->element = e; }
    /// Creates enriched part of B matrix.
    void XfemElementInterface_createEnrBmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    /// Partitions the element into patches by a triangulation.
    void XfemElementInterface_partitionElement(AList< Triangle > *answer, AList< FloatArray > *together);
    /// Updates integration rule based on the triangulation.
    void XfemElementInterface_updateIntegrationRule();
    /// Helpful routine to put the nodes for triangulation together, should be in protected members probably.
    void XfemElementInterface_prepareNodesForDelaunay(AList< FloatArray > *answer1, AList< FloatArray > *answer2);
protected:
    Element *element;
};
} // end namespace oofem
#endif // xfemelementinterface_h
