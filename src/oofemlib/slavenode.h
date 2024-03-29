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

#ifndef slavenode_h
#define slavenode_h

#include "node.h"

namespace oofem {

/**
 * Class implementing slave node connected to other nodes (masters) using predetermined weights.
 * Hanging node possess no degrees of freedom - all values are interpolated from corresponding master dofs.
 *
 * The contributions of hanging node are localized directly to master related equations.
 * The node can not have its own boundary or initial conditions,
 * they are determined completely from master dof conditions except for dofs of master type.
 * @see{HangingNode}
 */
class SlaveNode : public Node
{
protected:
    /// Master nodes for all dofs.
    IntArray masterDofManagers;
    /// Weights for each master node.
    FloatArray masterWeights;

public:
    /**
     * Constructor. Creates a hanging node with number n, belonging to aDomain.
     * @param n Node number in domain aDomain.
     * @param aDomain Domain to which node belongs.
     */
    SlaveNode(int n, Domain *aDomain) : Node(n, aDomain) { }
    /// Destructor.
    virtual ~SlaveNode(void) { }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int checkConsistency();
    virtual bool isDofTypeCompatible(dofType type) const { return ( type == DT_master || type == DT_slave ); }
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);

    virtual const char *giveClassName() const { return "SlaveNode"; }
    virtual classType giveClassID() const { return SlaveNodeClass; }
};
} // end namespace oofem
#endif // slavenode_h
