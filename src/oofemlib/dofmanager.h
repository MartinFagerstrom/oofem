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

#ifndef dofmanager_h
#define dofmanager_h

#ifndef __MAKEDEPEND
 #include <cstdio>
#endif

#include "alist.h"
#include "femcmpnn.h"
#include "intarray.h"
#include "equationid.h"
#include "valuemodetype.h"
#include "doftype.h"
#include "dofiditem.h"
#include "contextioresulttype.h"

namespace oofem {

class DataStream;
class Dof;
class Domain;
class EntityRenumberingFunctor;
class FloatMatrix;
class PrimaryField;
class UnknownNumberingScheme;

#ifdef __PARALLEL_MODE
class CommunicationBuffer;
/// In parallel mode, this type indicates the mode of DofManager.
enum dofManagerParallelMode {
    DofManager_local, /**< DofManager is local, there are no contribution from other domains to this DofManager.*/
    DofManager_shared, /**< DofManager is shared by neighboring partitions,
                            it is necessary to sum contributions from all contributing domains.
                            Typical for node cut algorithm. */
    DofManager_remote, /**< DofManager in active domain is only mirror of some remote DofManager.
                            It is necessary to copy remote values into local ones.
                            Typical for element cut.*/
    DofManager_null, /**< DofManager in active domain is shared only by remote elements (these are only
                          introduced for nonlocal constitutive model to allow effective local averaging, so only local material
                          value to be averaged are transferred for these remote elements). Null nodes are therefore used only
                          for computing real integration point coordinates of remote elements and there is no reason to maintain
                          their unknowns (they have no equation number assigned).*/
};

#endif

class NodalLoad;
class TimeStep;
class FloatArray;
class IntArray;

/**
 * Base class for dof managers. Dof manager is an abstraction for object possessing degrees of
 * freedom. Dof managers (respectively derived classes like nodes or sides) are usually attributes of
 * elements and are maintained by domain. Degrees of freedom belonging to dof manager are
 * stored in 'dofArray'. Dof manager also
 * maintain a list of applied loads it is subjected to.
 * Number and physical meaning of dofs can be specified by user in input file
 * (see input file description). If it is not specified, default values are obtained from domain,
 * based on domain type of problem.
 *
 * Tasks:
 * - managing its degrees of freedom (giveDof).
 * - calculating its nodal vector.
 * - printing and updating at end of step .
 * - managing its swapping to and from disk.
 * - managing the list of subjected loads (typically, concentrated forces and moments).
 * - constructing the location array for assembling loads (also used by elements for assembling matrices).
 */
class DofManager : public FEMComponent
{
protected:
    /// Total number of DOFs.
    int numberOfDofs;
    /// Array of DOFs.
    Dof **dofArray;
    /// List of applied loads.
    IntArray loadArray;
    /**
     * Indicates if dofManager is boundary (true boundary or on boundary between regions) or interior.
     * This information is required by some recovery techniques.
     */
    bool isBoundaryFlag;
    /// Flag indicating whether receiver has slave DOFs.
    bool hasSlaveDofs;
#if defined( __PARALLEL_MODE ) || defined( __ENABLE_COMPONENT_LABELS )
    /**
     * In parallel mode, globalNumber contains globally unique DoFManager number.
     * The component number, inherited from FEMComponent class contains
     * local domain number.
     */
    int globalNumber;
#endif
#ifdef __PARALLEL_MODE
    dofManagerParallelMode parallel_mode;
    /**
     * List of partition sharing the shared dof manager or
     * remote partition containing remote dofmanager counterpart.
     */
    IntArray partitions;
#endif

public:
    /**
     * Constructor. Creates DofManager with given number belonging to domain aDomain.
     * @param n DofManager's number in domain
     * @param aDomain reference to DofManager's domain
     */
    DofManager(int n, Domain *aDomain);
    /// Destructor.
    virtual ~DofManager();

    /**@name Dof management methods */
    //@{

    /**
     * Returns reference (pointer) to i-th dof of receiver.
     * Index of Dof with required physical meaning can be obtained by invoking
     * method findDofWithDofId.
     * @param i The index of the DOF.
     * @return The requested DOF.
     * @see DofManager::findDofWithDofId
     */
    Dof *giveDof(int i) const;
    /**
     * Returns DOF with given dofID; issues error if not present.
     * @param dofID The ID for the requested DOF.
     * @return The requested DOF.
     * @see DofIDItem
     */
    Dof *giveDofWithID(int dofID) const;
    /// @return Total number of dofs managed by receiver.
    int giveNumberOfDofs() const;

    /// Renumbers all contained DOFs.
    void askNewEquationNumbers(TimeStep *tStep);
    /**
     * Returns the number of primary dofs on which receiver dofs (given in dofArray) depend on.
     * If receiver has only primary dofs, the answer is the size of dofArray.
     * @param dofArray Array with indices to DOFs.
     * @return Number of primary dofs from given array.
     */
    int giveNumberOfPrimaryMasterDofs(IntArray &dofArray) const;
    /**
     * Returns location array (array containing for each requested dof related equation number) for
     * given numbering scheme.
     * @param dofIDArry Array containing dof mask. This mask containing DofIDItem values
     * (they describe physical meaning of dofs, see cltypes.h) is used to extract only
     * required values. If dof with requested physical meaning does not exist in receiver,
     * an error is generated and execution exits.
     * @param locationArray Return parameter containing required equation numbers.
     * @param s Determines the equation numbering scheme.
     * @see Element::giveDofManDofIDMask.
     */
    virtual void giveLocationArray(const IntArray &dofIDArry, IntArray &locationArray,
                                           const UnknownNumberingScheme &s) const;
    /**
     * Returns full location array of receiver containing equation numbers of all dofs
     * of receiver. Their order is specific to every DofManager. Mainly used at EngngModel level
     * to assemble DofManager contribution (typically load vector).
     * @param locationArray Complete location array of receiver.
     * @param s Determines the equation numbering scheme.
     */
    virtual void giveCompleteLocationArray(IntArray &locationArray, const UnknownNumberingScheme &s) const;
    /**
     * Returns DOFs numbers of receiver with required physical meaning.
     * @param dofIDArray Array containing DofIDItem-type values (this is enumeration
     * identifying physical meaning of particular DOF, see cltypes.h).
     * @param answer Array with DOF numbers. They are ordered according to dofIDArry.
     * @see DofIDItem
     * @see cltypes.h
     */
    void giveDofArray(const IntArray &dofIDArray, IntArray &answer) const;
    /**
     * Finds index of DOF with required physical meaning of receiver. This index can be different
     * for different DOFManagers (user can alter dof order and type in input file).
     * @param dofID Physical meaning of DOF.
     * @return Index of requested DOF. If such DOF with dofID does not exists, returns zero.
     */
    int findDofWithDofId(DofIDItem dofID) const;
    /**
     * Assembles the vector of unknowns in global c.s for given dofs of receiver.
     * @param answer Result (in global c.s.)
     * @param dofMask Array containing DOF mask. This mask containing DofIDItem values
     * (they describe physical meaning of dofs, see cltypes.h) is used to extract only
     * required values. If dof with requested physical meaning does not exist in receiver,
     * an error is generated and execution exits.
     * @param type Physical meaning of  unknown.
     * @param mode Mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @param stepN Time step when unknown requested. See documentation of particular EngngModel
     * class for valid stepN values (most implementation can return only values for current
     * and possibly for previous time step).
     * @see Dof::giveUnknown
     */
    virtual void giveUnknownVector(FloatArray &answer, const IntArray &dofMask,
                                   EquationID type, ValueModeType mode, TimeStep *stepN);
    /**
     * Assembles the vector of unknowns of given filed in global c.s for given dofs of receiver.
     * @param answer Result (in nodal cs.)
     * @param dofMask Array containing DOF mask. This mask containing DofIDItem values
     * (they describe physical meaning of dofs, see cltypes.h) is used to extract only
     * required values. If dof with requested physical meaning does not exist in receiver,
     * an error is generated and execution exits.
     * @param field Primary field.
     * @param mode Mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @param stepN Time step when unknown requested. See documentation of particular EngngModel
     * class for valid StepN values (most implementation can return only values for current
     * and possibly for previous time step).
     * @see Dof::giveUnknown
     */
    virtual void giveUnknownVector(FloatArray &answer, const IntArray &dofMask,
                                   PrimaryField &field, ValueModeType mode, TimeStep *stepN);

    /**
     * Assembles the vector of prescribed unknowns in nodal c.s for given dofs of receiver.
     * This vector may have size different from number of dofs requested,
     * because some dofs may depend on other dofs. Default implementation uses
     * Dof::giveBcValue and Dof::hasBc service.
     * @param answer Result (in nodal cs.)
     * @param dofMask Array containing dof mask. This mask containing DofIDItem values
     * (they describe physical meaning of dofs, see cltypes.h) is used to extract only
     * required values. If dof with requested physical meaning does not exist in receiver,
     * an error is generated and execution exits.
     * @param mode Mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @param stepN Time step when unknown requested. See documentation of particular EngngModel
     * class for valid StepN values (most implementation can return only values for current
     * and possibly for previous time step).
     * @see Dof::giveBcValue
     * @see Dof::hasBc
     * @todo Remove all usage of this. Just ask for the unknown vector instead, they are the same.
     */
    virtual void givePrescribedUnknownVector(FloatArray &answer, const IntArray &dofMask,
                                             ValueModeType mode, TimeStep *stepN);
    //@}

    /**@name Transformation functions
     * The governing equations can be assembled not only in global coordinate system, but
     * also in user-defined local coordinate system of each DOF manager. Methods in this section
     * introduce necessary transformation methods, allowing receiver DOFs to be expressed in
     * their own local c.s. or to be dependent on other DOFs on other DOF manager (to implement
     * slave or rigid arm nodes etc.). The method for computing global c.s to receiver c.s
     * transformation matrix is provided.
     */
    //@{
    /**
     * Computes receiver transformation matrix from global CS to dofManager specific
     * coordinate system; @f$ u_g = R\cdot u_m @f$.
     * @param answer Computed transformation matrix. It has generally dofIDArry.size rows and
     * if loc is obtained using giveLocationArray(dofIDArry, loc) call, loc.giveSize() columns.
     * This is because this transformation should generally include not only transformation to
     * dof manager local coordinate system, but receiver dofs can be expressed using
     * dofs of another dofManager (In this case, square answer is produced only if all
     * dof transformation is required).
     * @param dofIDArry Array containing DofIDItem-type values (this is enumeration
     * identifying physical meaning of particular DOF, see cltypes.h) for which transformation matrix is
     * assembled. if dofIDArry is NULL, then all receiver dofs are assumed.
     * @return True if transformation is needed, false otherwise.
     */
    bool computeM2GTransformation(FloatMatrix &answer, const IntArray &dofIDArry);
    /**
     * Computes transformation matrix from global c.s. to DOF-manager specific c.s; @f$ u_g = Q\cdot u_l @f$.
     * @param answer Computed transformation matrix.
     * @param dofIDArry Array containing DofIDItem-type values for which transformation matrix is
     * assembled. If dofIDArry is NULL, then all receiver DOFs are assumed.
     * @return True is transformation is needed, false otherwise.
     */
    virtual bool computeL2GTransformation(FloatMatrix &answer, const IntArray &dofIDArry);
    /**
     * Computes transformation matrix from local DOFs to master DOFs; @f$ u_l = M\cdot u_m @f$.
     * @param answer Computed transformation matrix.
     * @param dofIDArry Array containing DofIDItem-type values for which transformation matrix is
     * assembled. If dofIDArry is NULL, then all receiver DOFs are assumed.
     * @return True is transformation is needed, false otherwise.
     */
    virtual bool computeM2LTransformation(FloatMatrix &answer, const IntArray &dofIDArry);
    /**
     * Indicates, whether dofManager requires the transformation.
     * @return Nonzero if transformation is necessary, even for single dof.
     */
    virtual bool requiresTransformation();
    //@}

    /**@name Load related functions */
    //@{
    /**
     * Computes the load vector of receiver in given time.
     * @param answer Load vector.
     * @param stepN Time step when answer is computed.
     * @param mode Determines response mode.
     */
    virtual void computeLoadVectorAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    /**
     * Returns the array containing applied loadings of the receiver
     * @return Array with indices to the applied loads.
     */
    IntArray *giveLoadArray();
    /**
     * Sets the array of applied loadings of the receiver
     * @param load Array with indices to the applied loads.
     */
    void setLoadArray(IntArray &load);
    //@}

    /**@name Position query functions */
    //@{
    virtual bool hasCoordinates() { return false; }
    /// @return The i-th coordinate.
    virtual double giveCoordinate(int i) { return 0.0; }
    /// @return Pointer to node coordinate array.
    virtual FloatArray *giveCoordinates() { return NULL; }
    //@}

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    /**
     * Updates receiver after equilibrium in time step has been reached.
     * @param tStep Active time step.
     */
    void updateYourself(TimeStep *tStep);

    // Miscellaneous
    /// @return True if dofmanager is on boundary.
    bool isBoundary() { return isBoundaryFlag; }
    /**
     * Sets the boundary flag.
     * @param isBoundary Determines if receiver is on the boundary.
     */
    void setBoundaryFlag(bool isBoundary) { this->isBoundaryFlag = isBoundary; }
    /// @return True if receiver contains slave dofs
    virtual int hasAnySlaveDofs();
    /**
     * Returns true if the receiver is linked (its slave DOFs depend on master values) to some other dof managers.
     * In this case, the masters array should contain the list of masters.
     * In both serial and parallel modes, local numbers are be provided.
     * @param masters Indices of dof managers which receiver has slaves to.
     * @return If receiver contains only primary DOFs, false is returned.
     */
    virtual bool giveMasterDofMans(IntArray &masters);

    virtual const char *giveClassName() const { return "DofManager"; }
    virtual classType giveClassID() const { return DofManagerClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void printYourself();
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /// Returns true if dof of given type is allowed to be associated to receiver
    virtual bool isDofTypeCompatible(dofType type) const { return false; }
    /**
     * Checks internal data consistency in node.
     * Current implementation checks (when receiver has slave dofs) if receiver has the same
     * coordinate system as master dofManager of slave dof.
     * @return Nonzero if receiver check is o.k.
     */
    virtual int checkConsistency();
    /**
     * Local renumbering support. For some tasks (parallel load balancing, for example) it is necessary to
     * renumber the entities. The various FEM components (such as nodes or elements) typically contain
     * links to other entities in terms of their local numbers, etc. This service allows to update
     * these relations to reflect updated numbering. The renumbering function is passed, which is supposed
     * to return an updated number of specified entity type based on old number.
     */
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);

    /**@name Advanced functions */
    //@{
    /**
     * Sets number of dofs of the receiver; Deallocates existing DOFs;
     * Resizes the dofArray accordingly
     */
    void setNumberOfDofs(int _ndofs);
    /** Sets i-th DOF of receiver to given DOF */
    void setDof(int i, Dof *dof);
    //@}

    /**
     * Adds the given Dof into the receiver. The dofID of scheduled DOF should not be present
     * in receiver as multiple DOFs with same DofID are not allowed. The given DOF is appended
     * at the end of the dofArray.
     * @param dof
     */
    void appendDof(Dof *dof);

    /**
     * Removes Dof with given id from dofArray.
     * @param id
     */
    void removeDof(DofIDItem id);

    /**
     * Checks if receiver contains dof with given ID.
     * @param id Dof ID to check for.
     * @return True if receiver has dof with given id.
     * @see DofIDItem
     */
    bool hasDofID(DofIDItem id);

#ifdef __OOFEG
    virtual void drawYourself(oofegGraphicContext &context) { }
#endif

#if defined( __PARALLEL_MODE ) || defined( __ENABLE_COMPONENT_LABELS )
    /// @return Receivers globally unique number.
    int giveGlobalNumber() const { return globalNumber; }
    int giveLabel() const { return globalNumber; }
    /**
     * Sets receiver global number.
     * @param number New global number for receiver.
     */
    void setGlobalNumber(int number) { globalNumber = number; }
#endif
#ifdef __PARALLEL_MODE
    /**
     * Return dofManagerParallelMode of receiver. Defined for __Parallel_Mode only.
     */
    dofManagerParallelMode giveParallelMode() const { return parallel_mode; }
    /** Sets parallel mode of receiver */
    void setParallelMode(dofManagerParallelMode _mode) { parallel_mode = _mode; }
    /**
     * Packs specific  DOF Manager's dofs unknowns into communication buffer.
     * @param buff Communication buffer to pack data.
     * @param type Physical meaning of  unknown.
     * @param mode Mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @param stepN Time step when unknown requested. See documentation of particular EngngModel
     * class for valid stepN values (most implementation can return only values for current
     * and possibly for previous time step).
     * @return Nonzero if successful
     */
    int packDOFsUnknowns(CommunicationBuffer &buff, EquationID type, ValueModeType mode, TimeStep *stepN);
    /**
     * Returns partition list of receiver.
     * @return Partition array.
     */
    const IntArray *givePartitionList()  { return & partitions; }
    /** Sets receiver's partition list. */
    void setPartitionList(const IntArray *_p) { partitions = * _p; }
    /// Removes given partition from receiver list.
    void removePartitionFromList(int _part) {
        int _pos = partitions.findFirstIndexOf(_part);
        if ( _pos ) { partitions.erase(_pos); } }
    /// Merges receiver partition list with given lists.
    void mergePartitionList(IntArray &_p);
    /**
     * Returns number of partitions sharing given receiver (=number of shared partitions + local one).
     */
    int givePartitionsConnectivitySize();
    /// Returns true if receiver is locally maintained.
    bool isLocal();
    /// Returns true if receiver is shared.
    bool isShared() { if ( parallel_mode == DofManager_shared ) { return true; } else { return false; } }
#endif

protected:
    virtual IRResultType resolveDofIDArray(InputRecord *ir, IntArray &dofIDArry);

    /// Computes transformation matrix between DOFs in nodal c.s. and master DOFs.
    void computeSlaveDofTransformation(FloatMatrix &answer, const IntArray *dofMask);
    IntArray *giveCompleteGlobalDofIDArray() const;
};
} // end namespace oofem
#endif // dofmanager_h

