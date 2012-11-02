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

#ifndef trdirshell_h
#define trdirshell_h

#include "cct.h"
#include "feinterpol2d.h"

namespace oofem {
#ifndef __CHARTENSOR
 #define __CHARTENSOR
enum CharTensor {
    LocalStrainTensor,
    GlobalStrainTensor,
    LocalCurvatureTensor,
    GlobalCurvatureTensor,

    LocalForceTensor,
    GlobalForceTensor,
    LocalMomentumTensor,
    GlobalMomentumTensor
};
#endif

class FEI3dShellTrQuad;

/**
 * This class represent a 7 parameter shell element. 
 * Each node has 63 degrees of freedom (displ. vec., director vec., inhomogeneous thickness strain ).
 *
 * @author Jim Brouzoulis
 * @date 2012-10-xx
 */
class TrDirShell : public NLStructuralElement
{
protected:
    int numberOfGaussPoints;	
	static FEI3dShellTrQuad interpolation;

    static bool __initialized;
	static IntArray ordering_x;
	static IntArray ordering_m;
	static IntArray ordering_gam;

    static bool initOrdering() {
        ordering_x.setValues(18, 1, 2, 3, 8, 9, 10, 15, 16, 17, 22, 23, 24, 29, 30 ,31, 36, 37, 38);
		ordering_m.setValues(18, 4, 5, 6, 11, 12, 13, 18, 19, 20, 25, 26, 27, 32, 33 ,34, 39, 40, 41);
		ordering_gam.setValues(6, 7, 14, 21, 28, 35, 42);
        //edge_ordering [ 0 ].setValues(6,  1, 2, 4, 5, 10, 11);
        //edge_ordering [ 1 ].setValues(6,  4, 5, 7, 8, 12, 13);
        //edge_ordering [ 2 ].setValues(6,  7, 8, 1, 2, 14, 15);
        return true;
    }

	FloatArray initialNodeDirectors[6];
	void setupInitialNodeDirectors();
	FloatArray &giveInitialNodeDirector(int i){return this->initialNodeDirectors[i-1];};

public:
    TrDirShell(int n, Domain *d);	// constructor
    virtual ~TrDirShell() { }		// destructor -> declaring as virtual will make each subclass call their respective destr.
	virtual int computeNumberOfDofs(EquationID ut) { return 42; } // 6*7dofs
	virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

	virtual void TrDirShell :: computeGaussPoints();
	virtual IRResultType initializeFrom(InputRecord *ir);

	virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
	
	// Interpolation functions
	virtual void evalN(FloatArray &answer, const FloatArray &lcoords);
    virtual void evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords);
	void   giveDerivativeXi(FloatArray &n, const FloatArray &lcoords);
    void   giveDerivativeEta(FloatArray &n, const FloatArray &lcoords);


	// Base vectors
	void evalInitialDirectorAt(GaussPoint *gp, FloatArray &answer);
//	void evalInitialNodeDirectors(FloatMatrix &answer);
	

	void evalInitialCovarBaseVectorsAt(GaussPoint *gp, FloatArray &G1, FloatArray &G2, FloatArray &G3);
	void evalInitialContravarBaseVectorsAt(GaussPoint *gp, FloatArray &G1, FloatArray &G2, FloatArray &G3);

	void giveDualBase(const FloatArray &G1, const FloatArray &G2, const FloatArray &G3, FloatArray &g1, FloatArray &g2, FloatArray &g3 );
	void evalCovarBaseVectorsAt(GaussPoint *gp, FloatArray &g1, FloatArray &g2, FloatArray &g3, TimeStep *tStep);
	void evalContravarBaseVectorsAt(GaussPoint *gp, FloatArray &g1, FloatArray &g2, FloatArray &g3, TimeStep *tStep);


		// Stress and strain
	void computeFAt(GaussPoint *gp, FloatMatrix &answer);
	void computeCovarStressAt(GaussPoint *gp, FloatArray &answer);
	void transInitialCartesianToInitialContravar(GaussPoint *gp, const FloatArray &VoightMatrix, FloatArray &answer);


    // definition & identification
    virtual const char *giveClassName() const { return "TrDirShell"; }
    virtual classType giveClassID() const { return TrDirShellClass; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_2; }
	virtual MaterialMode giveMaterialMode() { return _3dMat; }

	virtual integrationDomain  giveIntegrationDomain() { return _Triangle; } // write new wedge-like type

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

	virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS);
	virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);

	virtual void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer,
                                             MatResponseMode rMode, GaussPoint *gp,
                                             TimeStep *tStep);
	
	virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    void computeConvectiveMassForce(FloatArray &answer, TimeStep *tStep);
    void computeThicknessMappingCoeff(FloatArray &answer);
 /**
     * Computes numerically stiffness matrix of receiver. Default implementation computes element stiffness using
     * @f$ K=\int_v B^{\mathrm{T}} D B \mathrm{d}V @f$ formulae, where @f$ B @f$ is element geometric matrix and @f$ D @f$ is material stiffness matrix.
     * No geometrical nonlinearity is taken into account. NUmerical integration procedure uses integrationRulesArray
     * for numerical integration. Support for reduced or selected integration is implemented. The individual integration
     * rules are assumed to correspond to different terms from which the overall matrix is assembled.
     *
     * If numberOfIntegrationRules is equal to
     * 1, the full integration of all coefficients is performed. Otherwise, integration is performed using following rules.
     * Each integration rule can specify start and end strain index of strain vector components for which is valid.
     * It is necessary to ensure that these start and end indexes, dividing geometrical matrix into blocks,
     * are not overlapping and that each strain component is included.
     *
     * Then stiffness matrix is obtained as summation of integrals @f$ I_{ij}=\int_v B^{\mathrm{T}}_i D_{ij} B_j \mathrm{d}V @f$
     * where @f$ B_i @f$ is i-th block of geometrical matrix and @f$ D_{ij} @f$ is corresponding constitutive sub-matrix.
     * The geometrical matrix is obtained using computeBmatrixAt service and the constitutive matrix is obtained using
     * computeConstitutiveMatrixAt service.
     * The @f$ I_{ij} @f$ integral is evaluated using such integration rule, which is valid for i-th or j-th block
     * and has smaller number of integration points.
     *
     * For higher numerical performance, only one half of stiffness matrix is computed and answer is then symmetrized.
     * Therefore, if element matrix will be generally nonsymmetric, one must specialize this method.
     * Finally, the result is transformed into global coordinate system (or nodal coordinate system, if it is defined).
     *
     * @param answer Computed stiffness matrix (symmetric).
     * @param rMode Response mode.
     * @param tStep Time step.
     */
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

	// stress equivalent vector in nodes (vector of internal forces)
    // - mainly for nonLinear Analysis.
    /**
     * Returns equivalent nodal forces vectors. Useful for nonlinear analysis.
     * Default implementation computes result as @f$ F=\int_v B^{\mathrm{T}} \sigma \mathrm{d}V @f$, where @f$ \sigma @f$ is the
     * real element stress vector obtained using computeStressVector service (if useUpdatedGpRecord=0) or
     * (if useUpdatedGpRecord=1) from integration point status.
     * The geometric matrix is obtained using computeBmatrixAt service.
     * Integration is performed using default integration rule, which should produce always valid results,
     * assuming that strains used for computation of stresses are valid.
     * @param answer Internal nodal forces vector.
     * @param tStep Time step.
     * @param useUpdatedGpRecord If equal to zero, the stresses in integration points are computed (slow but safe), else if
     * nonzero the stresses are taken directly from integration point status (should be derived from StructuralMaterialStatus)
     * (fast, but engineering model must ensure valid status data in each integration point).
     */
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    void computeSectionalForces(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);

    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode);

};



} // end namespace oofem
#endif // cct3d_h
