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

#include "tr_dirshell.h"
#include "node.h"
#include "load.h"
#include "structuralms.h"
#include "mathfem.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspnt.h"
//#include "fei3dshelltrquad.h"

namespace oofem {
	//FEI3dShellTrQuad TrDirShell :: interpolation;

	IntArray TrDirShell :: ordering_x(18);
	IntArray TrDirShell :: ordering_m(18);
	IntArray TrDirShell :: ordering_gam(6);
	
	bool TrDirShell :: __initialized = TrDirShell :: initOrdering();


TrDirShell :: TrDirShell(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
{
	this->numberOfIntegrationRules = 1;
	this->numberOfDofMans = 6;
    this->numberOfGaussPoints = 1;
    this->computeGaussPoints();
	//FEI3dShellTrQuad TrDirShell :: interpolation(1, 2);
	//IntArray TrDirShell :: initOrdering();
	
	//this->setupInitialNodeDirectors(this->initialNodeDirectors);
	
}

IRResultType TrDirShell :: initializeFrom(InputRecord *ir)
{
    this->NLStructuralElement :: initializeFrom(ir);
	this->setupInitialNodeDirectors();
    return IRRT_OK;
}


int 
TrDirShell :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
	
	//FloatMatrix G;
	//interpolation.local2global(answer, lcoords, FEIElementGeometryWrapper(this));


	//this->evalContravarBaseVectors(G, lcoords, cellgeo, zeta)
	

    return 1;
    
}

void
TrDirShell :: evalN(FloatArray &answer, const FloatArray &lcoords)
{// copied from fei2dtrquad
    double l1 = lcoords.at(1);
    double l2 = lcoords.at(2);
    double l3 = 1. - l1 - l2;

    answer.resize(6);

    answer.at(1) = ( 2. * l1 - 1. ) * l1;
    answer.at(2) = ( 2. * l2 - 1. ) * l2;
    answer.at(3) = ( 2. * l3 - 1. ) * l3;
    answer.at(4) = 4. * l1 * l2;
    answer.at(5) = 4. * l2 * l3;
    answer.at(6) = 4. * l3 * l1;

    return;
}

void
TrDirShell :: evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords)
{	// Returns matrix with derivatives wrt local coordinates, mostly copied from fei2dtrquad
    answer.resize(6, 2);
    int i;
    FloatArray dndxi(6), dndeta(6);

    this->giveDerivativeXi(dndxi, lcoords);
    this->giveDerivativeEta(dndeta, lcoords);

    for ( i = 1; i <= 6; i++ ) {
        answer.at(i, 1) = dndxi.at(i); 
        answer.at(i, 2) = dndeta.at(i);
    }
	
}

void
TrDirShell :: giveDerivativeXi(FloatArray &n, const FloatArray &lc)
{
    double l1, l2, l3;

    l1 = lc.at(1);
    l2 = lc.at(2);
    l3 = 1.0 - l1 - l2;

    n.resize(6);

    n.at(1) = 4.0 * l1 - 1.0;
    n.at(2) = 0.0;
    n.at(3) = -1.0 * ( 4.0 * l3 - 1.0 );
    n.at(4) = 4.0 * l2;
    n.at(5) = -4.0 * l2;
    n.at(6) = 4.0 * l3 - 4.0 * l1;
}

void
TrDirShell :: giveDerivativeEta(FloatArray &n, const FloatArray &lc)
{
    double l1, l2, l3;

    l1 = lc.at(1);
    l2 = lc.at(2);
    l3 = 1.0 - l1 - l2;

    n.resize(6);

    n.at(1) = 0.0;
    n.at(2) = 4.0 * l2 - 1.0;
    n.at(3) = -1.0 * ( 4.0 * l3 - 1.0 );
    n.at(4) = 4.0 * l1;
    n.at(5) = 4.0 * l3 - 4.0 * l2;
    n.at(6) = -4.0 * l1;
}




void 
TrDirShell :: evalInitialCovarBaseVectorsAt(GaussPoint *gp, FloatArray &G1, FloatArray &G2, FloatArray &G3)
{
	int i;
	double x, y, z, Mx, My, Mz, zeta;
	//FloatArray lcoords = *gp->giveCoordinates();
    FloatArray lcoords = *gp->giveCoordinates();
	zeta = lcoords.at(3);
	FloatArray N, M;
	FloatMatrix dNdxi, Mmat;
	
	// In plane base vectors
	this->evaldNdxi(dNdxi, lcoords);
	//this->evalInitialNodeDirectors(Mmat);
	G1.resize(3); G2.resize(3);
	G1.zero(); G2.zero();
	for ( i = 1; i <= 6; i++ ) {
		FloatArray *nodeI = this->giveNode(i)->giveCoordinates();
		x = nodeI->at(1); y = nodeI->at(2); z = nodeI->at(3);

		M=this->giveInitialNodeDirector(i);
		Mx = M.at(1); My = M.at(2); Mz = M.at(3);

        G1.at(1) += dNdxi.at(i,1)*( x + zeta*Mx ); 
		G1.at(2) += dNdxi.at(i,1)*( y + zeta*My ); 
		G1.at(3) += dNdxi.at(i,1)*( z + zeta*Mz ); 
        G2.at(1) += dNdxi.at(i,2)*( x + zeta*Mx ); 
		G2.at(2) += dNdxi.at(i,2)*( y + zeta*My ); 
		G2.at(3) += dNdxi.at(i,2)*( z + zeta*Mz ); 
    }
	// Out of plane base vector = director
	this->evalInitialDirectorAt(gp, G3); // G3=M
	
}

void
TrDirShell :: evalInitialContravarBaseVectorsAt(GaussPoint * gp, FloatArray &Gcon1, FloatArray &Gcon2, FloatArray &Gcon3)
{	
	FloatArray Gcov1, Gcov2, Gcov3;
	this->evalInitialCovarBaseVectorsAt(gp, Gcov1, Gcov2, Gcov3);
	this->giveDualBase(Gcov1, Gcov2, Gcov3, Gcon1, Gcon2, Gcon3);
}


void
TrDirShell :: evalContravarBaseVectorsAt(GaussPoint *gp, FloatArray &gcon1, FloatArray &gcon2, FloatArray &gcon3, TimeStep *tStep)
{	
	FloatArray gcov1, gcov2, gcov3;
	this->evalCovarBaseVectorsAt(gp, gcov1, gcov2, gcov3, tStep);
	this->giveDualBase(gcov1, gcov2, gcov3, gcon1, gcon2, gcon3);
}

void
TrDirShell :: giveDualBase(const FloatArray &G1, const FloatArray &G2, const FloatArray &G3, FloatArray &g1, FloatArray &g2, FloatArray &g3 )
{	
	FloatMatrix gmat, ginv, test;
	
	gmat.resize(3,3);
	gmat.at(1,1) = G1.dotProduct(G1); gmat.at(1,2) = G1.dotProduct(G2); gmat.at(1,3) = G1.dotProduct(G3);
	gmat.at(2,2) = G2.dotProduct(G2); gmat.at(2,3) = G2.dotProduct(G3); gmat.at(3,3) = G3.dotProduct(G3);
	gmat.symmetrized();
	
	ginv.beInverseOf(gmat);

	g1.resize(3); g1.zero(); g2.resize(3); g2.zero(); g3.resize(3); g3.zero();

	g1.add(ginv.at(1,1),G1); g1.add(ginv.at(1,2),G2); g1.add(ginv.at(1,3),G3); 
	g2.add(ginv.at(2,1),G1); g2.add(ginv.at(2,2),G2); g2.add(ginv.at(2,3),G3);
	g3.add(ginv.at(3,1),G1); g3.add(ginv.at(3,2),G2); g3.add(ginv.at(3,3),G3);

	// Add test of orthogonality, should be a diagonal matrix. 
	/*test.resize(3,3);
	test.at(1,1) = G1.dotProduct(g1); test.at(1,2) = G1.dotProduct(g2); test.at(1,3) = G1.dotProduct(g3);
	test.at(2,2) = G2.dotProduct(g2); test.at(2,3) = G2.dotProduct(g3); test.at(3,3) = G3.dotProduct(g3);
	test.symmetrized();
	test.printYourself(); */
}

void 
TrDirShell :: evalInitialDirectorAt(GaussPoint *gp, FloatArray &answer)
{	// Interpolates between the node directors
	FloatArray lcoords = *gp->giveCoordinates();
	FloatArray N;
	this->evalN(N, lcoords);
	
	answer.resize(3); answer.zero();
	for (int i = 1; i <= 6; i++ ) {
		answer.add( N.at(i), this->giveInitialNodeDirector(i) ); 
	}

}


void 
TrDirShell :: setupInitialNodeDirectors()
{	// If the directors are not present in the input file, then they should be approximated as the normal to the initial surface.
	FloatMatrix dNdxi;
	FloatArray M, G1, G2, lcoords, nodeLocalXiCoords, nodeLocalEtaCoords;

	// Compute directors as normals to the surface

	// Set the local coordinates for the element nodes
	nodeLocalXiCoords.setValues( 6, 1., 0., 0., .5, 0., .5); // corner nodes then midnodes, uncertain of node numbering
	nodeLocalEtaCoords.setValues(6, 0., 1., 0., .5, .5, 0.);
	lcoords.resize(2); G1.resize(3); G2.resize(3); M.resize(3);		


	for (int node = 1; node <= 6; node++ ) {

		this->initialNodeDirectors[node-1].resize(3);
		this->initialNodeDirectors[node-1].zero();

		lcoords.at(1) = nodeLocalXiCoords.at(node); 
		lcoords.at(2) = nodeLocalEtaCoords.at(node);
		//lcoords.printYourself();

		this->evaldNdxi(dNdxi, lcoords);
		
		G1.zero(); G2.zero(); M.zero();
		for ( int i = 1; i <= 6; i++){	// base vectors of the initial surface
			FloatArray *nodeI = this->giveNode(i)->giveCoordinates();
			G1.add( dNdxi.at(i,1), *nodeI );
			G2.add( dNdxi.at(i,2), *nodeI );
		}
		M.beVectorProductOf(G1,G2);
		M.normalize();

		this->initialNodeDirectors[node-1].add(M); 	
		//this->initialNodeDirectors[node-1].printYourself(); 	
	}
	
}

void
TrDirShell :: evalCovarBaseVectorsAt(GaussPoint *gp, FloatArray &g1, FloatArray &g2, FloatArray &g3, TimeStep *tStep)
{
	double zeta, gam;
	FloatArray lcoords = *gp->giveCoordinates();
	zeta = lcoords.at(3);
	FloatArray a;
	FloatMatrix B;

	this->computeBmatrixAt(gp, B, 1, ALL_STRAINS);
	//this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, a);
	this->giveUpdatedSolutionVector(a, tStep);
	FloatArray eps;	      // generalized strain
	eps.beProductOf(B,a); // [dxdxi, dmdxi, m, dgamdxi, gam]^T

	
	FloatArray dxdxi, m, dmdxi, dgamdxi, test;
	dxdxi.setValues(6, eps.at(1), eps.at(2), eps.at(3), eps.at(4), eps.at(5), eps.at(6) );
	dxdxi.printYourself();
	m.setValues(3, eps.at(13), eps.at(14), eps.at(15) );
	m.printYourself();
	dmdxi.setValues(6, eps.at(7), eps.at(8), eps.at(9), eps.at(10), eps.at(11), eps.at(12) );
	dmdxi.printYourself();
	gam = eps.at(18);
	dgamdxi.setValues(2, eps.at(16), eps.at(17) );
	dgamdxi.printYourself();

	g1.resize(3), g2.resize(3), g3.resize(3);
	double fac1 = ( zeta + 0.5*gam*zeta*zeta );
	double fac2 = ( 0.5*zeta*zeta );
	double fac3 = ( 1.0 + zeta*gam );
	
	// first compute as initial base vectors and then update, gi = Gi + u,i
	//this->evalInitialCovarBaseVectorsAt(gp, g1, g2, g3);
	g1.at(1) = dxdxi.at(1) + fac1*dmdxi.at(1) + fac2*m.at(1)*dgamdxi.at(1);
	g1.at(2) = dxdxi.at(2) + fac1*dmdxi.at(2) + fac2*m.at(2)*dgamdxi.at(1);
	g1.at(3) = dxdxi.at(3) + fac1*dmdxi.at(3) + fac2*m.at(3)*dgamdxi.at(1);

	g2.at(1) = dxdxi.at(4) + fac1*dmdxi.at(4) + fac2*m.at(1)*dgamdxi.at(2);
	g2.at(2) = dxdxi.at(5) + fac1*dmdxi.at(5) + fac2*m.at(2)*dgamdxi.at(2);
	g2.at(3) = dxdxi.at(6) + fac1*dmdxi.at(6) + fac2*m.at(3)*dgamdxi.at(2);

	g3.at(1) = fac3*m.at(1); 
	g3.at(2) = fac3*m.at(2); 
	g3.at(3) = fac3*m.at(3);


	// Extract x, m, gam from solution vector
	/*
	this->ordering_x.printYourself();
	this->ordering_m.printYourself();
	this->ordering_gam.printYourself();
	a_x.resize(18); a_m.resize(18);
	for ( i = 1; i <= 18; i++ ) {
		a_x.at(i) = a.at(ordering_x.at(i));
		a_m.at(i) = a.at(ordering_m.at(i));
	}
	a_gam.resize(6);
	for ( i = 1; i <= 6; i++ ) {
		a_gam.at(i) = a.at(ordering_gam.at(i));
	}
	*/



}


void
TrDirShell :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    answer.setValues(7, D_u, D_v, D_w, w_u, w_v, w_w, gam);
}


void
TrDirShell :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body loads, at stepN.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    
	/*
	double dens, dV, load;
    GaussPoint *gp = NULL;
    FloatArray force;
    FloatMatrix T;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        _error("computeBodyLoadVectorAt: unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, stepN, mode);

    if ( force.giveSize() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

        dens = this->giveMaterial()->give('d', gp);
        dV   = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness);

        answer.resize(18);
        answer.zero();

        load = force.at(1) * dens * dV / 3.0;
        answer.at(1) = load;
        answer.at(7) = load;
        answer.at(13) = load;

        load = force.at(2) * dens * dV / 3.0;
        answer.at(2) = load;
        answer.at(8) = load;
        answer.at(14) = load;

        load = force.at(3) * dens * dV / 3.0;
        answer.at(3) = load;
        answer.at(9) = load;
        answer.at(15) = load;

        // transform result from global cs to local element cs.
        if ( this->computeGtoLRotationMatrix(T) ) {
            answer.rotatedWith(T, 'n');
        }
    } else {
        answer.resize(0);          // nil resultant
    }
	*/
	answer.resize(0);
}


void
TrDirShell :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    /*
	int i, j;
    GaussPoint *gp;
    FloatArray v;

#if defined ( __PARALLEL_MODE ) || defined ( __ENABLE_COMPONENT_LABELS )
    fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );
#else
    fprintf(file, "element %d :\n", number);
#endif

	
    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        for ( j = 0; j < integrationRulesArray [ i ]->getNumberOfIntegrationPoints(); j++ ) {
            gp = integrationRulesArray [ i ]->getIntegrationPoint(j);

            // gp   -> printOutputAt(file,stepN) ;
			fprintf( file, "  GP %2d.%-2d :", i + 1, gp->giveNumber() );

            this->giveIPValue(v, gp, IST_ShellStrainCurvatureTensor, tStep);
            fprintf(file, "  strains ");
            fprintf( file,
                    " % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e ",
                    v.at(1), v.at(2), v.at(3),  2. * v.at(4), 2. * v.at(5), 2. * v.at(6),
                    v.at(7), v.at(8), v.at(9),  2. * v.at(10), 2. * v.at(11), 2. * v.at(12) );

            this->giveIPValue(v, gp, IST_ShellForceMomentumTensor, tStep);
            fprintf(file, "\n              stresses");
            fprintf( file,
                    " % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e ",
                    v.at(1), v.at(2), v.at(3),  v.at(4), v.at(5), v.at(6),
                    v.at(7), v.at(8), v.at(9),  v.at(10), v.at(11), v.at(12) );

            fprintf(file, "\n");
        }
    }
	*/
}

void
TrDirShell :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
	answer.resize(7, 42);
    answer.zero();

    FloatArray lcoords = *gp->giveCoordinates();
	FloatArray N;
	FloatMatrix dNdxi;
	this->evalN(N, lcoords);


	/*    18   18    6
  	   3 [N_x   0    0
	   3   0   N_m   0
	   1   0    0  N_gmm ]

	*/	
	int i, j;
	for( i = 1, j = 0; i<=6; i++, j+=3  ){ 
		answer.at(1,1+j) = N.at(i);
		answer.at(2,2+j) = N.at(i);
		answer.at(3,3+j) = N.at(i);
		answer.at(4,18+1+j) = N.at(i);
		answer.at(5,18+2+j) = N.at(i);
		answer.at(6,18+3+j) = N.at(i);
	}

	for( i = 1; i<=6; i++){ 
	    answer.at(7,36+i) = N.at(i);
	}

	
}


void
TrDirShell :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li , int ui)
/* Returns the  matrix {B} of the receiver, evaluated at aGaussPoint. Such that
   B*a = [dxbar_dxi, dwdxi, w, dgamdxi, gam]^T, where a is the vector of unknowns
*/
{
	FloatArray lcoords = *gp->giveCoordinates();
	answer.resize(18, 42);
    answer.zero();
	FloatArray N;
	FloatMatrix dNdxi;
	this->evalN(N, lcoords);
	this->evaldNdxi(dNdxi, lcoords);

	/*    18   18   6
  	   6 [B_u   0   0
	   6   0   B_w  0
	   3   0   N_w  0
	   2   0    0  B_gam 
	   1   0    0  N_gam] 
	*/
	
	int i, j, pos;
	// First column
	for( i = 1, j = 0; i<=6; i++, j+=3  ){ 
		answer.at(1,1+j) = dNdxi.at(i,1);
		answer.at(2,2+j) = dNdxi.at(i,1);
		answer.at(3,3+j) = dNdxi.at(i,1);
		answer.at(4,1+j) = dNdxi.at(i,2);
		answer.at(5,2+j) = dNdxi.at(i,2);
		answer.at(6,3+j) = dNdxi.at(i,2);
	}

	// Second column
	pos =18;
	for( i = 1, j = 0; i<=6; i++, j+=3  ){ 
		answer.at(6+1,pos+1+j) = dNdxi.at(i,1);
		answer.at(6+2,pos+2+j) = dNdxi.at(i,1);
		answer.at(6+3,pos+3+j) = dNdxi.at(i,1);
		answer.at(6+4,pos+1+j) = dNdxi.at(i,2);
		answer.at(6+5,pos+2+j) = dNdxi.at(i,2);
		answer.at(6+6,pos+3+j) = dNdxi.at(i,2);
		answer.at(6+7,pos+1+j) = N.at(i);
		answer.at(6+8,pos+2+j) = N.at(i);
		answer.at(6+9,pos+3+j) = N.at(i);
	}

	// Third column
	pos =36;
	for( i = 1, j = 0; i<=6; i++, j+=1  ){ 
		answer.at(15+1,pos+1+j) = dNdxi.at(i,1);
		answer.at(15+2,pos+1+j) = dNdxi.at(i,2);
		answer.at(15+3,pos+1+j) = N.at(i);
	}


}



void
TrDirShell :: computeConstitutiveMatrixAt(FloatMatrix &answer,
                                                 MatResponseMode rMode, GaussPoint *gp,
                                                 TimeStep *tStep)
// Returns the  material matrix {E} of the receiver.
// type of matrix is determined by this->giveMaterialMode()
// rMode parameter determines type of stiffness matrix to be requested
// (tangent, secant, ...)
{
	answer.resize(18, 18);
    answer.zero();
    //( ( StructuralCrossSection * ) this->giveCrossSection() )
    //->giveCharMaterialStiffnessMatrix(answer, rMode, gp, tStep);
}


void 
TrDirShell :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
		integrationRulesArray = new IntegrationRule * [ 1 ];
		integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Wedge, this->numberOfGaussPoints, _3dMat);
    }

}




void
TrDirShell :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                            TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    // Hi-jack routine. Do whatever I want  :)

	IntegrationRule *iRule = integrationRulesArray [ 0 ];
	GaussPoint *gp;
	gp = iRule->getIntegrationPoint(0);

	FloatMatrix Mmat;
	FloatArray lcoords, m, G1, G2, G3, g1, g2, g3;
	
	lcoords.resize(3);
	lcoords.setValues(3, .3333, .3333, 0.); // xi1, xi2, zeta



	answer.resize(42,42);
	answer.beUnitMatrix();
}


void 
TrDirShell :: computeFAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *stepN){
	
	FloatArray gcov1, gcov2, gcov3, Gcon1, Gcon2, Gcon3;
	this->evalCovarBaseVectorsAt(gp, gcov1, gcov2, gcov3, stepN);
	this->evalInitialContravarBaseVectorsAt(gp, Gcon1, Gcon2, Gcon3);
	gcov1.printYourself();
	gcov2.printYourself();
	gcov3.printYourself();
	Gcon1.printYourself();
	Gcon2.printYourself();
	Gcon3.printYourself();
	FloatMatrix F, F1, F2, F3;
	answer.resize(3,3); answer.zero();
	F1.beDyadicProductOf(gcov1,Gcon1);
	F2.beDyadicProductOf(gcov2,Gcon2);
	F3.beDyadicProductOf(gcov3,Gcon3);

	answer.add(F1);	answer.add(F2);	answer.add(F3);
}


void
TrDirShell :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step stepN. The nature of these strains depends
// on the element's type.
{

	// Deformation gradient
	FloatMatrix F;
	this->computeFAt(gp, F, stepN);
	F.printYourself();
	// Green-Lagrange strain tensor
	FloatMatrix temp, E;
	temp.beTProductOf(F,F);
	temp.at(1,1) += -1;
	temp.at(2,2) += -1;
	temp.at(3,3) += -1;
	E.add(0.5, temp);
	
	//E.printYourself();

	//answer.resize(6);
    //answer.zero();
	// Voight order: 11, 22, 33, 23, 13, 12
	answer.beReducedVectorForm(E);
	answer.printYourself();
	//answer.setValues(6, E.at(1,1), E.at(2,2), E.at(3,3), E.at(2,3), E.at(1,3), E.at(1,2) ); 

}


void
TrDirShell :: transInitialCartesianToInitialContravar(GaussPoint *gp, const FloatArray &VoightMatrix, FloatArray &answer){

    FloatArray Gcov1, Gcov2, Gcov3;
    this->evalInitialContravarBaseVectorsAt(gp, Gcov1, Gcov2, Gcov3);

	FloatArray E1(3), E2(3), E3(3);
    E1.zero(); E2.zero(); E3.zero();
    E1.at(1) = E2.at(2) = E3.at(1) = 1.;

	// Transformation matrix
    // Todo: Simplify! E1-E3 are cartesian vectors so the dot products will only be components of G1-G3
    FloatMatrix EG(3,3), EGt(3,3);
    EG.at(1,1) = E1.dotProduct(Gcov1); EG.at(1,2) = E1.dotProduct(Gcov2); EG.at(1,3) = E1.dotProduct(Gcov3);
    EG.at(2,2) = E2.dotProduct(Gcov2); EG.at(2,3) = E2.dotProduct(Gcov3); EG.at(3,3) = E3.dotProduct(Gcov3);
    EG.symmetrized();
	EGt.beTranspositionOf(EG);
	
	// transform according to: EG^T*mat*EG
    FloatMatrix temp(3,3), mat(3,3);
	mat.beMatrixForm(VoightMatrix);
    temp.beProductTOf(EG,mat);
    temp.beProductOf(temp,EG);
	answer.beReducedVectorForm(temp);

}


void
TrDirShell :: giveUpdatedSolutionVector(FloatArray &answer, TimeStep *tStep){
	
	FloatArray *Xi, Mi,temp;
	this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, temp);
	temp.resize(42);
	for( int i = 1, j = 0; i <= 6; i++, j += 7 ){
		Xi = this->giveNode(i)->giveCoordinates();
		Mi = this->giveInitialNodeDirector(i);
		temp.at(1+j) += Xi->at(1);
		temp.at(2+j) += Xi->at(2);
		temp.at(3+j) += Xi->at(3);
		temp.at(4+j) += Mi.at(1);
		temp.at(5+j) += Mi.at(2);
		temp.at(6+j) += Mi.at(3);
		//answer.at(7+j) += 0, gam(t=0)=0 so no update necessary
	}
	answer.resize(42);
	for( int i = 1; i <= 18; i++ ){
		answer.at(i)    = temp.at( ordering_x.at(i) );
		answer.at(i+18) = temp.at( ordering_m.at(i) );
	}
	for( int i = 1; i <= 6; i++ ){
		answer.at(i+36)    = temp.at( ordering_gam.at(i) );
	}

	// group same fields together
}


void
TrDirShell :: giveInternalForcesVector(FloatArray &answer,
                                              TimeStep *tStep, int useUpdatedGpRecord)
//
// returns nodal representation of real internal forces - necessary only for
// non-linear analysis.
// if useGpRecord == 1 then data stored in gp->giveStressVector() are used
// instead computing stressVector through this->ComputeStressVector();
// this must be done after you want internal forces after element->updateYourself()
// has been called for the same time step.
//
{
    GaussPoint *gp;
    Material *mat = this->giveMaterial();

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    FloatArray lcoords;
    double zeta;
    FloatMatrix B, Bt;
    FloatArray cartStressVector, contravarStressVector, sectionalForces, BF, a;
    double dV;

    // do not resize answer to computeNumberOfDofs(EID_MomentumBalance)
    // as this is valid only if receiver has no nodes with slaves
    // zero answer will resize accordingly when adding first contribution
    answer.resize(42);
	answer.zero();

    
	for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        lcoords = *gp->giveCoordinates();
        zeta = lcoords.at(3); 


        

        if ( useUpdatedGpRecord == 1 ) {
            cartStressVector = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )
                                ->giveStressVector();
        } else {
            this->computeStressVector(cartStressVector, gp, tStep);
        }


        //
        // updates gp stress and strain record  acording to current
        // increment of displacement
        //
        if ( cartStressVector.giveSize() == 0 ) {
            break;
        }
        
        FloatArray f(18), g1, g2, g3, S1g(3), S2g(3), S3g(3), m(3), dm1(3), dm2(3);
        this->transInitialCartesianToInitialContravar(gp, cartStressVector, contravarStressVector);
        this->evalCovarBaseVectorsAt(gp, g1, g2, g3, tStep);
        FloatMatrix S; 
		S.beMatrixForm(cartStressVector);
        double fac1, fac2, fac3, gam, dg1, dg2;

		this->computeBmatrixAt(gp, B);
		this->giveUpdatedSolutionVector(a,tStep);
		FloatArray eps;	      // generalized strain
		eps.beProductOf(B,a); // [dxdxi, dmdxi, m, dgamdxi, gam]^T

		FloatArray dxdxi, dmdxi, dgamdxi;
		dm1.setValues( 3, eps.at(7), eps.at(8), eps.at(9) );
		dm2.setValues( 3, eps.at(10), eps.at(11), eps.at(12) );
		m.setValues(3, eps.at(13), eps.at(14), eps.at(15) );
		dg1 = eps.at(16);
		dg2 = eps.at(17);
		gam = eps.at(18);


        fac1 = zeta + 0.5*gam*zeta*zeta;
        fac2 = 1. + gam*zeta;
        fac3 = 0.5*zeta*zeta;

		// S1g =S(1,i)*g_j, etc
        for( int j=1; j<=3; j++){
            S1g.at(j) = S.at(1,1)*g1.at(j) + S.at(1,2)*g2.at(j) + S.at(1,3)*g3.at(j);
            S2g.at(j) = S.at(2,1)*g1.at(j) + S.at(2,2)*g2.at(j) + S.at(2,3)*g3.at(j);
            S3g.at(j) = S.at(3,1)*g1.at(j) + S.at(3,2)*g2.at(j) + S.at(3,3)*g3.at(j);
        }

        /* The following expressions are sectional forces if integrated over the thickness. However, 
           now they are integrated over the volume.*/

        // Normal forces - N1, N2 
        f.at(1) = S1g.at(1); f.at(2) = S1g.at(2); f.at(3) = S1g.at(3);
        f.at(4) = S2g.at(1); f.at(5) = S2g.at(2); f.at(6) = S2g.at(3);

        // Moments - M1, M2
        f.at(7)  = fac1 * S1g.at(1);
        f.at(8)  = fac1 * S1g.at(2);
        f.at(9)  = fac1 * S1g.at(3);
        f.at(10) = fac1 * S2g.at(1);
        f.at(11) = fac1 * S2g.at(2);
        f.at(12) = fac1 * S2g.at(3);

        // Shear force - T
        f.at(13) = fac2 * S3g.at(1) + fac3 * ( S1g.at(1)*dg1 + S2g.at(1)*dg2) ;
        f.at(14) = fac2 * S3g.at(2) + fac3 * ( S1g.at(2)*dg1 + S2g.at(2)*dg2) ;
        f.at(15) = fac2 * S3g.at(3) + fac3 * ( S1g.at(3)*dg1 + S2g.at(3)*dg2) ;

        // Ms1, Ms2
        f.at(16) = fac3 * m.dotProduct(S1g);
        f.at(17) = fac3 * m.dotProduct(S1g);

        // Ts
        f.at(18) = fac3 * ( dm1.dotProduct(S1g) + dm2.dotProduct(S2g) + zeta * m.dotProduct(S3g)  ) ;
            
        // compute nodal representation of internal forces using f = B^T*f*dV
        Bt.beTranspositionOf(B);
        dV = this->computeVolumeAround(gp);
        BF.beProductOf(Bt, f);
        BF.times(dV);
        answer.add(BF);
    }

    // if inactive update state, but no contribution to global system
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
	

}

void 
TrDirShell :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep){
	FloatMatrix mConsistent;
	this->computeMassMatrix(mConsistent, tStep);
	// Reduce to lumped form 
	// Todo: add algorithm for this
	//answer.resize(mConsistent);
	answer.resize(42,42); answer.zero();
	for( int i = 1; i<=42; i++){
		answer.at(i,i) = mConsistent.at(i,i);
	}

}


void 
TrDirShell :: computeMassMatrix(FloatMatrix &answer, TimeStep *tStep){
	// Analytically integrated over the thickness. Constant density assumed.
	IntegrationRule *iRule = integrationRulesArray [ 0 ];
	GaussPoint *gp;

	//------------------------------
    FloatMatrix N, Nt, Ntm, NtmN, mass;
    FloatArray a, unknowns, m(3);
    double gam, dV;
    
	answer.resize(42,42); answer.zero();
	for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);

		this->giveUpdatedSolutionVector(a, tStep);
		this->computeNmatrixAt(gp, N);
		unknowns.beProductOf(N,a); // [x, m, gam]^T
		m.setValues(3, unknowns.at(4), unknowns.at(5), unknowns.at(6) );
		gam = unknowns.at(7);
		//------------------------------

		// Consistent mass matrix
		/*     3    3    1
  		   3 [a*I  b*I   c*m      [A  B  C
		   3       d*I   e*m    =     D  E
		   q  sym       f*m.m]     sym   F]
		*/
		double a1, a2, a3;
		FloatArray coeff;
		this->computeThicknessMappingCoeff(coeff);
		a1 = coeff.at(1); a2 = coeff.at(2); a3 = coeff.at(3);


		double h, h2, h3, h5, fac1, fac2, fac3, fac4, fac5, fac6, gam2, rho; 
		rho = this->giveMaterial()->give('d', gp);
		h = this->giveCrossSection()->give(CS_Thickness);
		h2 = h*h; h3 = h2*h; h5 = h2*h3;
		gam2 = gam*gam;

		mass.resize(7,7);
		fac1 = a3*h + (a1*h3)/12.;
		fac2 = h3*(40.*a2 + 20.*a3*gam + 3.*a1*h2*gam)/480.;
		fac3 = h3*(20.*a3 + 3.*a1*h2)/480.;
		fac4 = (28.*a3*h3*(80. + 3.*h2*gam2) + 3.*h5*(112.*a2*gam + a1*(112. + 5.*h2*gam2)))/26880.;
		fac5 = h*(56.*a2 + 28.*a3*gam + 5.*a1*h2*gam)/8960.;
		fac6 = h5*(28.*a3 + 5.*a1*h2)/8960.;
		mass.at(1,1) = answer.at(2,2) = answer.at(3,3) = fac1; // A
		mass.at(1,4) = answer.at(2,5) = answer.at(3,6) = fac2; // B
		mass.at(1,7) = fac3 * m.at(1);  answer.at(2,7) = fac3 * m.at(2); answer.at(3,7) = fac3 * m.at(3); // C
		mass.at(4,4) = answer.at(5,5) = answer.at(6,6) = fac4; // D
		mass.at(4,7) = fac5 * m.at(1);  answer.at(5,7) = fac5 * m.at(2); answer.at(6,7) = fac5 * m.at(3); // E
		mass.at(7,7) = fac6 * m.dotProduct(m); // F

        
        Nt.beTranspositionOf(N);
        dV = this->computeVolumeAround(gp)/h ; // Should be an area!
        Ntm.beProductOf(Nt, mass);
		NtmN.beProductOf(Ntm, N);
        NtmN.times(dV);
        answer.add(NtmN);
    
	}
	

/* Salars Matlab code
function [mass,cmass]= prbf(rho,h,gam,gamd,m,md,PHI1,PHI2,M,M1,M2,bc)

%     consistent mass matrix

sc=norm(cross(PHI1,PHI2));
a1=dot(M,cross(M1,M2))/sc;
a2=dot(M,cross(PHI1,M2) + cross(M1,PHI2))/sc;
a3=dot(M,cross(PHI1,PHI2))/sc;
eye3=eye(3);

% mass matrix

macc=[...
    (a3*h + (a1*h^3)/12.)*eye3,(h^3*(40*a2 + 20*a3*gam + 3*a1*h^2*gam))/480.*eye3,(h^3*(20*a3 + 3*a1*h^2)*m)/480.;...
    (h^3*(40*a2 + 20*a3*gam + 3*a1*h^2*gam))/480.*eye3,(28*a3*h^3*(80 + 3*h^2*gam^2) + 3*h^5*(112*a2*gam + a1*(112 + 5*h^2*gam^2)))/26880.*eye3,...
    (h^5*m*(56*a2 + 28*a3*gam + 5*a1*h^2*gam))/8960.;...
    (h^3*(20*a3 + 3*a1*h^2)*m')/480.,(h^5*m'*(56*a2 + 28*a3*gam + 5*a1*h^2*gam))/8960.,(h^5*(28*a3 + 5*a1*h^2)*m'*m)/8960.];

mass=rho*bc'*macc*bc;

*/
}


void 
TrDirShell :: computeConvectiveMassForce(FloatArray &answer, TimeStep *tStep){
	// Analytically integrated over the thickness. Constant density assumed.
	IntegrationRule *iRule = integrationRulesArray [ 0 ];
	GaussPoint *gp;
	gp = iRule->getIntegrationPoint(0);

	//------------------------------
    FloatMatrix N;
    FloatArray a, unknowns, m;
    double gam;
    this->computeNmatrixAt(gp, N);
	this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, a);
   
    unknowns.beProductOf(N,a); // [x, m, gam]^T
	m.setValues(3, unknowns.at(4), unknowns.at(5), unknowns.at(6) );
	gam = unknowns.at(7);


    double a1, a2, a3, h, h2, h3, h5, fac1, fac2, fac3, gam2, rho; 
    FloatArray coeff;
    this->computeThicknessMappingCoeff(coeff);
    a1 = coeff.at(1); a2 = coeff.at(2); a3 = coeff.at(3);

    
    rho = this->giveMaterial()->give('d', gp);
    h = this->giveCrossSection()->give(CS_Thickness);
    h2 = h*h; h3 = h2*h; h5 = h2*h3;
    gam2 = gam*gam;

    // Convective mass "matrix"
	/*     3
  	   3 [ a*m
	   3   b*m
	   1  c*m.dm]*dgam
	*/
    // \dot{m} and \dot{gamma} should be obtained from solution
    FloatArray md; md.resize(3); m.zero();
    double gamd=0;
    answer.resize(7);
    fac1 = rho*h3*(20.*a3 + 3.*a1*h2)/240.;
    fac2 = rho*h5*(56.*a2 + 28.*a3*gam + 5.*a1*h2*gam)/4480.;
    fac3 = rho*h5*(28.*a3 + 5.*a1*h2)/4480.;
    answer.at(1) = fac1*md.at(1)*gamd;
    answer.at(2) = fac1*md.at(2)*gamd;
    answer.at(3) = fac1*md.at(3)*gamd;
    answer.at(4) = fac2*md.at(1)*gamd;
    answer.at(5) = fac2*md.at(2)*gamd;
    answer.at(6) = fac2*md.at(3)*gamd;
    answer.at(7) = fac3*m.dotProduct(md)*gamd;


/* Salars Matlab code

%convective mass

cmc=[...
    (h^3*(20*a3 + 3*a1*h^2)*md*gamd)/240.;...
    (h^5*md*(56*a2 + 28*a3*gam + 5*a1*h^2*gam)*gamd)/4480.;...
    (h^5*(28*a3 + 5*a1*h^2)*m'*md*gamd)/4480.];
     
cmass=rho*bc'*cmc;
*/
}

void
TrDirShell :: computeThicknessMappingCoeff(FloatArray &answer){
	//thickness jacobian = ratio between volume and area: j0 = a3 + a2*zeta^2 + a1 * zeta
	// Returns array with a1-a3, used in expression for analytical integration of mass matrix.
    IntegrationRule *iRule = integrationRulesArray [ 0 ];
	GaussPoint *gp;
	gp = iRule->getIntegrationPoint(0);
    FloatArray lcoords = *gp->giveCoordinates();
	
    FloatMatrix dNdxi;
    this->evaldNdxi(dNdxi, lcoords);

	FloatArray M, dM1(3), dM2(3), dX1(3), dX2(3);
    this->evalInitialDirectorAt(gp, M); 
    //dX1.resize(3); dX2.resize(3); dM1.resize(3); dM2.resize(3);

	for (int i = 1; i <= 6; i++ ) {
		double x, y, z, Mx, My, Mz;
        FloatArray *nodeI = this->giveNode(i)->giveCoordinates();
		x = nodeI->at(1); y = nodeI->at(2); z = nodeI->at(3);
        
		M=this->giveInitialNodeDirector(i);
		Mx = M.at(1); My = M.at(2); Mz = M.at(3);

        dX1.at(1) += dNdxi.at(i,1)* x; 
		dX1.at(2) += dNdxi.at(i,1)* y; 
		dX1.at(3) += dNdxi.at(i,1)* z;
        dM1.at(1) += dNdxi.at(i,1)* Mx; 
		dM1.at(2) += dNdxi.at(i,1)* My; 
		dM1.at(3) += dNdxi.at(i,1)* Mz; 

        dX2.at(1) += dNdxi.at(i,2)* x; 
		dX2.at(2) += dNdxi.at(i,2)* y; 
		dX2.at(3) += dNdxi.at(i,2)* z;
        dM2.at(1) += dNdxi.at(i,2)* Mx; 
		dM2.at(2) += dNdxi.at(i,2)* My; 
		dM2.at(3) += dNdxi.at(i,2)* Mz; 

    }

    double sc;
    FloatArray temp, temp2;
    temp.beVectorProductOf(dX1,dX2);
    sc = temp.computeNorm();
    answer.resize(3);
    answer.at(3) = M.dotProduct(temp)/sc;
    
    temp.beVectorProductOf(dX1,dM2); 
    temp2.beVectorProductOf(dX2,dM1); 
    temp.add(temp2);
    answer.at(2) = M.dotProduct(temp)/sc;

    temp.beVectorProductOf(dM1,dM2); 
    answer.at(1) = M.dotProduct(temp)/sc;

}

double 
TrDirShell :: computeVolumeAround(GaussPoint *gp){
	FloatArray G1, G2, G3, temp;
	double detJ;
	this->evalInitialCovarBaseVectorsAt(gp,G1, G2, G3);
	temp.beVectorProductOf(G1, G2);
	detJ = temp.dotProduct(G3) * this->giveCrossSection()->give(CS_Thickness)*0.5;
    return detJ * gp->giveWeight() ;
}




} // end namespace oofem
