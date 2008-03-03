/* $Header: /home/cvs/bp/oofem/oofemlib/src/feinterpol2d.h,v 1.1 2003/04/06 14:08:24 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2002   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

//   ************************************
//   *** CLASS FEInterpolation2d ***
//   ************************************


#ifndef feinterpol2d_h


#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "feinterpol.h"

/*
 Class representing a general abstraction for surface finite element interpolation class.
*/
class FEInterpolation2d : public FEInterpolation
{

protected:
public:
 FEInterpolation2d (int o) : FEInterpolation (o) {}

 /**@name Edge interpolation services */
 //@{
 virtual void computeLocalEdgeMapping (IntArray& edgeNodes, int iedge) = 0;
 void computeEdgeMapping (IntArray& edgeNodes, IntArray& elemNodes, int iedge) {
   int i, size;
   IntArray ln;
   this->computeLocalEdgeMapping (ln, iedge);
   size = ln.giveSize();
   edgeNodes.resize(size);
   for (i=1; i<=size; i++) edgeNodes.at(i)=elemNodes.at(ln.at(i));
 }
 /** 
  Evaluates the array of edge interpolation functions (shape functions) at given point.
  @param answer contains resulting array of evaluated interpolation functions
  @param lcoords array containing (local) coordinates 
  @param time time
  */
 virtual void edgeEvalN (FloatArray& answer, const FloatArray& lcoords, double time) = 0;
 /** 
  Evaluates the matrix of derivatives of edge interpolation functions (shape functions) at given point.
  These derivatives are in global coordinate system (where the nodal coordinates are defined) 
  @param answer contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj
  @param iedge determines the edge number
  @param nodes array of node numbers (for the whole element) defining the interpolation geometry
  @param lcoords array containing (local) coordinates 
  @param time time
  */
 virtual void edgeEvaldNdx (FloatMatrix&answer, int iedge, 
                            Domain* d, IntArray& nodes, const FloatArray& lcoords, double time) = 0;
 /** 
  Evaluates the matrix of derivatives of edge interpolation functions (shape functions) at given point.
  These derivatives are in global coordinate system (where the nodal coordinates are defined) 
  @param answer contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj
  @param iedge determines the edge number
  @param coords coordinates of nodes defining the interpolation geometry (for the whole element)
  @param lcoords array containing (local) coordinates 
  @param time time
  */
 virtual void edgeEvaldNdx (FloatMatrix&answer, int iedge, 
                            const FloatArray** coords, const FloatArray& lcoords, double time) = 0;
 /** 
  Evaluates edge global coordinates from given local ones
  These derivatives are in global coordinate system (where the nodal coordinates are defined) 
  @param answer contains resulting global coordinates
  @param iedge determines edge number
  @param nodes array of node numbers (for the whole element) defining the interpolation geometry
  @param lcoords array containing (local) coordinates 
  @param time time
  */
 virtual void edgeLocal2global (FloatArray& answer, int iedge,
                 Domain* d, IntArray& nodes, const FloatArray& lcoords, double time) = 0;
 /** 
  Evaluates edge global coordinates from given local ones
  These derivatives are in global coordinate system (where the nodal coordinates are defined) 
  @param answer contains resulting global coordinates
  @param iedge determines edge number
  @param coords coordinates of nodes defining the interpolation geometry (for the whole element) 
  @param lcoords array containing (local) coordinates 
  @param time time
  */
 virtual void edgeLocal2global (FloatArray& answer, int iedge,
                                const FloatArray** coords, const FloatArray& lcoords, double time) = 0;
 /**
  Evaluates the edge jacobian of transformation between local and global coordinates.
  */
 virtual double edgeGiveTransformationJacobian (int iedge, const FloatArray **coords, const FloatArray& lcoords, 
                                                double time) = 0;
 /**
  Evaluates the edge jacobian of transformation between local and global coordinates.
  */
 virtual double edgeGiveTransformationJacobian (int iedge, Domain* d, IntArray& nodes, const FloatArray& lcoords, 
                                                double time) = 0;
 //@}

} ;




#define feinterpol2d_h
#endif





