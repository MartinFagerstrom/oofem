/* $Header: /home/cvs/bp/oofem/oofemlib/src/matrix.C,v 1.5.4.1 2004/04/05 15:19:43 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



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
/*
 The original idea for this class comes from 
  Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 PhD Thesis, EPFL, Lausanne, 1992.
*/


//   file MATRIX.CC

#include "matrix.h"
#include "error.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <stdlib.h>
#endif


void  Matrix :: checkBounds (int i,int j) const
   // Checks that the receiver includes a position (i,j).
{
   if (i<=0) {
     OOFEM_ERROR2 ("Matrix::checkBounds : matrix error on rows : %d < 0",i) ;
   }
   if (j<=0) {
     OOFEM_ERROR2 ("Matrix::checkBounds : matrix error on columns : %d < 0",j) ;
   }
   if (i>nRows) {
     OOFEM_ERROR3 ("Matrix::checkBounds : matrix error on rows : %d > %d",i,nRows) ;
   }
   if (j>nColumns) {
     OOFEM_ERROR3 ("Matrix::checkBounds : matrix error on columns : %d > %d",j,nColumns) ;
   }
}

