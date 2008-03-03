/* $Header: /home/cvs/bp/oofem/oofemlib/src/inputrecord.C,v 1.4.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//
// file inputrecord.cc
//

#include "inputrecord.h"
#ifndef __MAKEDEPEND
#include <ctype.h>
#endif

InputRecord::InputRecord ()  
{}

InputRecord::InputRecord (const InputRecord& src)
{}



InputRecord&
InputRecord::operator= (const InputRecord& src)
{
 return *this;
}

const char*
InputRecord::strerror (IRResultType rt)
{
 switch(rt) {
 case IRRT_NOTFOUND:
  return "Missing Keyword"; // string literal is statically allocated, return safe
 case IRRT_BAD_FORMAT:
  return "Bad format";
 default:
  return "Unknown error";
 }
}

void 
InputRecord::report_error (const char* _class, const char* proc, const InputFieldType fieldID, const char* kwd,
              IRResultType result, const char* file, int line)
{
  __OOFEM_ERROR6 (file, line, "Input error: \"%s\", field keyword \"%s\" (fieldID=%d)\n%s::%s", strerror(result), kwd, fieldID, _class, proc);
}
