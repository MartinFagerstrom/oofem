/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/processcomm.C,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
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

#ifdef __PARALLEL_MODE

#include "cltypes.h"
#include "processcomm.h"
#include "intarray.h"

#include "combuff.h"
#include "dyncombuff.h"

#ifdef __USE_MPI
#ifndef __MAKEDEPEND
#include "mpi.h"
#endif
#endif


ProcessCommunicatorBuff :: ProcessCommunicatorBuff (CommBuffType t) 
{
  if (t == CBT_static) {
    send_buff = new StaticCommunicationBuffer (MPI_COMM_WORLD);
    recv_buff = new StaticCommunicationBuffer (MPI_COMM_WORLD);
  } else {
    send_buff = new DynamicCommunicationBuffer (MPI_COMM_WORLD);
    recv_buff = new DynamicCommunicationBuffer (MPI_COMM_WORLD);
  }
}


ProcessCommunicator :: ProcessCommunicator (EngngModel* d, ProcessCommunicatorBuff* b, int rank, CommunicatorMode m)
: toSend(), toReceive()
{
 this->localProblem = d;
 this->rank         = rank;
 this->pcBuffer =  b;
 this->mode         = m;
}


int
ProcessCommunicator :: initSend (int tag)
{
 int result = 1;
 if (!toSend.isEmpty() || (this->mode == CommMode_Dynamic)) {
//  fprintf (stderr, "\nPNlDEIDynamicDomainComunicator :: initExchange: sending to %d",rank);
  result = giveProcessCommunicatorBuff()->initSend (this->rank, tag);
 } else giveProcessCommunicatorBuff()->initSendBuff ();
 return result;
}


int
ProcessCommunicator :: initReceive (int tag)
{
 int result = 1;
 if (!toReceive.isEmpty() || (this->mode == CommMode_Dynamic)) {
//  fprintf (stderr, "\nPNlDEIDynamicDomainComunicator :: initExchange: recv from %d",rank);
  result &= giveProcessCommunicatorBuff()->initReceive (this->rank, tag);
 } else giveProcessCommunicatorBuff()->initRecvBuff ();

 return result;
}


int
ProcessCommunicator :: initExchange (int tag)
{
 int result = 1;
 result &= initSend(tag);
 result &= initReceive (tag);

 return result;
}


void
ProcessCommunicator::clearBuffers ()
{
 giveProcessCommunicatorBuff()->init ();
}


#endif