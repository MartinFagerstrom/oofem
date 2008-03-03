/* $Header: $ */
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

#ifndef fastmarchingmethod_h
#include "domain.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
#include <vector>
#include <list>
#include <queue>
#endif
/**
   Fast Marching Method for unstructured grids. 
   Used to solve Eikonal equation and especially to construct 
   signed distance function.
 */
class FastMarchingMethod 
{
 protected:
  
  /**
     Type describing node status for fast marching method.
     FMM_Status_FAR - nodes not yet visited 
     FMM_Status_TRIAL - trial nodes, candidates for known (accepted) 
     FMM_Status_KNOWN - accepted nodes
     FMM_Status_KNOWN_BOUNDARY - boundary nodes, from which the front will not propagate
   */
  enum FNM_Status_Type {FMM_Status_FAR, FMM_Status_TRIAL, FMM_Status_KNOWN, FMM_Status_KNOWN_BOUNDARY};
  /// DofManager Fast Marching data record;
  class FMM_DofmanRecord {
  public:
    FNM_Status_Type status;
  };

  /// Array of DofManager records
  std::vector<FMM_DofmanRecord> dmanRecords;
  /// Pointer to working set of dmanValues
  const FloatArray *dmanValuesPtr;
  /// Delegate of FMM_DofmanRecord; stored in priority queue
  //class FMM_DofmanRecordDelegate {
  //public:FMM_DofmanRecordDelegaFMM_DofmanRecordDelegatete
  //  int id;
  //};
  
  struct FMM_DofmanRecordDelegate_greater {
    const FloatArray** dmanValuesPtrRef;
    FMM_DofmanRecordDelegate_greater (const FloatArray** dv) : dmanValuesPtrRef(dv)
    {}
    bool operator() (const int& p, const int& q) const 
      {return (fabs((*dmanValuesPtrRef)->at(p)) > fabs((*dmanValuesPtrRef)->at(q)));}
  };

  /// Domain
  Domain* domain;

  /// Priority queue for trial T values
  std::priority_queue<int, std::vector<int>, FMM_DofmanRecordDelegate_greater > dmanTrialQueue;

 public:
  /** Constructor. Takes two two arguments. Creates 
      MaterialInterface instance with given number and belonging to given domain.
      @param n component number in particular domain. For instance, can represent  
      node number in particular domain.
      @param d domain to which component belongs to 
  */
  FastMarchingMethod (Domain* d):dmanTrialQueue(FMM_DofmanRecordDelegate_greater(&this->dmanValuesPtr)) { domain = d;}
  ~FastMarchingMethod() {}

  /** solution of problem. I/O param dmanValues on input will contain boundary
      values for those dofnam, that are known; on output will contain solution.
      paramemeter bcDofMans is a list containing IDs (numbers) of those 
      dofmans, for which boundary value is known. If this number is positive, 
      then the front will propagate from this dofman, if negative, then the front 
      will not propagate from this dofman (usefull, when one needs to construct 
      "one sided" solution).
      F is the front propagation speed.
  */
  void solve (FloatArray &dmanValues, const std::list<int>& bcDofMans, double F);

  // identification
  const char* giveClassName () const { return "FastMarchingMethod";}
  classType giveClassID ()     const { return FastMarchingMethodClass;}


 protected:
  /// initialize receiver
  void initialize (FloatArray &dmanValues, const std::list<int>& bcDofMans, double F);

  // updates the distace of trial node with given id)
  void updateTrialValue (FloatArray &dmanValues, int id, double F);

  /// get the trial point with smallest T; zero if emty
  int  getSmallestTrialDofMan ();
};

#define fastmarchingmethod_h
#endif