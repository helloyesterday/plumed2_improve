/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_tools_NeighborList2D_h
#define __PLUMED_tools_NeighborList2D_h

#include "Vector.h"
#include "AtomNumber.h"

#include <vector>

namespace PLMD {

class Pbc;

/// \ingroup TOOLBOX
/// A class that implements neighbor lists from two lists or a single list of atoms
class NeighborList2D
{
  bool reduced;
  bool do_pair_,do_pbc_,twolists_;
  const PLMD::Pbc* pbc_;
  std::vector<PLMD::AtomNumber> fullatomlist_,requestlist_;
  std::vector<std::pair<unsigned,unsigned> > neighbors_;
  double distance_;
// added by Yi Isaac Yang
  unsigned axis_id,map_id1,map_id2;
  double axis_dis;
// end add
  unsigned stride_,nlist0_,nlist1_,nallpairs_,lastupdate_;
/// Initialize the neighbor list with all possible pairs
  void initialize();
/// Return the pair of indexes in the positions array
/// of the two atoms forming the i-th pair among all possible pairs
  std::pair<unsigned,unsigned> getIndexPair(unsigned i);
/// Extract the list of atoms from the current list of close pairs
  void setRequestList();
public:
  NeighborList2D(const std::vector<PLMD::AtomNumber>& list0,
               const std::vector<PLMD::AtomNumber>& list1,
               const bool& do_pair, const bool& do_pbc, const PLMD::Pbc& pbc,
               const double& distance=1.0e+30,
               const unsigned& _axis_id=2,const double& _axis_cut=1e+30,
               const unsigned& stride=0);
  NeighborList2D(const std::vector<PLMD::AtomNumber>& list0, const bool& do_pbc,
               const PLMD::Pbc& pbc, const double& distance=1.0e+30,
               const unsigned& _axis_id=2,const double& _axis_cut=1e+30,
               const unsigned& stride=0);
/// Return the list of all atoms. These are needed to rebuild the neighbor list.
  std::vector<PLMD::AtomNumber>& getFullAtomList();
/// Update the indexes in the neighbor list to match the
/// ordering in the new positions array
/// and return the new list of atoms that must be requested to the main code
  std::vector<PLMD::AtomNumber>& getReducedAtomList();
/// Update the neighbor list and prepare the new
/// list of atoms that will be requested to the main code
  void update(const std::vector<PLMD::Vector>& positions);
/// Get the update stride of the neighbor list
  unsigned getStride() const;
/// Get the last step in which the neighbor list was updated
  unsigned getLastUpdate() const;
/// Set the step of the last update
  void setLastUpdate(unsigned step);
/// Get the size of the neighbor list
  unsigned size() const;
/// Get the i-th pair of the neighbor list
  std::pair<unsigned,unsigned> getClosePair(unsigned i) const;
/// Get the list of neighbors of the i-th atom
  std::vector<unsigned> getNeighbors(unsigned i);
  ~NeighborList2D() {}
};

}

#endif
