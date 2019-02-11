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
#include "NeighborList2DParallel.h"
#include "Vector.h"
#include "Pbc.h"
#include "AtomNumber.h"
#include "Tools.h"
#include <vector>
#include <algorithm>
#include "Communicator.h"
#include "Log.h"

namespace PLMD {
using namespace std;

NeighborList2DParallel::NeighborList2DParallel(const vector<AtomNumber>& list0, const vector<AtomNumber>& list1,
                           const bool& do_pair, const bool& do_pbc, const Pbc& pbc, Communicator& cc, Log& log,
                           const unsigned& _axis_id,const double& distance,const double& _axis_dis,
                           const bool& do_full_list, const int& stride, const double& skin, const double& _axis_skin): 
  do_pair_(do_pair), do_pbc_(do_pbc), pbc_(&pbc),axis_id(_axis_id),
  distance_(distance),axis_dis(_axis_dis),
  stride_(stride), mycomm(cc), mylog(log),
  skin_(skin),axis_skin(_axis_skin), do_full_list_(do_full_list)
{
// store full list of atoms needed
  if(axis_id==0) {map_id1=1;map_id2=2;}
  else if(axis_id==1) {map_id1=0;map_id2=2;}
  else if(axis_id==2) {map_id1=0;map_id2=1;}
  else plumed_merror("axis_id must be 0, 1 or 2");
  fullatomlist_=list0;
  fullatomlist_.insert(fullatomlist_.end(),list1.begin(),list1.end());
  nlist0_=list0.size();
  nlist1_=list1.size();
  twolists_=true;
  if(!do_pair) {
    nallpairs_=nlist0_*nlist1_;
  } else {
    plumed_assert(nlist0_==nlist1_);
    nallpairs_=nlist0_;
  }
  lastupdate_=0;
  positions_old_.resize(fullatomlist_.size());
  dangerousBuilds_=0;
  firsttime_=true;
  numberOfBuilds_=0;
  avgTotalNeighbors_=0.;
  maxLoadImbalance_=2.;
  avgLoadImbalance_=0.;
  // Setup mpi
  mpi_rank=mycomm.Get_rank();
  mpi_stride=mycomm.Get_size();
}

NeighborList2DParallel::NeighborList2DParallel(const vector<AtomNumber>& list0, const bool& do_pbc,
                           const Pbc& pbc, Communicator& cc, Log& log, const unsigned& _axis_id,
                           const double& distance, const double& _axis_dis,
                           const bool& do_full_list, const int& stride,
                           const double& skin, const double& _axis_skin):
  do_pbc_(do_pbc), pbc_(&pbc),axis_id(_axis_id),
  distance_(distance),axis_dis(_axis_dis),
  stride_(stride), mycomm(cc), mylog(log),
  skin_(skin) , axis_skin(_axis_skin), do_full_list_(do_full_list)
{
  if(axis_id==0) {map_id1=1;map_id2=2;}
  else if(axis_id==1) {map_id1=0;map_id2=2;}
  else if(axis_id==2) {map_id1=0;map_id2=1;}
  else plumed_merror("axis_id must be 0, 1 or 2");
  fullatomlist_=list0;
  nlist0_=list0.size();
  twolists_=false;
  nallpairs_=nlist0_*(nlist0_-1)/2;
  lastupdate_=0;
  positions_old_.resize(fullatomlist_.size());
  dangerousBuilds_=0;
  firsttime_=true;
  numberOfBuilds_=0;
  avgTotalNeighbors_=0.;
  maxLoadImbalance_=2.;
  avgLoadImbalance_=0.;
  // Setup mpi
  mpi_rank=mycomm.Get_rank();
  mpi_stride=mycomm.Get_size();
}

vector<AtomNumber>& NeighborList2DParallel::getFullAtomList() {
  return fullatomlist_;
}

bool NeighborList2DParallel::isListStillGood(const vector<Vector>& positions) {
  bool flag=true;
  plumed_assert(positions.size()==fullatomlist_.size());
  double skin2=skin_*skin_;
  for(unsigned int i=0;i<fullatomlist_.size();i++) {
    Vector distance;
    if(do_pbc_) {
       distance=pbc_->distance(positions[i],positions_old_[i]);
    } else {
       distance=delta(positions[i],positions_old_[i]);
    }
    double mod2_dis = distance[map_id1]*distance[map_id1]+distance[map_id2]*distance[map_id2];
    if (mod2_dis>skin2||distance[axis_id]>axis_skin) {
       flag=false;
       break;
    }
  }
  return flag;
}

void NeighborList2DParallel::printStats() {
  mylog.printf("Neighbor list 2D statistics\n");
  if(do_full_list_)
    mylog.printf("  with full neighbor list\n");
  else
    mylog.printf("  with half neighbor list\n");
  mylog.printf("  with mapping axis %d and update step %d\n", int(axis_id),int(stride_));
  mylog.printf("  with plane cutoff %f and axis cutoff %f\n", distance_,axis_dis);
  mylog.printf("  with plane skin %f and axis skin %f\n", skin_,axis_skin);
  mylog.printf("  Total # of neighbors = %f \n", avgTotalNeighbors_);
  mylog.printf("  Ave neighs/atom = %f \n", avgTotalNeighbors_ /(double) nlist0_);
  mylog.printf("  Neighbor list builds = %d \n",numberOfBuilds_);
  mylog.printf("  Dangerous builds = %d \n",dangerousBuilds_);
  mylog.printf("  Average load imbalance (min/max) = %f \n",avgLoadImbalance_);
  mylog.printf("  Maximum load imbalance (min/max) = %f \n",maxLoadImbalance_);
  if (do_linked_list_) {
     mylog.printf("  Number of bins in linked list = %d %d %d \n", nbinx, nbiny, nbinz);
  }
}

void NeighborList2DParallel::update(const vector<Vector>& positions) {
  // clear previous list
  neighbors_.clear();
  local_atoms_.clear();
  // check if positions array has the correct length
  plumed_assert(positions.size()==fullatomlist_.size());
  // Prepare linked lists
  prepareLinkedList();
  // Decide whether to do linked lists or the N^2 calculation
  do_linked_list_=true;
  if ((2*sx+1) >= nbinx || (2*sy+1) >= nbiny || (2*sz+1) >= nbinz) do_linked_list_=false;
  if (do_linked_list_) {
     if (!do_full_list_) updateHalfListWithLinkedList(positions);
     else updateFullListWithLinkedList(positions);
  } else {
     if (!do_full_list_) updateHalfList(positions);
     else updateFullList(positions);
  }
  gatherStats(positions);
  // Store positions for checking
  for(unsigned int i=0;i<fullatomlist_.size();i++) {
     positions_old_[i]=positions[i];
  }
}

void NeighborList2DParallel::updateFullListWithLinkedList(const vector<Vector>& positions) {
  const double d2=distance_*distance_;
  // Set bin head of chain to -1.0
  binhead = std::vector<int>(nbinx*nbiny*nbinz,-1.0);
  if (!twolists_) {
    neighbors_.resize(nlist0_);
    // Construct linked list and assign head of chain
    bins.resize(nlist0_);
    for(unsigned int i=0;i<nlist0_;i+=1) {
       unsigned ibin = coord2bin(positions[i]);
       bins[i] = binhead[ibin];
       binhead[ibin] = i;
    }
    // Calculate
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       unsigned atombinx, atombiny, atombinz, atombin;
       atombin = coord2bin(positions[i],atombinx, atombiny, atombinz);
       for(unsigned k=0; k < nstencil; k++) {
          // Loop over atoms in the k-th neighboring cell
          unsigned kbin = neighborCell(k,atombinx, atombiny, atombinz) ;
          for (int j = binhead[kbin]; j >= 0; j = bins[j]) {
			 Vector dis=distance(position_i,positions[j]);
             double value=dis[map_id1]*dis[map_id1]+dis[map_id2]*dis[map_id2];
             double ax_value=std::fabs(dis[axis_id]);
             if(value<=d2&&ax_value<=axis_dis) neighbors_[i].push_back(j);
          }
       }
    }
  } else if (twolists_ && !do_pair_) {
    neighbors_.resize(nlist0_+nlist1_);
    // Construct linked list and assign head of chain
    bins.resize(nlist1_);
    for(unsigned int i=0;i<nlist1_;i++) {
       unsigned ibin = coord2bin(positions[i+nlist0_]);
       bins[i] = binhead[ibin];
       binhead[ibin] = i;
    }
    // Calculate
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       unsigned atombinx, atombiny, atombinz, atombin;
       atombin = coord2bin(positions[i],atombinx, atombiny, atombinz);
       for(unsigned k=0; k < nstencil; k++) {
          // Loop over atoms in the k-th neighboring cell
          unsigned kbin = neighborCell(k,atombinx, atombiny, atombinz) ;
          for (int j = binhead[kbin]; j >= 0; j = bins[j]) {
             Vector dis=distance(position_i,positions[j+nlist0_]);
             double value=dis[map_id1]*dis[map_id1]+dis[map_id2]*dis[map_id2];
             double ax_value=std::fabs(dis[axis_id]);
             if(value<=d2&&ax_value<=axis_dis) {
                neighbors_[i].push_back(j+nlist0_);
                neighbors_[j+nlist0_].push_back(i);
             }
          }
       }
    }
  }
}

void NeighborList2DParallel::updateHalfListWithLinkedList(const vector<Vector>& positions) {
  const double d2=distance_*distance_;
  neighbors_.resize(nlist0_);
  // Set bin head of chain to -1.0
  binhead = std::vector<int>(nbinx*nbiny*nbinz,-1.0);
  if (!twolists_) {
    // Construct linked list and assign head of chain
    bins.resize(nlist0_);
    for(unsigned int i=0;i<nlist0_;i+=1) {
       unsigned ibin = coord2bin(positions[i]);
       bins[i] = binhead[ibin];
       binhead[ibin] = i;
    }
    // Calculate
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       unsigned atombinx, atombiny, atombinz, atombin;
       atombin = coord2bin(positions[i],atombinx, atombiny, atombinz);
       for(unsigned k=0; k < nstencil; k++) {
          // Loop over atoms in the k-th neighboring cell
          unsigned kbin = neighborCell(k,atombinx, atombiny, atombinz) ;
          for (int j = binhead[kbin]; j >= 0; j = bins[j]) {
             if (j>=i) continue;
             Vector dis=distance(position_i,positions[j]);
             double value=dis[map_id1]*dis[map_id1]+dis[map_id2]*dis[map_id2];
             double ax_value=std::fabs(dis[axis_id]);
             if(value<=d2&&ax_value<=axis_dis) neighbors_[i].push_back(j);
          }
       }
    }
  } else if (twolists_ && !do_pair_) {
    // Construct linked list and assign head of chain
    bins.resize(nlist1_);
    for(unsigned int i=0;i<nlist1_;i++) {
       unsigned ibin = coord2bin(positions[i+nlist0_]);
       bins[i] = binhead[ibin];
       binhead[ibin] = i;
    }
    // Calculate
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       unsigned atombinx, atombiny, atombinz, atombin;
       atombin = coord2bin(positions[i],atombinx, atombiny, atombinz);
       for(unsigned k=0; k < nstencil; k++) {
          // Loop over atoms in the k-th neighboring cell
          unsigned kbin = neighborCell(k,atombinx, atombiny, atombinz) ;
          for (int j = binhead[kbin]; j >= 0; j = bins[j]) {
             Vector dis=distance(position_i,positions[j+nlist0_]);
             double value=dis[map_id1]*dis[map_id1]+dis[map_id2]*dis[map_id2];
             double ax_value=std::fabs(dis[axis_id]);
             if(value<=d2&&ax_value<=axis_dis) {
                neighbors_[i].push_back(j+nlist0_);
             }
          }
       }
    }
  }
}


void NeighborList2DParallel::updateFullList(const vector<Vector>& positions) {
  const double d2=distance_*distance_;
  if (!twolists_) {
    neighbors_.resize(nlist0_);
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       for(unsigned int j=0;j<nlist0_;j+=1) {
          Vector dis=distance(position_i,positions[j]);
          double value=dis[map_id1]*dis[map_id1]+dis[map_id2]*dis[map_id2];
          double ax_value=std::fabs(dis[axis_id]);
          if(value<=d2&&ax_value<=axis_dis) neighbors_[i].push_back(j);
       }
    }
  } else if(twolists_ && do_pair_) {
    neighbors_.resize(nlist0_+nlist1_);
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector dis=distance(positions[i],positions[nlist0_+i]);
       double value=dis[map_id1]*dis[map_id1]+dis[map_id2]*dis[map_id2];
       double ax_value=std::fabs(dis[axis_id]);
       if(value<=d2&&ax_value<=axis_dis) {
          neighbors_[i].push_back(nlist0_+i);
          neighbors_[nlist0_].push_back(i);
       }
    }
  } else if (twolists_ && !do_pair_) {
    neighbors_.resize(nlist0_+nlist1_);
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       for(unsigned int j=0;j<nlist1_;j+=1) {
          Vector dis=distance(position_i,positions[nlist0_+j]);
          double value=dis[map_id1]*dis[map_id1]+dis[map_id2]*dis[map_id2];
          double ax_value=std::fabs(dis[axis_id]);
          if(value<=d2&&ax_value<=axis_dis) {
             neighbors_[i].push_back(nlist0_+j);
             neighbors_[nlist0_+j].push_back(i);
          }
       }
    }
  }
}

void NeighborList2DParallel::updateHalfList(const vector<Vector>& positions) {
  const double d2=distance_*distance_;
  neighbors_.resize(nlist0_);
  if (!twolists_) {
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       for(unsigned int j=i+1;j<nlist0_;j+=1) {
          Vector dis=distance(position_i,positions[j]);
          double value=dis[map_id1]*dis[map_id1]+dis[map_id2]*dis[map_id2];
          double ax_value=std::fabs(dis[axis_id]);
          if(value<=d2&&ax_value<=axis_dis) neighbors_[i].push_back(j);
       }
    }
  } else if(twolists_ && do_pair_) {
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector dis=distance(positions[i],positions[nlist0_+i]);
       double value=dis[map_id1]*dis[map_id1]+dis[map_id2]*dis[map_id2];
       double ax_value=std::fabs(dis[axis_id]);
       if(value<=d2&&ax_value<=axis_dis) neighbors_[i].push_back(nlist0_+i);
    }
  } else if (twolists_ && !do_pair_) {
    for(unsigned int i=mpi_rank;i<nlist0_;i+=mpi_stride) {
       local_atoms_.push_back(i);
       Vector position_i=positions[i];
       for(unsigned int j=0;j<nlist1_;j++) {
          Vector dis=distance(position_i,positions[nlist0_+j]);
          double value=dis[map_id1]*dis[map_id1]+dis[map_id2]*dis[map_id2];
          double ax_value=std::fabs(dis[axis_id]);
          if(value<=d2&&ax_value<=axis_dis) {
             neighbors_[i].push_back(nlist0_+j);
          }
       }
    }
  }
}

void NeighborList2DParallel::gatherStats(const vector<Vector>& positions) {
  // Check if rebuilt was dangerous
  if (!firsttime_ && !isListStillGood(positions)) {
     dangerousBuilds_++;
  }
  firsttime_=false;
  numberOfBuilds_++;
  unsigned neighNum=0;
  for(unsigned int i=0;i<nlist0_;i++) {
     neighNum += neighbors_[i].size();
  }
  unsigned allNeighNum=0;
  std::vector<unsigned> neighbors_ranks_(mycomm.Get_size());
  mycomm.Allgather(&neighNum,1,&neighbors_ranks_[0],1);
  for(unsigned int i=0;i<mycomm.Get_size();i+=1) allNeighNum+=neighbors_ranks_[i];
  auto min_element_ = *std::min_element(neighbors_ranks_.begin(), neighbors_ranks_.end());
  auto max_element_ = *std::max_element(neighbors_ranks_.begin(), neighbors_ranks_.end());
  double loadImbalance=min_element_ / (double) max_element_;
  if (maxLoadImbalance_>loadImbalance) maxLoadImbalance_=loadImbalance;
  avgLoadImbalance_ += (loadImbalance-avgLoadImbalance_)/numberOfBuilds_;
  avgTotalNeighbors_ += (allNeighNum-avgTotalNeighbors_)/numberOfBuilds_;
}

unsigned NeighborList2DParallel::getNumberOfLocalAtoms() const {
  return local_atoms_.size();
}

unsigned NeighborList2DParallel::getIndexOfLocalAtom(unsigned i) const {
  return local_atoms_[i];
}

int NeighborList2DParallel::getStride() const {
  return stride_;
}

unsigned NeighborList2DParallel::getLastUpdate() const {
  return lastupdate_;
}

void NeighborList2DParallel::setLastUpdate(unsigned step) {
  lastupdate_=step;
}

// This function is obsolete and should be changed
unsigned NeighborList2DParallel::size() const {
  return neighbors_.size();
}

// This function is obsolete and should be changed
pair<unsigned,unsigned> NeighborList2DParallel::getClosePair(unsigned i) const {
  return pair<unsigned,unsigned>(i,i);
}

vector<unsigned> NeighborList2DParallel::getNeighbors(unsigned index) {
  return neighbors_[index];
}

Vector NeighborList2DParallel::distance(const Vector& position1, const Vector& position2) {
  Vector distance;
  if(do_pbc_) {
    distance=pbc_->distance(position1,position2);
  } else {
    distance=delta(position1,position2);
  }
  return distance;
}

// Linked lists stuff

void NeighborList2DParallel::prepareLinkedList() {
  // Determine optimal number of cells in each direction
  double binsize_optimal = 0.5*distance_;
  double binsize_axis = 0.5*axis_dis;
  double binsizeinv = 1.0/binsize_optimal;
  double binaxinv = 1.0/binsize_axis;
  Tensor bbox=pbc_->getBox();
  nbinx = static_cast<int> (bbox[0][0]*binsizeinv);
  nbiny = static_cast<int> (bbox[1][1]*binsizeinv);
  nbinz = static_cast<int> (bbox[2][2]*binsizeinv);
  
  if(axis_id==0) {nbinx = static_cast<int> (bbox[0][0]*binaxinv);}
  else if(axis_id==1) {nbiny = static_cast<int> (bbox[1][1]*binaxinv);}
  else if(axis_id==2) {nbinz = static_cast<int> (bbox[2][2]*binaxinv);}
  
  if (nbinx == 0) nbinx = 1;
  if (nbiny == 0) nbiny = 1;
  if (nbinz == 0) nbinz = 1;
  // Determine cell sizes
  binsizex = 1.0 / nbinx;
  binsizey = 1.0 / nbiny;
  binsizez = 1.0 / nbinz;
  bininvx = 1.0 / binsizex;
  bininvy = 1.0 / binsizey;
  bininvz = 1.0 / binsizez;
  // sx,sy,sz = max range of stencil in each dim
  Vector scaled_distance_=pbc_->realToScaled(Vector(distance_,distance_,distance_));
  scaled_distance_[axis_id]=axis_dis;
  sx = static_cast<int> (scaled_distance_[0]*bininvx);
  if (sx*binsizex < scaled_distance_[0]) sx++;
  sy = static_cast<int> (scaled_distance_[1]*bininvy);
  if (sy*binsizey < scaled_distance_[1]) sy++;
  sz = static_cast<int> (scaled_distance_[2]*bininvz);
  if (sz*binsizez < scaled_distance_[2]) sz++;
  // Create stencil
  nstencil = (2*sx+1)*(2*sy+1)*(2*sz+1);
  stencilx.resize(nstencil);
  stencily.resize(nstencil);
  stencilz.resize(nstencil);
  unsigned k=0;
  for(int ix=-sx;ix<=sx;ix++) {
     for(int iy=-sy;iy<=sy;iy++) {
        for(int iz=-sz;iz<=sz;iz++) {
           stencilx[k]=ix;
           stencily[k]=iy;
           stencilz[k]=iz;
           k++;
        }
     }
  }
   // clear previous list
  bins.clear();
  binhead.clear();
}

unsigned NeighborList2DParallel::neighborCell(const unsigned& k, const unsigned& atombinx, const unsigned& atombiny, const unsigned& atombinz) {
  int ibinx = atombinx+stencilx[k];
  if (ibinx<0) ibinx += nbinx; 
  if (ibinx>=nbinx) ibinx -= nbinx;
  int ibiny = atombiny+stencily[k];
  if (ibiny<0) ibiny += nbiny; 
  if (ibiny>=nbiny) ibiny -= nbiny;
  int ibinz = atombinz+stencilz[k];
  if (ibinz<0) ibinz += nbinz; 
  if (ibinz>=nbinz) ibinz -= nbinz;
  unsigned neighbor = ibinx + ibiny*nbinx + ibinz*nbinx*nbiny;
  return neighbor; 
}

unsigned NeighborList2DParallel::coord2bin(const Vector& position) {
  Vector pbc_position=pbc_->distance(Vector(0.0,0.0,0.0),position);
  Vector scaled_position=pbc_->realToScaled(pbc_position);
  if (scaled_position[0]<0.) scaled_position[0] += 1;
  if (scaled_position[1]<0.) scaled_position[1] += 1;
  if (scaled_position[2]<0.) scaled_position[2] += 1;
  unsigned ibinx, ibiny, ibinz, ibin;
  ibinx = static_cast<int> (scaled_position[0]*bininvx);
  ibiny = static_cast<int> (scaled_position[1]*bininvy);
  ibinz = static_cast<int> (scaled_position[2]*bininvz);
  ibin = ibinx + ibiny*nbinx + ibinz*nbinx*nbiny;
  return ibin;
}

unsigned NeighborList2DParallel::coord2bin(const Vector& position, unsigned& ibinx, unsigned& ibiny, unsigned& ibinz) {
  Vector pbc_position=pbc_->distance(Vector(0.0,0.0,0.0),position);
  Vector scaled_position=pbc_->realToScaled(pbc_position);
  if (scaled_position[0]<0.) scaled_position[0] += 1;
  if (scaled_position[1]<0.) scaled_position[1] += 1;
  if (scaled_position[2]<0.) scaled_position[2] += 1;
  unsigned ibin;
  ibinx = static_cast<int> (scaled_position[0]*bininvx);
  ibiny = static_cast<int> (scaled_position[1]*bininvy);
  ibinz = static_cast<int> (scaled_position[2]*bininvz);
  ibin = ibinx + ibiny*nbinx + ibinz*nbinx*nbiny;
  return ibin;
}

}
