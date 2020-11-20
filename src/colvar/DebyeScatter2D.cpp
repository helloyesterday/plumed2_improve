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
#include "Colvar.h"
#include "tools/NeighborList2DParallel.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "ActionRegister.h"

#include <string>
#include <memory>
#include <tr1/cmath>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR COORDINATION
/*
Calculate coordination numbers.

This keyword can be used to calculate the number of contacts between two groups of atoms
and is defined as
\f[
\sum_{i\in A} \sum_{i\in B} s_{ij}
\f]
where \f$s_{ij}\f$ is 1 if the contact between atoms \f$i\f$ and \f$j\f$ is formed,
zero otherwise.
In practise, \f$s_{ij}\f$ is replaced with a switching function to make it differentiable.
The default switching function is:
\f[
s_{ij} = \frac{ 1 - \left(\frac{{\bf r}_{ij}-d_0}{r_0}\right)^n } { 1 - \left(\frac{{\bf r}_{ij}-d_0}{r_0}\right)^m }
\f]
but it can be changed using the optional SWITCH option.

To make your calculation faster you can use a neighbor list, which makes it that only a
relevant subset of the pairwise distance are calculated at every step.

If GROUPB is empty, it will sum the \f$\frac{N(N-1)}{2}\f$ pairs in GROUPA. This avoids computing
twice permuted indexes (e.g. pair (i,j) and (j,i)) thus running at twice the speed.

Notice that if there are common atoms between GROUPA and GROUPB the switching function should be
equal to one. These "self contacts" are discarded by plumed (since version 2.1),
so that they actually count as "zero".


\par Examples

The following example instructs plumed to calculate the total coordination number of the atoms in group 1-10 with the atoms in group 20-100.  For atoms 1-10 coordination numbers are calculated that count the number of atoms from the second group that are within 0.3 nm of the central atom.  A neighbour list is used to make this calculation faster, this neighbour list is updated every 100 steps.
\plumedfile
COORDINATION GROUPA=1-10 GROUPB=20-100 R_0=0.3 NLIST NL_CUTOFF=0.5 NL_STRIDE=100
\endplumedfile

The following is a dummy example which should compute the value 0 because the self interaction
of atom 1 is skipped. Notice that in plumed 2.0 "self interactions" were not skipped, and the
same calculation should return 1.
\plumedfile
c: COORDINATION GROUPA=1 GROUPB=1 R_0=0.3
PRINT ARG=c STRIDE=10
\endplumedfile

Here's an example that shows what happens when providing COORDINATION with
a single group:
\plumedfile
# define some huge group:
group: GROUP ATOMS=1-1000
# Here's coordination of a group against itself:
c1: COORDINATION GROUPA=group GROUPB=group R_0=0.3
# Here's coordination within a single group:
x: COORDINATION GROUPA=group R_0=0.3
# This is just multiplying times 2 the variable x:
c2: COMBINE ARG=x COEFFICIENTS=2

# the two variables c1 and c2 should be identical, but the calculation of c2 is twice faster
# since it runs on half of the pairs.
PRINT ARG=c1,c2 STRIDE=10
\endplumedfile



*/
//+ENDPLUMEDOC

class DebyeScatter2D : public Colvar {
  bool pbc;
  bool serial;
  NeighborList2DParallel* nl;
  bool invalidateList;
  bool firsttime;
  bool usew;
  bool use_grid;
  bool doneigh;
  bool nl_full_list;
  
// Low communication variant
  //~ bool doLowComm;
  
  double theta;
  double maxr;
  double lambda;
  double deltad;
  double q;
  double fij;
  double fij2;
  double axis_maxr;

  unsigned qbin;
  unsigned qnum;
  unsigned axis_id;
  unsigned map_id1,map_id2;
  std::vector<double> preintensity;
  std::vector<double> prederivsf;
  
  std::string map_axis;
  
  inline void calc_value(const Vector&,double&,Vector&,Vector&,Tensor&);
  inline void calc_value(const Vector&,double&,Vector&,Tensor&);

public:
  explicit DebyeScatter2D(const ActionOptions&);
  ~DebyeScatter2D();
  static void registerKeywords( Keywords& keys );
// active methods:
  virtual void calculate();
  virtual void prepare();
  inline double pairing(double distance,double&dfunc)const;
};

PLUMED_REGISTER_ACTION(DebyeScatter2D,"DEBYE_SCATTER_2D")

void DebyeScatter2D::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.addFlag("USE_WINDOW",false,"Use a window function to revise the finite box");
  keys.addFlag("GRID_INTENSITY",false,"Use a grid to caluculate the intensity by distance");
  //~ keys.addFlag("LOW_COMM",false,"Use an algorithm with less communication between processors");
  keys.add("atoms","ATOMS","First list of atoms");
  keys.add("compulsory","THETA","Diffraction angles");
  keys.add("compulsory","LAMBDA","0.15406","Wavelength of the incident wave");
  keys.add("compulsory","MAXR","2.7","Maximum distance for the radial distribution function ");
  keys.add("compulsory","MAP_AXIS","the coordinate axis (x,y or z) to map the atom coordinates to calculate 2D XRD");
  keys.add("compulsory","AXIS_MAXR","1.6","the cutoff for mapping the coordinate axis");
  keys.add("compulsory","NL_SKIN","0.08","skin of neighbor list");
  keys.add("compulsory","NL_AXIS_SKIN","0.045","skin of neighbor list at axis");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
  keys.add("optional","NL_AXIS_CUTOFF","The cutoff for the mapping axis on the neighbour list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
  keys.add("optional","QBIN","Number of grid bins of intensity when use GRID_INTENSITY");
}

DebyeScatter2D::DebyeScatter2D(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  invalidateList(true),
  firsttime(true)
{

  parseFlag("SERIAL",serial);

  vector<AtomNumber> atoms_lista;
  parseAtomList("ATOMS",atoms_lista);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
    
  parse("LAMBDA",lambda);
  parse("THETA",theta);
  parse("MAXR",maxr);
  parse("QBIN",qbin);
  
  parse("MAP_AXIS",map_axis);
  parse("AXIS_MAXR",axis_maxr);

  if(map_axis=="x"||map_axis=="X"||map_axis=="0")
  {
	  map_axis="X";
	  axis_id=0;
	  map_id1=1;
	  map_id2=2;
  }
  else if(map_axis=="y"||map_axis=="Y"||map_axis=="1")
  {
	  map_axis="Y";
	  axis_id=1;
	  map_id1=0;
	  map_id2=2;
  }
  else if(map_axis=="z"||map_axis=="Z"||map_axis=="2")
  {
	  map_axis="Z";
	  axis_id=2;
	  map_id1=0;
	  map_id2=1;
  }
  else
    plumed_merror("Cannot understand MAP_AXIS \""+map_axis+"\", please use x, y or z");

// neighbor list stuff
  doneigh=false;
  nl_full_list=false;
  double nl_skin;
  double nl_axis_skin;
  int nl_st=0;
  parseFlag("NLIST",doneigh);
  parse("NL_SKIN",nl_skin);
  parse("NL_AXIS_SKIN",nl_axis_skin);
  
  double nl_cut=maxr+3*nl_skin;
  double nl_axis_cut=axis_maxr+3*nl_axis_skin;
  if(doneigh) {
    parse("NL_CUTOFF",nl_cut);
    if(nl_cut<=0.0) error("NL_CUTOFF ("+std::to_string(nl_cut)+") should be explicitly specified and positive");
    if(nl_cut<maxr) error("NL_CUTOFF ("+std::to_string(nl_cut)+") should not be smaller than MAXR");
    parse("NL_STRIDE",nl_st);
    if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");
    parse("NL_AXIS_CUTOFF",nl_axis_cut);
    if(nl_axis_cut<axis_maxr) error("NL_AXIS_CUTOFF should not be smaller than AXIS_MAXR");
  }
  
  parseFlag("USE_WINDOW",usew);
  parseFlag("GRID_INTENSITY",use_grid);
  
  //~ doLowComm=false;
  //~ parseFlag("LOW_COMM",doLowComm);
  //~ if (doLowComm) {
     //~ nl_full_list=true;
     //~ if (!doneigh) error("LOW_COMM can only be used with neighbor lists");
  //~ }

  addValueWithDerivatives(); setNotPeriodic();

  checkRead();
  
  // Setup neighbor list
  if (doneigh) {
    nl= new NeighborList2DParallel(atoms_lista,pbc,getPbc(),comm,log,axis_id,nl_cut,nl_axis_cut,nl_full_list,nl_st,nl_skin,nl_axis_skin);
    requestAtoms(nl->getFullAtomList());
  } else {
    requestAtoms(atoms_lista);
  }
  
  qnum=getNumberOfAtoms();
  
  deltad=maxr/qbin;
  
  q=4*pi*std::sin(theta*pi/180.0)/lambda;
  fij= 3.0485*exp(-0.132771*(q/(4*pi))*(q/(4*pi))) + 
              2.2868*exp(-0.057011*(q/(4*pi))*(q/(4*pi))) +
              1.5463*exp(-0.003239*(q/(4*pi))*(q/(4*pi))) + 
              0.8670*exp(-0.329089*(q/(4*pi))*(q/(4*pi))) + 0.2508;
  fij2 = fij*fij;

  if(use_grid)
  {
    plumed_massert(qbin>0,"QBIN must be larger than 0!");
    preintensity.assign(qbin,0);
    prederivsf.assign(qbin,0);
    for(unsigned j=0;j!=qbin;++j)
    {
      double r = deltad*(j+0.5);
      double qr = q*r;
      
      double v_ij  = std::tr1::cyl_bessel_j(0,qr);
      double dv_ij = -1.0*q*std::tr1::cyl_bessel_j(1,qr);
      
      double intensity=1;
      double devrf=0;
      //~ if(usew)
      //~ {
        //~ double w_ij=1;
        //~ double dw_ij=0;
      	//~ double pi_rc= pi*r/maxr;
      	//~ w_ij  = std::tr1::cyl_bessel_j(0,pi_rc);
      	//~ dw_ij = (pi/maxr)*std::tr1::cyl_bessel_j(1,pi_rc);
      	//~ intensity =  = 2.0 * fij2 * v_ij * w_ij / qnum;
        //~ devrf= 2.0 * fij2 * (dv_ij * w_ij + v_ij * dw_ij) / qnum;
      //~ }
      //~ else
      //~ {
        intensity = fij2 * v_ij / qnum;
        devrf = fij2 * dv_ij / qnum;;
	  //~ }

      preintensity[j] = intensity;
      prederivsf[j] = devrf;
	}
  }

  log.printf("  ATOMS:\n");
  for(unsigned int i=0; i<atoms_lista.size(); ++i) {
    if ( (i+1) % 25 == 0 ) log.printf("  \n");
    log.printf("  %d", atoms_lista[i].serial());
  }
  log.printf("  \n");
  log.printf("  with CV atom number %d\n",int(qnum),getNumberOfAtoms());
  log.printf("  with mapping axis: %s(%d)\n",map_axis.c_str(),int(axis_id));
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  if(usew) log.printf("  using window function to revise finite box\n");
  //~ if (doLowComm)log.printf("  Using the low communication variant of the algorithm");
  if(doneigh) {
    log.printf("  using neighbor lists with\n");
    log.printf("  cutoff %f, axis cutoff %f, and skin %f\n",nl_cut,nl_axis_cut,nl_skin);
    if(nl_st>=0){
      log.printf("  update every %d steps\n",nl_st);
    } else {
      log.printf("  checking every step for dangerous builds and rebuilding as needed\n");
    }
    if (nl_full_list) {
      log.printf("  using a full neighbor list\n");
    } else {
      log.printf("  using a half neighbor list\n");
    }
  }

  log.printf("  with X-ray wave length %f, angle %f and Q value %f.\n",lambda,theta,q);
  log.printf("  with atomic form factors %f.\n",fij);
  log.printf("  with plane upper limied %f and axis upper limited %f.\n", maxr,axis_maxr);
  if(use_grid) log.printf("  with integrated number %d and interval %f.\n", qbin,deltad);
  //~ for(unsigned j=0;j!=qbin;++j)
  //~ {
	//~ double r = deltad*(j+0.5);
    //~ double qr = q*r;
	//~ log.printf("    %d. with distance %f, QR %f, value %e and Derivative %e.\n", int(j),r,qr,preintensity[j],prederivsf[j]);
  //~ }
}

DebyeScatter2D::~DebyeScatter2D()
{
  if (doneigh) {
     nl->printStats();
     delete nl;
  }
}

void DebyeScatter2D::prepare() {
  if(doneigh && nl->getStride()>0){
    if(firsttime) {
      invalidateList=true;
      firsttime=false;
    } else if ( (nl->getStride()>=0) &&  (getStep()%nl->getStride()==0) ){
      invalidateList=true;
    } else if ( (nl->getStride()<0) && !(nl->isListStillGood(getPositions())) ){
      invalidateList=true;
    } else {
      invalidateList=false;
    }
  }
}

// calculator
void DebyeScatter2D::calculate()
{
  double fin_intensity=0.;
  Tensor virial;
  vector<Vector> deriv(getNumberOfAtoms());

  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) {
    stride=1;
    rank=0;
  } else {
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  unsigned nt=OpenMP::getNumThreads();

  unsigned nnum=qnum*(qnum-1)/2;
  if(doneigh)
    nnum=nl->size();
   
  const unsigned nn=nnum;

  if(nt*stride*10>nn) nt=nn/stride/10;
  if(nt==0)nt=1;
  
  double sum=0;
  #pragma omp parallel num_threads(nt)
  {
    std::vector<Vector> omp_deriv(getPositions().size());
    Tensor omp_virial;

    //~ if (doneigh && !doLowComm)
    if (doneigh)
    {
      if(invalidateList)
        nl->update(getPositions());
      
      // Loop over all atoms
      #pragma omp for reduction(+:fin_intensity) nowait
      for(unsigned int i=0;i<nl->getNumberOfLocalAtoms();++i)
      {
         std::vector<unsigned> neighbors;
         unsigned i0=nl->getIndexOfLocalAtom(i);
         neighbors=nl->getNeighbors(i0);
         Vector position_i0=getPosition(i0);
         // Loop over neighbors
         for(unsigned int j=0;j!=neighbors.size();++j)
         {
            unsigned i1=neighbors[j];
            if(getAbsoluteIndex(i0)==getAbsoluteIndex(i1)) continue;
            
            Vector distance=pbcDistance(position_i0,getPosition(i1));
            if(fabs(distance[axis_id])>axis_maxr) continue;
            
            sum+=1.0;
            distance[axis_id]=0;
            if(nt>1)
              calc_value(distance,fin_intensity,omp_deriv[i0],omp_deriv[i1],omp_virial);
            else
              calc_value(distance,fin_intensity,deriv[i0],deriv[i1],virial);
         }
      }
    }
    //~ else if (doneigh && doLowComm)
    //~ {
      //~ if(invalidateList)
        //~ nl->update(getPositions());
      
      //~ // Loop over all atoms
      //~ #pragma omp for reduction(+:fin_intensity) nowait
      //~ for(unsigned int i=0;i<nl->getNumberOfLocalAtoms();++i)
      //~ {
         //~ std::vector<unsigned> neighbors;
         //~ unsigned i0=nl->getIndexOfLocalAtom(i);
         //~ neighbors=nl->getNeighbors(i0);
         //~ Vector position_i0=getPosition(i0);
         //~ // Loop over neighbors
         //~ for(unsigned int j=0;j!=neighbors.size();++j)
         //~ {
            //~ unsigned i1=neighbors[j];
            //~ if(getAbsoluteIndex(i0)==getAbsoluteIndex(i1)) continue;
            //~ Vector distance=pbcDistance(position_i0,getPosition(i1));
            //~ if(fabs(distance[axis_id])>axis_maxr) continue;
            //~ sum+=1.0;
            //~ distance[axis_id]=0;
            //~ if(nt>1)
              //~ calc_value(distance,fin_intensity,omp_deriv[i0],omp_virial);
            //~ else
              //~ calc_value(distance,fin_intensity,deriv[i0],virial);
         //~ }
      //~ }
    //~ }
    else
    {
	  #pragma omp for reduction(+:fin_intensity) nowait
      for(unsigned int i=rank;i<getNumberOfAtoms();i+=stride)
      {
        for(unsigned int j=0;j<i;++j)
        {
            Vector distance=pbcDistance(getPosition(i),getPosition(j));
            if(fabs(distance[axis_id])>axis_maxr) continue;
            sum+=1.0;
            distance[axis_id]=0;
            if(nt>1)
              calc_value(distance,fin_intensity,omp_deriv[i],omp_deriv[j],omp_virial);
            else
              calc_value(distance,fin_intensity,deriv[i],deriv[j],virial);
        }
      }
    }
    
    comm.Barrier();

    #pragma omp critical
    if(nt>1) {
      for(unsigned i=0; i<getPositions().size(); i++) deriv[i]+=omp_deriv[i];
      virial+=omp_virial;
    }
    
  }

  if(!serial) {
	//~ if(!doLowComm)
	//~ {
       comm.Sum(fin_intensity);
       comm.Sum(sum);
       if(!deriv.empty()) comm.Sum(&deriv[0][0],3*deriv.size());
       comm.Sum(virial);
    //~ }
  }
  
  fin_intensity+=fij2;

  for(unsigned i=0; i<deriv.size(); ++i) setAtomsDerivatives(i,deriv[i]);
  setValue           (fin_intensity);
  setBoxDerivatives  (virial);

}

inline void DebyeScatter2D::calc_value(const Vector& distance,double& value,Vector& deriv0,Vector& deriv1,Tensor& virial)
{
  double dfunc=0.;
  double dis_mod=std::sqrt(distance[map_id1]*distance[map_id1]+distance[map_id2]*distance[map_id2]);
  double val = pairing(dis_mod, dfunc);
  
  //~ val*=2;
  //~ dfunc*=2;

  Vector norm_dis(distance/dis_mod);
  Vector dd(dfunc*norm_dis);
  Tensor vv(dd,distance);
    
  if(nl_full_list)
  {
    value+=val;
    deriv0-=dd;
    deriv1+=dd;
    virial-=vv;
  }
  else
  {
    value+=val*2;
    deriv0-=dd*2;
    deriv1+=dd*2;
    virial-=vv*2;
  }
}

//~ inline void DebyeScatter2D::calc_value(const Vector& distance,double& value,Vector& deriv,Tensor& virial)
//~ {
  //~ double dfunc=0.;
  //~ double dis_mod=std::sqrt(distance[map_id1]*distance[map_id1]+distance[map_id2]*distance[map_id2]);
  //~ value += pairing(dis_mod, dfunc);

  //~ Vector norm_dis(distance/dis_mod);
  //~ Vector dd(dfunc*norm_dis);
  //~ Tensor vv(dd,distance);
    
  //~ deriv-=2*dd;
  //~ virial-=vv;
//~ }

inline double DebyeScatter2D::pairing(double distance,double&dfunc) const
{
  double r=distance;
  if(r>maxr)
  {
	  dfunc=0;
	  return 0;
  }
  else if(r>0)
  {
	if(use_grid)
	{
      unsigned q_id=floor(r/deltad);
      dfunc=prederivsf[q_id];
      return preintensity[q_id];
	}
	else
	{
      double qr = q*r;
      double v_ij  = std::tr1::cyl_bessel_j(0,qr);
      double dv_ij = -1.0*q*std::tr1::cyl_bessel_j(1,qr);
      
      double intensity=1;
      double devrf=0;
      //~ if(usew)
      //~ {
        //~ double w_ij=1;
        //~ double dw_ij=0;
      	//~ double pi_rc= pi*r/maxr;
      	//~ w_ij  = std::tr1::cyl_bessel_j(0,pi_rc);
      	//~ dw_ij = (pi/maxr)*std::tr1::cyl_bessel_j(1,pi_rc);
      	//~ intensity =  = 2.0 * fij2 * v_ij * w_ij / qnum;
        //~ devrf= 2.0 * fij2 * (dv_ij * w_ij + v_ij * dw_ij) / qnum;
      //~ }
      //~ else
      //~ {
        intensity = fij2 * v_ij / qnum;
        devrf = fij2 * dv_ij / qnum;;
	  //~ }
      
      dfunc = devrf;
      return intensity;
	}
  }
  else
  {
	  dfunc=0;
	  return fij2;
  }
}

}

}
