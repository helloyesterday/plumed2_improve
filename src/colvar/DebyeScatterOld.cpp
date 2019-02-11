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
#include "tools/NeighborList.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "ActionRegister.h"

#include <string>
#include <memory>

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

class DebyeScatterOld : public Colvar {
  bool pbc;
  bool serial;
  std::unique_ptr<NeighborList> nl;
  bool invalidateList;
  bool firsttime;
  bool usew;
  bool use_grid;
  
  double theta;
  double maxr;
  double lambda;
  double deltad;
  double q;
  double fij;
  double fij2;

  unsigned qbin;
  unsigned qnum;
  std::vector<double> preintensity;
  std::vector<double> prederivsf;

public:
  explicit DebyeScatterOld(const ActionOptions&);
  ~DebyeScatterOld();
  static void registerKeywords( Keywords& keys );
// active methods:
  virtual void calculate();
  virtual void prepare();
  virtual double pairing(double distance,double&dfunc,unsigned i,unsigned j)const;
};

PLUMED_REGISTER_ACTION(DebyeScatterOld,"DEBYE_SCATTER_OLD")

void DebyeScatterOld::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.addFlag("USE_WINDOW",false,"Use a window function to revise the finite box");
  keys.addFlag("GRID_INTENSITY",false,"Use a grid to caluculate the intensity by distance");
  keys.add("atoms","ATOMS","First list of atoms");
  keys.add("compulsory","THETA","Diffraction angles");
  keys.add("compulsory","LAMBDA","0.15406","Wavelength of the incident wave");
  keys.add("compulsory","MAXR","3.0","Maximum distance for the radial distribution function ");
  keys.add("optional","QBIN","Number of grid bins of intensity when use GRID_INTENSITY");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
}

DebyeScatterOld::DebyeScatterOld(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  invalidateList(true),
  firsttime(true)
{

  parseFlag("SERIAL",serial);

  vector<AtomNumber> ga_lista;
  parseAtomList("ATOMS",ga_lista);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
    
  parse("LAMBDA",lambda);
  parse("THETA",theta);
  parse("MAXR",maxr);
  parse("QBIN",qbin);

// neighbor list stuff
  bool doneigh=false;
  double nl_cut=0.0;
  int nl_st=0;
  parseFlag("NLIST",doneigh);
  if(doneigh) {
    parse("NL_CUTOFF",nl_cut);
    if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
    if(nl_cut<maxr) error("NL_CUTOFF should not be smaller than MAXR");
    parse("NL_STRIDE",nl_st);
    if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");
  }
  
  parseFlag("USE_WINDOW",usew);
  parseFlag("GRID_INTENSITY",use_grid);

  addValueWithDerivatives(); setNotPeriodic();

  if(doneigh)  nl.reset( new NeighborList(ga_lista,pbc,getPbc(),nl_cut,nl_st) );
  else         nl.reset( new NeighborList(ga_lista,pbc,getPbc()) );

  requestAtoms(nl->getFullAtomList());
  
  qnum=ga_lista.size();

  checkRead();
  
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
      double r2 = r*r;
      double sinqr = std::sin(qr);
      double cosqr = std::cos(qr);
      
      double v_ij = sinqr/qr;
      double dv_ij = (qr*cosqr-sinqr)/(q*r2);
      
      double intensity=1;
      double devrf=0;
      if(usew)
      {
        double pi_rc= pi*r/maxr;
        double sin_pi_rc = std::sin(pi_rc);
        double cos_pi_rc = std::cos(pi_rc);
        double w_ij = sin_pi_rc/pi_rc;
        double dw_ij = (pi_rc*cos_pi_rc-sin_pi_rc)/(pi*r2/maxr);
      
        intensity = 2.0 * fij2 * v_ij * w_ij / qnum;
        devrf = 2.0 * fij2 * (dv_ij * w_ij + v_ij * dw_ij) / qnum;
      }
      else
      {
        intensity = 2.0 * fij2 * v_ij/ qnum;
        devrf = 2.0 * fij2 * dv_ij / qnum;
	  }
		preintensity[j] = intensity;
		prederivsf[j] = devrf;
	}
  }

  log.printf("  ATOMS:\n");
  for(unsigned int i=0; i<ga_lista.size(); ++i) {
    if ( (i+1) % 25 == 0 ) log.printf("  \n");
    log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n");
  log.printf("  with CV atom number %d\n",int(qnum),getNumberOfAtoms());
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  if(usew) log.printf("  using window function to revise finite box\n");
  if(doneigh) {
    log.printf("  using neighbor lists with\n");
    log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);
  }

  log.printf("  with X-ray wave length %f, angle %f and Q value %f.\n",lambda,theta,q);
  log.printf("  with atomic form factors %f.\n",fij);
  log.printf("  with upper limied %f.\n", maxr);
  if(use_grid) log.printf("  with integrated number %d and interval %f.\n", qbin,deltad);
  //~ for(unsigned j=0;j!=qbin;++j)
  //~ {
	//~ double r = deltad*(j+0.5);
    //~ double qr = q*r;
	//~ log.printf("    %d. with distance %f, QR %f, value %e and Derivative %e.\n", int(j),r,qr,preintensity[j],prederivsf[j]);
  //~ }
}

DebyeScatterOld::~DebyeScatterOld()
{
}

void DebyeScatterOld::prepare() {
  if(nl->getStride()>0) {
    if(firsttime || (getStep()%nl->getStride()==0)) {
      requestAtoms(nl->getFullAtomList());
      invalidateList=true;
      firsttime=false;
    } else {
      requestAtoms(nl->getReducedAtomList());
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}

// calculator
void DebyeScatterOld::calculate()
{

  double fin_intensity=0.;
  Tensor virial;
  vector<Vector> deriv(getNumberOfAtoms());
// deriv.resize(getPositions().size());

  if(nl->getStride()>0 && invalidateList) {
    nl->update(getPositions());
  }

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

  const unsigned nn=nl->size();

  if(nt*stride*10>nn) nt=nn/stride/10;
  if(nt==0)nt=1;
  
  double sum=0;

  #pragma omp parallel num_threads(nt)
  {
    std::vector<Vector> omp_deriv(getPositions().size());
    Tensor omp_virial;

    #pragma omp for reduction(+:fin_intensity) nowait
    for(unsigned int i=rank; i<nn; i+=stride) {

      Vector distance;
      unsigned i0=nl->getClosePair(i).first;
      unsigned i1=nl->getClosePair(i).second;

      if(getAbsoluteIndex(i0)==getAbsoluteIndex(i1))
		  continue;

      if(pbc) {
        distance=pbcDistance(getPosition(i0),getPosition(i1));
      } else {
        distance=delta(getPosition(i0),getPosition(i1));
      }
      
      sum+=1.0;

      double dfunc=0.;
      double dis_mod=distance.modulo();
      fin_intensity += pairing(dis_mod, dfunc,i0,i1);

      Vector norm_dis(distance/dis_mod);
      Vector dd(dfunc*norm_dis);
      Tensor vv(dd,distance);
      if(nt>1) {
        omp_deriv[i0]-=dd;
        omp_deriv[i1]+=dd;
        omp_virial-=vv;
      } else {
        deriv[i0]-=dd;
        deriv[i1]+=dd;
        virial-=vv;
      }

    }
    #pragma omp critical
    if(nt>1) {
      for(unsigned i=0; i<getPositions().size(); i++) deriv[i]+=omp_deriv[i];
      virial+=omp_virial;
    }
  }

  if(!serial) {
    comm.Sum(fin_intensity);
    comm.Sum(sum);
    if(!deriv.empty()) comm.Sum(&deriv[0][0],3*deriv.size());
    comm.Sum(virial);
  }
  
  fin_intensity+=fij2;

  for(unsigned i=0; i<deriv.size(); ++i) setAtomsDerivatives(i,deriv[i]);
  setValue           (fin_intensity);
  setBoxDerivatives  (virial);

}

double DebyeScatterOld::pairing(double distance,double&dfunc,unsigned i,unsigned j)const {
  (void) i; // avoid warnings
  (void) j; // avoid warnings
  
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
      double r2 = r*r;
      double sinqr = std::sin(qr);
      double cosqr = std::cos(qr);
      
      double v_ij = sinqr/qr;
      double dv_ij = (qr*cosqr-sinqr)/(q*r2);
      
      double intensity=0;
      if(usew)
      {
        double pi_rc= pi*r/maxr;
        double sin_pi_rc = std::sin(pi_rc);
        double cos_pi_rc = std::cos(pi_rc);
        double w_ij = sin_pi_rc/pi_rc;
        double dw_ij = (pi_rc*cos_pi_rc-sin_pi_rc)/(pi*r2/maxr);
        
        intensity = 2.0 * fij2 * v_ij * w_ij / qnum;
        dfunc = 2.0 * fij2 * (dv_ij * w_ij + v_ij * dw_ij) / qnum;
      }
      else
      {
		intensity = 2.0 * fij2 * v_ij/ qnum;
        dfunc = 2.0 * fij2 * dv_ij / qnum;
	  }
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
