/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#include "ActionRegister.h"
#include "tools/NeighborList.h"
#include "tools/Communicator.h"
#include "tools/Tools.h"

#include <time.h>
#include <string>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR Structure Funtion S(Q)
/*
Calculate the global pair entropy using the expression:
\f[
s=-2\pi\rho k_B \int\limits_0^{r_{\mathrm{max}}} \left [ g(r) \ln g(r) - g(r) + 1 \right ] r^2 dr .
\f]
where \f$ g(r) $\f is the pair distribution function and \f$ r_{\mathrm{max}} $\f is a cutoff in the integration (MAXR).
For the integration the interval from 0 to  \f$ r_{\mathrm{max}} $\f is partitioned in NHIST equal intervals. 
To make the calculation of \f$ g(r) $\f differentiable, the following function is used:
\f[
g(r) = \frac{1}{4 \pi \rho r^2} \sum\limits_{j} \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(r-r_{ij})^2/(2\sigma^2)} ,
\f]
where \f$ \rho $\f is the density and \f$ sigma $\f is a broadening parameter (SIGMA).  
\par Example)
The following input tells plumed to calculate the pair entropy of atoms 1-250 with themselves.
\verbatim
PAIRENTROPY ...
 LABEL=s2
 GROUPA=1-250
 MAXR=0.65
 SIGMA=0.025
 NHIST=100
 NLIST
 NL_CUTOFF=0.75
 NL_STRIDE=10
... PAIRENTROPY
\endverbatim
*/
//+ENDPLUMEDOC

class IntensityO : public Colvar {
  bool pbc;
  bool serial;
  NeighborList *nl;
  bool invalidateList;
  bool firsttime;
  std::vector<double> theta;
  double q;
  double maxr;
  double   lamda;
  unsigned qhist;
  double rcut2;

  unsigned qbin;
  std::vector<double> preintensity;
  std::vector<double> prederivsf;

  std::vector<double> structureF;
  std::vector<double> finalsf;
  std::vector<double> value;
  std::vector<Value*> valueSF;


public:
  explicit IntensityO(const ActionOptions&);
  ~IntensityO();
// active methods:
  virtual void calculate();
  virtual void prepare();
  static void registerKeywords( Keywords& keys );
};



PLUMED_REGISTER_ACTION(IntensityO,"INTENSITYO")

void IntensityO::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
  keys.add("atoms","GROUPA","First list of atoms");
  keys.add("atoms","GROUPB","Second list of atoms (if empty, N*(N-1)/2 pairs in GROUPA are counted)");
  keys.add("compulsory","LAMDA","0.15406","Wavelength of the incident wave");
  keys.add("compulsory","MAXR","1","Maximum distance for the radial distribution function ");
  keys.add("compulsory","THETA","1","Diffraction angles");
  keys.add("compulsory","NHIST","1","Number of bins in the rdf ");
  keys.add("compulsory","QHIST","1","Number of bins of Q ");
  keys.add("compulsory","QBIN","1","Number of bins of Q ");


  keys.addOutputComponent("SF","default","the calculated structure function");
  ActionWithValue::useCustomisableComponents(keys); //The components in the action will depend on the user. 
}

IntensityO::IntensityO(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false),
invalidateList(true),
firsttime(true)
{


  parseFlag("SERIAL",serial);

  vector<AtomNumber> ga_lista,gb_lista;
  parseAtomList("GROUPA",ga_lista);
  parseAtomList("GROUPB",gb_lista);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

// pair stuff
  bool dopair=false;
  parseFlag("PAIR",dopair);

// neighbor list stuff
  bool doneigh=false;
  double nl_cut=0.0;
  int nl_st=0;
  parseFlag("NLIST",doneigh);
  if(doneigh){
   parse("NL_CUTOFF",nl_cut);
   if(nl_cut<=0.0) error("NL_CUTOFF should be explicitly specified and positive");
   parse("NL_STRIDE",nl_st);
   if(nl_st<=0) error("NL_STRIDE should be explicitly specified and positive");
  }
  parse("QHIST",qhist);
  log.printf("Calculation of structure function of Q from 0 to %u. \n", qhist );

  parse("QBIN",qbin);
  preintensity.assign(qbin,0);
  prederivsf.assign(qbin,0);

  structureF.resize(qhist);
  finalsf.resize(qhist);
  valueSF.resize(qhist);
  value.resize(qhist);
 
  int pos_count=1;
  std::ostringstream oss;
  for (unsigned int k=0;k<qhist;k++)
  {
   pos_count++;
   oss.str("");
   oss<<"SF["<<k<<"]"; 
   addComponentWithDerivatives(oss.str());  componentIsNotPeriodic(oss.str());  valueSF[k]=getPntrToComponent(oss.str());
  }


  if(gb_lista.size()>0){
    if(doneigh)  nl= new NeighborList(ga_lista,gb_lista,dopair,pbc,getPbc(),nl_cut,nl_st);
    else         nl= new NeighborList(ga_lista,gb_lista,dopair,pbc,getPbc());
  } else {
    if(doneigh)  nl= new NeighborList(ga_lista,pbc,getPbc(),nl_cut,nl_st);
    else         nl= new NeighborList(ga_lista,pbc,getPbc());
  }

  requestAtoms(nl->getFullAtomList());

  log.printf("  between two groups of %u and %u atoms\n",static_cast<unsigned>(ga_lista.size()),static_cast<unsigned>(gb_lista.size()));
  log.printf("  first group:\n");
  for(unsigned int i=0;i<ga_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n  second group:\n");
  for(unsigned int i=0;i<gb_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", gb_lista[i].serial());
  }
  log.printf("  \n");
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  if(dopair) log.printf("  with PAIR option\n");
  if(doneigh){
   log.printf("  using neighbor lists with\n");
   log.printf("  update every %d steps and cutoff %f\n",nl_st,nl_cut);
  }


  parse("LAMDA",lamda);

  parseVector("THETA",theta);
//  log.printf("targeted diffraction theta %f,%f,%f. \n", theta[0],theta[1],theta[2] );
  parse("MAXR",maxr);
  log.printf("Integration in the interval from 0. to %f nm. \n", maxr );
  
  if(doneigh){
    if(nl_cut<maxr) error("NL_CUTOFF should be larger than MAXR + 3*SIGMA");
  }

  checkRead();



  // Define heavily used expressions
}

IntensityO::~IntensityO(){
  delete nl;
}

void IntensityO::prepare(){
  if(nl->getStride()>0){
    if(firsttime || (getStep()%nl->getStride()==0)){
      requestAtoms(nl->getFullAtomList());
      invalidateList=true;
      firsttime=false;
    }else{
      requestAtoms(nl->getReducedAtomList());
      invalidateList=false;
      if(getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
    }
    if(getExchangeStep()) firsttime=true;
  }
}

// calculator
void IntensityO::calculate()
{

  // Setup neighbor list and parallelization
  if(nl->getStride()>0 && invalidateList){
    nl->update(getPositions());
  }


  Matrix<Vector> sfPrime(qhist,getNumberOfAtoms());
  vector<Tensor> sfVirial(qhist);

  // Loop over neighbors
  const unsigned nn=nl->size();
//   log.printf(" Distance nn %d \n",nn);

  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();

    double deltad=maxr/qbin;

  rcut2=maxr*maxr; 
    double d2;
    Vector distance;
    double distanceModulo;
    Vector distance_versor;
    double fij_O;
  for(unsigned int m=0;m<qhist;++m){
    value[m]=0.0;
    q=4*pi*std::sin(theta[m]*pi/180)/lamda;

    double q_A=4*pi*std::sin(theta[m]*pi/180)/(lamda*10);
         fij_O= 3.0485*exp(-13.2771*(q_A/(4*pi))*(q_A/(4*pi))) + 
                2.2868*exp(-5.7011*(q_A/(4*pi))*(q_A/(4*pi))) + 
                1.5463*exp(-0.3239*(q_A/(4*pi))*(q_A/(4*pi))) + 
                0.8670*exp(-32.9089*(q_A/(4*pi))*(q_A/(4*pi))) + 0.2508;
  //~ log.printf("  with X-ray wave length %f, angle %f and Q value %f.\n",lamda,theta[m],q);
  //~ log.printf("  with atomic form factors %f.\n",fij_O);
  for(unsigned int j=0;j<qbin;j+=1){
  double predistanceM= deltad*(j+0.5);
    preintensity[j]=2*std::sin(q*predistanceM)/(getNumberOfAtoms()*q*predistanceM);
    prederivsf[j]=-2/q*((q*std::cos(q*predistanceM)/predistanceM)-(std::sin(q*predistanceM)/(predistanceM*predistanceM)))/(getNumberOfAtoms());
  	//~ log.printf("    %d. with distance %f, QR %f, value %e and Derivative %e.\n", int(j),predistanceM,predistanceM*q,preintensity[j],prederivsf[j]);

   }

//        fij_O= 3.0485*exp(-13.2771*(q/(4*pi))*(q/(4*pi)))+2.2868*exp(-5.7011*(q/(4*pi))*(q/(4*pi)))+1.5463*exp(-0.3239*(q/(4*pi))*(q/(4*pi)))+0.867*exp(-32.9089*(q/(4*pi))*(q/(4*pi)))+0.2508;

//   for(unsigned int i=0;i<nn;++i) {
   for(unsigned int i=rank;i<nn;i+=stride) {
    unsigned i0=nl->getClosePair(i).first;
    unsigned i1=nl->getClosePair(i).second;
    if(getAbsoluteIndex(i0)==getAbsoluteIndex(i1)) continue;
    if(pbc){
     distance=pbcDistance(getPosition(i0),getPosition(i1));
    } else {
     distance=delta(getPosition(i0),getPosition(i1));
    }
      d2=distance[0]*distance[0]+distance[1]*distance[1]+distance[2]*distance[2]; 
     distanceModulo=std::sqrt(d2);
//     if(distanceModulo<maxr){
     distance_versor = distance / distanceModulo;
//    log.printf(" Distance %f\n",std::sin(q*distanceModulo));
//        double lf= ((1+std::cos(2*theta[m]*pi/180)*std::cos(2*theta[m]*pi/180))/(std::sin(theta[m]*pi/180)*std::sin(theta[m]*pi/180)*std::cos(theta[m]*pi/180)));
//        structureF[m]+=2*std::sin(q*distanceModulo)/(getNumberOfAtoms()*q*distanceModulo)*std::sin(pi*distanceModulo/maxr)*maxr/(pi*distanceModulo);
//          value[m]+=2*std::sin(q*distanceModulo)/(getNumberOfAtoms()*q*distanceModulo)*std::sin(pi*distanceModulo/maxr)*maxr/(pi*distanceModulo);

       Vector derivsf;
     if(distanceModulo<maxr){
     int q_num=distanceModulo/deltad;
      value[m]+=preintensity[q_num];
      derivsf=prederivsf[q_num]*distance_versor;
   }
       sfPrime[m][i0]+= derivsf;
       sfPrime[m][i1]+= -derivsf;
       Tensor vv(derivsf, distance);
        sfVirial[m] += vv ;
//   }
   }
       comm.Sum(value);
       comm.Sum(sfPrime);
       comm.Sum(sfVirial);

         structureF[m]=fij_O*fij_O*(1+value[m]);
//  log.printf(" fij_O %f \n",fij_O);

         for(unsigned i=0;i<getNumberOfAtoms();++i) setAtomsDerivatives(valueSF[m],i,fij_O*fij_O*sfPrime[m][i]);
         valueSF[m]->set(structureF[m]);
         setBoxDerivatives  (valueSF[m],fij_O*fij_O*sfVirial[m]);

  }


}

}
}
