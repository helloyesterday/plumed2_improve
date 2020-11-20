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


#include <tr1/cmath>
#include <iostream>


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

class TDIntensityO : public Colvar {
  bool pbc;
  bool serial;
  NeighborList *nl;
  bool invalidateList;
  bool firsttime;
  std::vector<double> theta;
  double   thetamin,thetamax;
  double q;
  double zlow;
  double zup;
  double   lamda;
  unsigned qhist;
  double rcut2;
  unsigned qbin;


  double maxr;
  std::vector<double> preintensity;
  std::vector<double> prederivsf;


  std::vector<double> structureF;
  std::vector<double> structureF2;
  std::vector<double> finalsf;
  std::vector<double> value;
  std::vector<Value*> valueSF;


public:
  explicit TDIntensityO(const ActionOptions&);
  ~TDIntensityO();
// active methods:
  virtual void calculate();
  virtual void prepare();
  static void registerKeywords( Keywords& keys );
};



PLUMED_REGISTER_ACTION(TDIntensityO,"TDINTENSITYO")

void TDIntensityO::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.addFlag("NLIST",false,"Use a neighbour list to speed up the calculation");
  keys.add("optional","NL_CUTOFF","The cutoff for the neighbour list");
  keys.add("optional","NL_STRIDE","The frequency with which we are updating the atoms in the neighbour list");
  keys.add("atoms","GROUPA","First list of atoms");
  keys.add("atoms","GROUPB","Second list of atoms (if empty, N*(N-1)/2 pairs in GROUPA are counted)");
  keys.add("compulsory","ZMAXL","0","Wavelength of the incident wave");
  keys.add("compulsory","ZMAXU","1","Wavelength of the incident wave");
  keys.add("compulsory","LAMDA","0.15406","Wavelength of the incident wave");
  keys.add("compulsory","THETA","1","Diffraction angles");
  keys.add("compulsory","NHIST","1","Number of bins in the rdf ");
  keys.add("compulsory","QHIST","1","Number of bins of Q ");
  keys.add("compulsory","MAXR","1","Maximum distance for the radial distribution function ");
  keys.add("compulsory","QBIN","1","Number of bins of Q ");


  keys.addOutputComponent("SF","default","the calculated structure function");
  ActionWithValue::useCustomisableComponents(keys); //The components in the action will depend on the user. 
}

TDIntensityO::TDIntensityO(const ActionOptions&ao):
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

  parse("MAXR",maxr);


  parse("QBIN",qbin);
  preintensity.assign(qbin,0);
  prederivsf.assign(qbin,0);


  structureF.resize(qhist);
  structureF2.resize(qhist);
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
  parse("ZMAXL",zlow);
  parse("ZMAXU",zup);

  parseVector("THETA",theta);


  checkRead();
}

TDIntensityO::~TDIntensityO(){
  delete nl;
}

void TDIntensityO::prepare(){
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
void TDIntensityO::calculate()
{

  // Setup neighbor list and parallelization
  if(nl->getStride()>0 && invalidateList){
    nl->update(getPositions());
  }


  Matrix<Vector> sfPrime(qhist,getNumberOfAtoms());
  vector<Tensor> sfVirial(qhist);


//  double thetastep=(thetamax-thetamin)/(qhist-1);
//     theta.resize(qhist);


  // Loop over neighbors
  const unsigned nn=nl->size();
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();

    int count1;
    int N_atom;
    double d2;
    Vector distance;
    double distanceModulo;
    Vector distance_versor;
    Vector derivsf;

    double deltad=maxr/qbin;
    double fij_O;
    N_atom=0;

  for(unsigned int m=0;m<qhist;++m){
    value[m]=0.0;
//    theta[m]=thetamin+thetastep*m;
    q=4*pi*std::sin(theta[m]*pi/180)/lamda;
    count1=0;




  for(unsigned int j=0;j<qbin;j+=1){
  double predistanceM= deltad*(j+0.5);
    preintensity[j]=std::tr1::cyl_bessel_j(0,q*predistanceM);
    prederivsf[j]=q*std::tr1::cyl_bessel_j(1,q*predistanceM);
   }


    double q_A=4*pi*std::sin(theta[m]*pi/180)/(lamda*10);
         fij_O= 3.0485*exp(-13.2771*(q_A/(4*pi))*(q_A/(4*pi)))+2.2868*exp(-5.7011*(q_A/(4*pi))*(q_A/(4*pi)))+1.5463*exp(-0.3239*(q_A/(4*pi))*(q_A/(4*pi)))+0.867*exp(-32.9089*(q_A/(4*pi))*(q_A/(4*pi)))+0.2508;

   for(unsigned int i=rank;i<nn;i+=stride) {
    unsigned i0=nl->getClosePair(i).first;
    unsigned i1=nl->getClosePair(i).second;
    Vector pos0=getPosition(i0);
    Vector pos1=getPosition(i1);

    if(pos0[2]>zlow && pos0[2]<zup && pos1[2]>zlow && pos1[2]<zup) {
    count1=count1+1;
    if(getAbsoluteIndex(i0)==getAbsoluteIndex(i1)) continue;
    if(pbc){
     distance=pbcDistance(getPosition(i0),getPosition(i1));
    } else {
     distance=delta(getPosition(i0),getPosition(i1));
    }
      d2=distance[0]*distance[0]+distance[1]*distance[1]; 
      distanceModulo=std::sqrt(d2);
      distance[2]=0.0;
      distance_versor = distance / distanceModulo;
      distance_versor[2]=0.0;

     if(distanceModulo<maxr){
     int q_num=distanceModulo/deltad;
      value[m]+=preintensity[q_num];
      derivsf=prederivsf[q_num]*distance_versor;
   }

      sfPrime[m][i0]+= derivsf;
      sfPrime[m][i1]+= -derivsf;
//      derivsf[0]=0;
//      derivsf[1]=0;
//      derivsf[2]=0;
      Tensor vv(derivsf, distance);
      sfVirial[m] += vv ;
   }
   }
       comm.Sum(value);
       comm.Sum(sfPrime);
       comm.Sum(sfVirial);
       comm.Sum(count1);
       N_atom=int(sqrt(2*count1));


         structureF[m]=fij_O*fij_O*(1+value[m])/N_atom;
         for(unsigned i=0;i<getNumberOfAtoms();++i) setAtomsDerivatives(valueSF[m],i,fij_O*fij_O*sfPrime[m][i]/N_atom);
//    log.printf(" 2theta and Intensity %f, %f\n",2*theta[m],structureF[m]);

         valueSF[m]->set(structureF[m]);
         setBoxDerivatives  (valueSF[m],fij_O*fij_O*sfVirial[m]/N_atom);


  }


}

}
}
