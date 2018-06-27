/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "tools/Torsion.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC COLVAR ALPHABETA2
/*
Measures a distance including pbc between the instantaneous values of a set of torsional angles and set of reference values.

This colvar calculates the following quantity.

\f[
s = \sum_i w_i \left[ \frac{\cos( \phi_i - \phi_i^{\textrm{Ref_B}} ) - \cos( \phi_i - \phi_i^{\textrm{Ref_A}} )}
{2 - \cos( \phi_i - \phi_i^{\textrm{Ref_A}} ) - \cos( \phi_i - \phi_i^{\textrm{Ref_B}} )} \right]   
\f]

where the \f$\phi_i\f$ values are the instantaneous values for the \ref TORSION angles of interest.
The \f$\phi_i^{\textrm{Ref}}\f$ values are the user-specified reference values for the torsional angles.

\par Examples

The following provides an example of the input for an alpha beta similarity.

\verbatim
ALPHABETA ...
ATOMS1=168,170,172,188 REFA1=3.14 REFB1=1.22
ATOMS2=170,172,188,190 REFA2=3.14 REFB2=1.22
ATOMS3=188,190,192,230 REFA3=3.14 REFB3=1.22
LABEL=ab
... ALPHABETA
PRINT ARG=ab FILE=colvar STRIDE=10
\endverbatim

Because all the reference values are the same we can calculate the same quantity using

\verbatim
ALPHABETA ...
ATOMS1=168,170,172,188 REFA1=3.14 REFB1=1.22
ATOMS2=170,172,188,190 REFA2=3.14 REFB2=1.22
ATOMS3=188,190,192,230 REFA3=3.14 REFB3=1.22
LABEL=ab
... ALPHABETA
PRINT ARG=ab FILE=colvar STRIDE=10
\endverbatim

Writing out the atoms involved in all the torsions in this way can be rather tedious. Thankfully if you are working with protein you
can avoid this by using the \ref MOLINFO command.  PLUMED uses the pdb file that you provide to this command to learn 
about the topology of the protein molecule.  This means that you can specify torsion angles using the following syntax:

\verbatim
MOLINFO MOLTYPE=protein STRUCTURE=myprotein.pdb
ALPHABETA ...
ATOMS1=@phi-3 REFA1=3.14 REFB1=1.22
ATOMS2=@psi-3 REFA2=3.14 REFB2=1.22
ATOMS3=@phi-4 REFA3=3.14 REFB3=1.22
LABEL=ab
... ALPHABETA 
PRINT ARG=ab FILE=colvar STRIDE=10
\endverbatim

Here, \@phi-3 tells plumed that you would like to calculate the \f$\phi\f$ angle in the third residue of the protein.  
Similarly \@psi-4 tells plumed that you want to calculate the \f$\psi\f$ angle of the 4th residue of the protein.


*/
//+ENDPLUMEDOC

class Torsional_Orderness : public MultiColvarBase {
private:
  std::vector<double> target;
  std::vector<double> var;
  std::vector<double> var2;
  std::vector<double> cos0;
  std::vector<double> sin0;
  std::vector<double> weights;
  double wnorm;
  bool need_normal;
  bool need_normal2;
public:
  static void registerKeywords( Keywords& keys );
  explicit Torsional_Orderness(const ActionOptions&);
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
  bool isPeriodic() { return false; }
};

PLUMED_REGISTER_ACTION(Torsional_Orderness,"TORSIONAL_ORDERNESS")

void Torsional_Orderness::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.add("numbered","ATOMS","the atoms involved in each of the alpha-beta variables you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one alpha-beta values will be "
           "calculated for each ATOM keyword you specify (all ATOM keywords should "
           "specify the indices of four atoms).  The eventual number of quantities calculated by this "
           "action will depend on what functions of the distribution you choose to calculate.");
  keys.reset_style("ATOMS","atoms");
  keys.add("numbered","REFERENCE","the reference values for each of the torsional angles.  If you use a single REFERENCE value the "
           "same reference value is used for all torsions");
  keys.reset_style("REFERENCE","compulsory");
  keys.add("numbered","SIGMA","the variance of the normal distribution, by default is 0.5");
  keys.add("numbered","WEIGHT","A weight value for a given contact, by default is 1.0 "
                               "You can either specify a global weight value using WEIGHT or one "
                               "weight value for each contact."); 
  keys.addFlag("NORMALIZE",false,"to normalize the weights as the summation of the values equal to 1");
  keys.addFlag("NORMALIZE2",false,"to normalize the weights as the summation of the squre of values equal to 1");
}

Torsional_Orderness::Torsional_Orderness(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
	
  // Read in the atoms
  std::vector<AtomNumber> all_atoms;
  readAtomsLikeKeyword( "ATOMS", 4, all_atoms );
  setupMultiColvarBase( all_atoms );
  // Resize target
  target.resize( getFullNumberOfTasks() );
  var.resize( getFullNumberOfTasks() ,0.5);
  var2.resize( getFullNumberOfTasks() ,0.25);
  cos0.resize( getFullNumberOfTasks() );
  sin0.resize( getFullNumberOfTasks() );
  //~ csin12.resize( getFullNumberOfTasks() );
  // Resize weights
  weights.assign( getFullNumberOfTasks() ,1);
  // Setup central atom indices
  std::vector<bool> catom_ind(4, false);
  catom_ind[1]=catom_ind[2]=true;
  setAtomsForCentralAtom( catom_ind );

  // Read in reference values
  for(unsigned i=0;i<target.size();++i){
     if( !parseNumbered( "REFERENCE", i+1, target[i] ) )
     {
		 string eno;
		 Tools::convert(int(i+1),eno);
		 error("Can't find REFERENCE"+eno+"!");
	 }
	 cos0[i]=cos(target[i]);
	 sin0[i]=sin(target[i]);
  }
  
  // Read in sigma
  unsigned nvar=0;
  for(unsigned i=0;i<var.size();++i){
     if( !parseNumbered( "SIGMA", i+1, var[i] ) ) break;
     var2[i]=var[i]*var[i];
     nvar++;
  }
  if(nvar!=0 && nvar!=var.size() ){
      error("found wrong number of SIGMA values");
  }
  
  
  // Read in weights
  unsigned nweights=0;
  for(unsigned i=0;i<weights.size();++i){
     if( !parseNumbered( "WEIGHT", i+1, weights[i] ) ) break;
     nweights++;
  }
  if(nweights!=0 && nweights!=weights.size() ){
      error("found wrong number of WEIGHT values");
  }
  
  // Normalization
  parseFlag("NORMALIZE",need_normal);
  parseFlag("NORMALIZE2",need_normal2);
  if(need_normal&&need_normal2)
      error("FLAG \"NORMALIZE\" and \"NORMALIZE2\" can be only accepted one!");
  if(need_normal)
  {
	  wnorm=0;
      for(unsigned i=0;i!=weights.size();++i)
         wnorm+=weights[i];
      for(unsigned i=0;i!=weights.size();++i)
         weights[i]/=wnorm;
  }
  else if(need_normal2)
  {
	  wnorm=norm(weights);
      for(unsigned i=0;i!=weights.size();++i)
         weights[i]/=wnorm;
  }

  // And setup the ActionWithVessel
  if( getNumberOfVessels()==0 ) {
    std::string fake_input;
    addVessel( "SUM", fake_input, -1 );  // -1 here means that this value will be named getLabel()
    readVesselKeywords();  // This makes sure resizing is done
  }
  // And check everything has been read in correctly
  checkRead();
}

double Torsional_Orderness::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  const Vector d0=getSeparation(myatoms.getPosition(1),myatoms.getPosition(0));
  const Vector d1=getSeparation(myatoms.getPosition(2),myatoms.getPosition(1));
  const Vector d2=getSeparation(myatoms.getPosition(3),myatoms.getPosition(2));

  Vector dd0,dd1,dd2;
  PLMD::Torsion t;
  const double value  = t.compute(d0,d1,d2,dd0,dd1,dd2);
  const double vcos   = cos(value);
  const double vcos0  = vcos - cos0[tindex];
  const double vsin   = sin(value);
  const double vsin0  = vsin - sin0[tindex];
  const double cvalue = exp(-(vcos0*vcos0+vsin0*vsin0)/var2[tindex])*weights[tindex];
  const double svalue = 2.0*(sin0[tindex]*vcos-cos0[tindex]*vsin)*cvalue/var2[tindex];

  dd0 *= svalue;
  dd1 *= svalue;
  dd2 *= svalue;

  addAtomDerivatives(1, 0, dd0, myatoms);
  addAtomDerivatives(1, 1, dd1-dd0, myatoms);
  addAtomDerivatives(1, 2, dd2-dd1, myatoms);
  addAtomDerivatives(1, 3, -dd2, myatoms);

  myatoms.addBoxDerivatives(1, -(extProduct(d0,dd0)+extProduct(d1,dd1)+extProduct(d2,dd2)));

  return cvalue;
}


}
}
