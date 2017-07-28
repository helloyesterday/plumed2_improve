/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
#include "MultiColvar.h"
#include "tools/Torsion.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace multicolvar{

//+PLUMEDOC COLVAR ALPHABETA2 
/*
Measures a distance including pbc between the instantaneous values of a set of torsional angles and set of reference values.

This colvar calculates the following quantity.

\f[
s = \frac{1}{2} \sum_i \left[ 1 + \cos( \phi_i - \phi_i^{\textrm{Ref}} ) \right]   
\f]

where the \f$\phi_i\f$ values are the instantaneous values for the \ref TORSION angles of interest.
The \f$\phi_i^{\textrm{Ref}}\f$ values are the user-specified reference values for the torsional angles.

\par Examples

The following provides an example of the input for an alpha beta similarity.

\verbatim
ALPHABETA ...
ATOMS1=168,170,172,188 REFERENCE1=3.14 
ATOMS2=170,172,188,190 REFERENCE2=3.14 
ATOMS3=188,190,192,230 REFERENCE3=3.14
LABEL=ab
... ALPHABETA
PRINT ARG=ab FILE=colvar STRIDE=10
\endverbatim

Because all the reference values are the same we can calculate the same quantity using

\verbatim
ALPHABETA ...
ATOMS1=168,170,172,188 REFERENCE=3.14 
ATOMS2=170,172,188,190 
ATOMS3=188,190,192,230 
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
ATOMS1=@phi-3 REFERENCE=3.14
ATOMS2=@psi-3
ATOMS3=@phi-4
LABEL=ab
... ALPHABETA 
PRINT ARG=ab FILE=colvar STRIDE=10
\endverbatim

Here, \@phi-3 tells plumed that you would like to calculate the \f$\phi\f$ angle in the third residue of the protein.  
Similarly \@psi-4 tells plumed that you want to calculate the \f$\psi\f$ angle of the 4th residue of the protein.


*/
//+ENDPLUMEDOC

class AlphaBeta2 : public MultiColvar {
private:
  std::vector<double> target1;
  std::vector<double> target2;
  std::vector<double> csin12;
public:
  static void registerKeywords( Keywords& keys );
  explicit AlphaBeta2(const ActionOptions&);
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(AlphaBeta2,"ALPHABETA2")

void AlphaBeta2::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  keys.use("ATOMS");
  keys.add("numbered","REFA","the reference values for each of the first torsional angles.");
  keys.add("numbered","REFB","the reference values for each of the second torsional angles.");
}

AlphaBeta2::AlphaBeta2(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  // Read in the atoms
  int natoms=4; std::vector<AtomNumber> all_atoms;
  readAtoms( natoms, all_atoms );
  // Resize target
  target1.resize( getFullNumberOfTasks() );
  target2.resize( getFullNumberOfTasks() );
  csin12.resize( getFullNumberOfTasks() );
  // Setup central atom indices
  std::vector<bool> catom_ind(4, false); 
  catom_ind[1]=catom_ind[2]=true;
  setAtomsForCentralAtom( catom_ind );

  // Read in reference values
  for(unsigned i=0;i<target1.size();++i){
     if( !parseNumbered( "REFA", i+1, target1[i] ) )
     {
		 string eno;
		 Tools::convert(int(i+1),eno);
		 error("Can't find REFA"+eno+"!");
	 }
     if( !parseNumbered( "REFB", i+1, target2[i] ) )
     {
		 string eno;
		 Tools::convert(int(i+1),eno);
		 error("Can't find REFB"+eno+"!");
	 }
	 csin12[i]=sin((target1[i]-target2[i])/2.0);
  }

  // And setup the ActionWithVessel
  if( getNumberOfVessels()==0 ){
     std::string fake_input;
     addVessel( "SUM", fake_input, -1 );  // -1 here means that this value will be named getLabel()
     readVesselKeywords();  // This makes sure resizing is done
  }

  // And check everything has been read in correctly
  checkRead();
}

double AlphaBeta2::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  const Vector d0=getSeparation(myatoms.getPosition(1),myatoms.getPosition(0));
  const Vector d1=getSeparation(myatoms.getPosition(2),myatoms.getPosition(1));
  const Vector d2=getSeparation(myatoms.getPosition(3),myatoms.getPosition(2));

  Vector dd0,dd1,dd2;
  PLMD::Torsion t;
  const double value  = t.compute(d0,d1,d2,dd0,dd1,dd2);
  const double vcos1  = cos(value-target1[tindex]);
  const double vcos2  = cos(value-target2[tindex]);
  const double vsin1  = sin((target1[tindex]-value)/2.0);
  const double vsin2  = sin((target2[tindex]-value)/2.0);
  const double v2cos  = 2.0-vcos1-vcos2;
  const double svalue = 8*csin12[tindex]*vsin1*vsin2/(v2cos*v2cos);
  const double cvalue = (vcos2-vcos1)/v2cos;

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
