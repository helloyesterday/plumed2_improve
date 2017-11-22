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
#include "tools/Pbc.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR HBENERGY
/*
Calculate a torsional angle.

This command can be used to compute the torsion between four atoms or alternatively
to calculate the angle between two vectors projected on the plane
orthogonal to an axis. 

\par Examples

This input tells plumed to print the torsional angle between atoms 1, 2, 3 and 4
on file COLVAR.
\verbatim
t: TORSION ATOMS=1,2,3,4
# this is an alternative, equivalent, definition:
# t: TORSION VECTOR1=2,1 AXIS=2,3 VECTOR2=3,4
PRINT ARG=t FILE=COLVAR
\endverbatim

If you are working with a protein you can specify the special named torsion angles \f$\phi\f$, \f$\psi\f$, \f$\omega\f$ and \f$\chi_1\f$
by using TORSION in combination with the \ref MOLINFO command.  This can be done by using the following 
syntax.

\verbatim
MOLINFO MOLTYPE=protein STRUCTURE=myprotein.pdb
t1: TORSION ATOMS=@phi-3
t2: TORSION ATOMS=@psi-4
PRINT ARG=t1,t2 FILE=colvar STRIDE=10
\endverbatim

Here, \@phi-3 tells plumed that you would like to calculate the \f$\phi\f$ angle in the third residue of the protein.
Similarly \@psi-4 tells plumed that you want to calculate the \f$\psi\f$ angle of the 4th residue of the protein.
*/
//+ENDPLUMEDOC
   
class HBEnergy : public Colvar {
  bool pbc;
  double k_hb;

public:
  explicit HBEnergy(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(HBEnergy,"HBENERGY")

void HBEnergy::registerKeywords(Keywords& keys){
   Colvar::registerKeywords( keys );
   keys.add("atoms","ATOMS","the four atoms involved in the hydrogen bond in the order of C O N H");
}

HBEnergy::HBEnergy(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
k_hb(11.66834)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  
  if(atoms.size()!=4)
    error("Number of specified atoms should be 4");

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  checkRead();

  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  
  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
}

// calculator
void HBEnergy::calculate(){

  // 0:C 1:O 2:N 3:H
  Vector v_ON,v_CH,v_OH,v_CN;
  
  if(pbc) makeWhole();
  
  v_ON=delta(getPosition(1),getPosition(2));
  v_CH=delta(getPosition(0),getPosition(3));
  v_OH=delta(getPosition(1),getPosition(3));
  v_CN=delta(getPosition(0),getPosition(2));
  
  double d_ON,d_CH,d_OH,d_CN;
  
  d_ON=v_ON.modulo();
  d_CH=v_CH.modulo();
  d_OH=v_OH.modulo();
  d_CN=v_CN.modulo();
  
  double value=k_hb*(1.0/d_ON+1.0/d_CH-1.0/d_OH-1.0/d_CN);
  
  double _d3_ON,_d3_CH,_d3_OH,_d3_CN;
  _d3_ON=k_hb/std::pow(d_ON,3);
  _d3_CH=k_hb/std::pow(d_CH,3);
  _d3_OH=k_hb/std::pow(d_OH,3);
  _d3_CN=k_hb/std::pow(d_CN,3);
  
  setAtomsDerivatives(0,+v_CH*_d3_CH+v_CN*_d3_CN);
  setAtomsDerivatives(1,+v_ON*_d3_ON+v_OH*_d3_OH);
  setAtomsDerivatives(2,-v_ON*_d3_ON-v_CN*_d3_CN);
  setAtomsDerivatives(3,-v_CH*_d3_CH-v_OH*_d3_OH);

  setBoxDerivativesNoPbc();
  setValue(value);
}

}
}



