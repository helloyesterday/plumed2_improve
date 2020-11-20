/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DIPOLE
/*
Calculate the dipole moment for a group of atoms.

When running with periodic boundary conditions, the atoms should be
in the proper periodic image. This is done automatically since PLUMED 2.5,
by considering the ordered list of atoms and rebuilding the molecule with a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED.

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

\par Examples

The following tells plumed to calculate the dipole of the group of atoms containing
the atoms from 1-10 and print it every 5 steps
\plumedfile
d: DIPOLE GROUP=1-10
PRINT FILE=output STRIDE=5 ARG=d
\endplumedfile

\attention
If the total charge Q of the group in non zero, then a charge Q/N will be subtracted to every atom,
where N is the number of atoms. This implies that the dipole (which for a charged system depends
on the position) is computed on the geometric center of the group.


*/
//+ENDPLUMEDOC

class Average_Dipole : public Colvar {
  vector<AtomNumber> ga_lista;
  bool components;
  bool nopbc;
  unsigned nq;
  unsigned nmol;
  double ctot;
  double cmol;
  std::vector<double> charges;
  std::vector<bool> calc_charge;
  Vector zero_vector;
public:
  explicit Average_Dipole(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(Average_Dipole,"AVERAGE_DIPOLE")

void Average_Dipole::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.add("atoms","GROUP","the group of atoms we are calculating the dipole moment for");
  keys.add("compulsory","ATOM_CHARGES","the charges of each atom in molecule");
  keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the dipole separately and store them as label.x, label.y and label.z");
  keys.addOutputComponent("x","COMPONENTS","the x-component of the dipole");
  keys.addOutputComponent("y","COMPONENTS","the y-component of the dipole");
  keys.addOutputComponent("z","COMPONENTS","the z-component of the dipole");
}

Average_Dipole::Average_Dipole(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  components(false),
  ctot(0),
  cmol(0)
{
  parseAtomList("GROUP",ga_lista);
  parseFlag("COMPONENTS",components);
  parseFlag("NOPBC",nopbc);
  std::vector<double> _charges;
  parseVector("ATOM_CHARGES",_charges);
  checkRead();
  if(components) {
    addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
    addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
    addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
  } else {
    addValueWithDerivatives(); setNotPeriodic();
  }
  zero_vector.zero();
  nq=_charges.size();
  for(unsigned i=0;i!=nq;++i)
      cmol+=_charges[i];

  charges=_charges;
  calc_charge.assign(nq,false);
  for(unsigned i=0;i!=nq;++i)
  {
     charges[i]-=cmol;
     if(fabs(charges[i])>1.0e-8)
       calc_charge[i]=true;
  }
  
  plumed_massert(getNumberOfAtoms()%nq==0,"the number of atoms can not be divisible by the number of molecules!");
  nmol=ga_lista.size()/nq;
  ctot=cmol*nmol;

  log.printf("  of %u atoms\n",static_cast<unsigned>(ga_lista.size()));
  for(unsigned int i=0; i<ga_lista.size(); ++i) {
    log.printf("  %d", ga_lista[i].serial());
  }
  log.printf("  \n");
  
  log.printf("  with number of molecules: %d\n",nmol);
  for(unsigned i=0;i!=nq;++i)
    log.printf("   Atom %d with charge %f\n",int(i+1),_charges[i]);
  log.printf("  with total charge in molecules: %f\n",cmol);
  
  if(nopbc) log.printf("  without periodic boundary conditions\n");
  else      log.printf("  using periodic boundary conditions\n");

  requestAtoms(ga_lista);
}

// calculator
void Average_Dipole::calculate()
{
  if(!nopbc) makeWhole();
  //~ double ctot=0.;
  unsigned N=getNumberOfAtoms();

  Vector dipje;
  dipje.zero();
  for(unsigned i=0;i!=nmol;++i)
  {
      for(unsigned j=0;j!=nq;++j)
      {
          if(calc_charge[j])
          {
              unsigned iatom=i*nq+j;
              Vector pos=getPosition(iatom);
              dipje += charges[j]*pos;
          }
      }
  }
  dipje/=nmol;

  if(!components) {
    double dipole = dipje.modulo();
    double idip = 1./dipole;

    for(unsigned i=0;i!=nmol;++i)
    {
        for(unsigned j=0;j!=nq;++j)
        {
            unsigned iatom=i*nq+j;
            if(calc_charge[j])
            {
                double dfunc=charges[j]*idip;
                setAtomsDerivatives(iatom,dfunc*dipje);
            }
            else
                setAtomsDerivatives(iatom,zero_vector);
        }
    }
    setBoxDerivativesNoPbc();
    setValue(dipole);
  } else {
    Value* valuex=getPntrToComponent("x");
    Value* valuey=getPntrToComponent("y");
    Value* valuez=getPntrToComponent("z");
    for(unsigned i=0; i<N; i++) {
      setAtomsDerivatives(valuex,i,charges[i]*Vector(1.0,0.0,0.0));
      setAtomsDerivatives(valuey,i,charges[i]*Vector(0.0,1.0,0.0));
      setAtomsDerivatives(valuez,i,charges[i]*Vector(0.0,0.0,1.0));
    }
    setBoxDerivativesNoPbc(valuex);
    setBoxDerivativesNoPbc(valuey);
    setBoxDerivativesNoPbc(valuez);
    valuex->set(dipje[0]);
    valuey->set(dipje[1]);
    valuez->set(dipje[2]);
  }
}

}
}
