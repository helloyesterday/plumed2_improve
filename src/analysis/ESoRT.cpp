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

#include "core/ActionPilot.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace analysis {

//+PLUMEDOC PRINTANALYSIS ESORT
/*
Does a committor analysis.

\par Examples

The following input monitors two torsional angles during a simulation,
defines two basins (A and B) as a function of the two torsions and
stops the simulation when it falls in one of the two. In the log
file will be shown the latest values for the CVs and the basin reached.
\plumedfile
TORSION ATOMS=1,2,3,4 LABEL=r1
TORSION ATOMS=2,3,4,5 LABEL=r2
ESORT ...
  ARG=r1,r2
  STRIDE=10
  BASIN_LL1=0.15,0.20
  BASIN_UL1=0.25,0.40
  BASIN_LL2=-0.25,-0.40
  BASIN_UL2=-0.15,-0.20
... ESORT
\endplumedfile

*/
//+ENDPLUMEDOC

class ESoRT :
  public ActionPilot,
  public ActionWithArguments
{
private:
  std::string file;
  OFile ofile;
  std::string fmt;
  std::vector<double> init_arg;
  std::vector< std::vector<double> > lowerlimits;
  std::vector< std::vector<double> > upperlimits;
  unsigned nbasins;
  unsigned basin;
  bool doNotStop;
  bool init_not_get;
public:
  static void registerKeywords( Keywords& keys );
  explicit ESoRT(const ActionOptions&ao);
  void calculate();
  void apply() {}
};

PLUMED_REGISTER_ACTION(ESoRT,"ESORT")

void ESoRT::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG");
  keys.add("numbered", "BASIN_LL","List of lower limits for basin #");
  keys.add("numbered", "BASIN_UL","List of upper limits for basin #");
  keys.reset_style("BASIN_LL","compulsory"); keys.reset_style("BASIN_UL","compulsory");
  keys.add("compulsory","STRIDE","1","the frequency with which the CVs are analysed");
  keys.add("optional","FILE","the name of the file on which to output the reached basin");
  keys.add("optional","FMT","the format that should be used to output real numbers");
  keys.addFlag("NOSTOP",false,"if true do not stop the simulation when reaching a basin but just keep track of it");
}

ESoRT::ESoRT(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithArguments(ao),
  fmt("%f"),
  basin(0),
  doNotStop(false),
  init_not_get(true)
{
  ofile.link(*this);
  parse("FILE",file);
  if(file.length()>0) {
    ofile.open(file);
    log.printf("  on file %s\n",file.c_str());
  } else {
    log.printf("  on plumed log file\n");
    ofile.link(log);
  }
  parse("FMT",fmt);
  fmt=" "+fmt;
  log.printf("  with format %s\n",fmt.c_str());
  init_arg.resize(getNumberOfArguments());

  for(unsigned b=1;; ++b ) {
    std::vector<double> tmpl, tmpu;
    parseNumberedVector("BASIN_LL", b, tmpl );
    parseNumberedVector("BASIN_UL", b, tmpu );
    if( tmpl.empty() && tmpu.empty() ) break;
    if( tmpl.size()!=getNumberOfArguments()) error("Wrong number of values for BASIN_LL: they should be equal to the number of arguments");
    if( tmpu.size()!=getNumberOfArguments()) error("Wrong number of values for BASIN_UL: they should be equal to the number of arguments");
    lowerlimits.push_back(tmpl);
    upperlimits.push_back(tmpu);
    nbasins=b;
  }

  parseFlag("NOSTOP", doNotStop);

  checkRead();


  for(unsigned b=0; b<nbasins; b++) {
    log.printf("  BASIN %u definition:\n", b+1);
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      if(lowerlimits[b][i]>upperlimits[b][i]) error("ESORT: UPPER bounds must always be greater than LOWER bounds");
      log.printf(" %f - %f\n", lowerlimits[b][i], upperlimits[b][i]);
    }
    if(doNotStop) log.printf(" COMMITOR will keep track of the visited basins without stopping the simulations\n");
  }

  for(unsigned i=0; i<getNumberOfArguments(); ++i) ofile.setupPrintValue( getPntrToArgument(i) );
}

void ESoRT::calculate() {
  std::vector<unsigned> inbasin;
  inbasin.assign (nbasins,1);
  
  if(init_not_get)
  {
	  for(unsigned i=0; i<getNumberOfArguments(); ++i)
	    init_arg[i]=getArgument(i);
	  init_not_get=false;
  }

  // check if current configuration belongs to a basin
  for(unsigned b=0; b<nbasins; ++b) {
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      if(getArgument(i)>lowerlimits[b][i]&&getArgument(i)<upperlimits[b][i]) {
        inbasin[b]*=1;
      } else {
        inbasin[b]*=0;
      }
    }
  }

  // check in which basin we are if any and if this is the same or a new one
  bool inonebasin = false;
  for(unsigned b=0; b<nbasins; ++b) {
    if(inbasin[b]==1) {
      if(basin!=(b+1)) {
        basin = b+1;
        ofile.fmtField(" %f");
        for(unsigned i=0; i<getNumberOfArguments(); i++) {
          ofile.fmtField(fmt);
          ofile.printField( "init_" + getPntrToArgument(i)->getName(), init_arg[i] );
          ofile.printField( getPntrToArgument(i), getArgument(i) );
        }
        ofile.printField("basin", static_cast<int> (b+1));
        ofile.printField("lifetime",getTime());
        ofile.printField();
      }
      inonebasin = true;
      break;
    }
  }
  if(!inonebasin) basin = 0;

  // then check if the simulation should be stopped
  if(inonebasin&&(!doNotStop)) {
    ofile.flush();
    plumed.stop();
  }
}

}
}
