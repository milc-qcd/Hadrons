/*
 * SpinTasteMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2021
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>
 * Author: Michael Lynch <michaellynch628@gmail.com>
 * Author: Carleton DeTar <detar@physics.utah.edu>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */

#ifndef Hadrons_MFermion_SpinTasteMILC_hpp_
#define Hadrons_MFermion_SpinTasteMILC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                 SpinTasteMILC                                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

template <typename FImpl>
class TSpinTasteMILC: public Module<NoPar>
{
public:
  FERM_TYPE_ALIASES(FImpl,);
  typedef std::map<Gamma::Algebra, LatticeComplex> PhaseMap;
  typedef std::function<LatticeComplex (Gamma::Algebra gamma)> GammaFn;
public:
  // constructor
  TSpinTasteMILC(const std::string name);
  // destructor
  virtual ~TSpinTasteMILC(void) {};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
protected:
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
private:
  bool hasPhase_{false};
  std::string phName_;
};

MODULE_REGISTER_TMP(SpinTaste, TSpinTasteMILC<STAGIMPL>, MFermion);

/******************************************************************************
 *                          TSpinTasteMILC implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSpinTasteMILC<FImpl>::TSpinTasteMILC(const std::string name)
  : Module<NoPar>(name)
, phName_ (name + "_sph")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSpinTasteMILC<FImpl>::getInput(void)
{
  std::vector<std::string> in;

  return in;
}

template <typename FImpl>
std::vector<std::string> TSpinTasteMILC<FImpl>::getOutput(void)
{
  std::vector<std::string> out = {getName()};

  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSpinTasteMILC<FImpl>::setup(void)
{
    envCreate(GammaFn, getName(), 1, nullptr);
    PhaseMap dummy;
    envCache(PhaseMap, phName_, 1, dummy);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSpinTasteMILC<FImpl>::execute(void)
{

  std::vector<Gamma::Algebra> keys = {
      Gamma::Algebra::Identity,
      Gamma::Algebra::GammaX,
      Gamma::Algebra::GammaY,
      Gamma::Algebra::GammaZ,
      Gamma::Algebra::Gamma5
  };


  auto &stag_phase = envGet(PhaseMap,phName_);

  for (const auto &key:keys) {
      stag_phase.insert({key,envGetGrid(LatticeComplex)});
      stag_phase.at(key) = 1.0;
  }

  Lattice<iScalar<vInteger> > x(env().getGrid()); LatticeCoordinate(x,0);
  Lattice<iScalar<vInteger> > y(env().getGrid()); LatticeCoordinate(y,1);
  Lattice<iScalar<vInteger> > z(env().getGrid()); LatticeCoordinate(z,2);
  Lattice<iScalar<vInteger> > t(env().getGrid()); LatticeCoordinate(t,3);
  Lattice<iScalar<vInteger> > lin_5(env().getGrid()); lin_5=x+y+z+t;
  
  stag_phase.at(Gamma::Algebra::Identity) = where( mod(lin_5,2)==(Integer)0, stag_phase.at(Gamma::Algebra::Identity), -stag_phase.at(Gamma::Algebra::Identity));
  stag_phase.at(Gamma::Algebra::GammaX)   = where( mod(x,2)==(Integer)0, stag_phase.at(Gamma::Algebra::GammaX), -stag_phase.at(Gamma::Algebra::GammaX));
  stag_phase.at(Gamma::Algebra::GammaY)   = where( mod(y,2)==(Integer)0, stag_phase.at(Gamma::Algebra::GammaY), -stag_phase.at(Gamma::Algebra::GammaY));
  stag_phase.at(Gamma::Algebra::GammaZ)   = where( mod(z,2)==(Integer)0, stag_phase.at(Gamma::Algebra::GammaZ), -stag_phase.at(Gamma::Algebra::GammaZ));

  auto spinOp = [this](Gamma::Algebra gamma) {

    LatticeComplex result(envGetGrid(LatticeComplex));

    auto &stagPh = envGet(PhaseMap, phName_);

    if (stagPh.find(gamma) == stagPh.end()) {
      std::string gammaStr = Gamma::name[gamma];
      HADRONS_ERROR(Implementation,"The gamma operator '" + gammaStr + "' is not supported for stag fields");
    }

    result = stagPh.at(gamma);

    return result;
  };

  envGet(GammaFn, getName()) = spinOp;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_SpinTasteMILC_hpp_
