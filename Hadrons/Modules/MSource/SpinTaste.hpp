/*
 * SpinTaste.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2021
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>
 * Author: Michael Lynch <ml11@illinois.edu>
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

#ifndef Hadrons_MSource_SpinTaste_hpp_
#define Hadrons_MSource_SpinTaste_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                 SpinTaste                                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class SpinTastePar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(SpinTastePar,
                                  std::string, q,
                                  unsigned int, t0,
                                  std::string, gamma_spin);
};

template <typename FImpl, typename TField = typename FImpl::PropagatorField>
class TSpinTaste: public Module<SpinTastePar>
{
public:
  FERM_TYPE_ALIASES(FImpl,);
public:
  // constructor
  TSpinTaste(const std::string name);
  // destructor
  virtual ~TSpinTaste(void) {};
  // dependency relation
  virtual std::vector<std::string> getInput(void);
  virtual std::vector<std::string> getOutput(void);
protected:
  // setup
  virtual void setup(void);
  // execution
  virtual void execute(void);
private:
  void spinTasteOp(TField &out, const TField &in);
  bool hasPhase_{false};
  std::string phName_;
};

MODULE_REGISTER_TMP(SpinTasteProp, TSpinTaste<STAGIMPL>, MSource);
MODULE_REGISTER_TMP(SpinTaste, ARG(TSpinTaste<STAGIMPL, STAGIMPL::FermionField>), MSource);

/******************************************************************************
 *                          TSpinTaste implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename TField>
TSpinTaste<FImpl, TField>::TSpinTaste(const std::string name)
  : Module<SpinTastePar>(name)
, phName_ (name + "_sph")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename TField>
std::vector<std::string> TSpinTaste<FImpl, TField>::getInput(void)
{
  std::vector<std::string> in = {par().q};

  return in;
}

template <typename FImpl, typename TField>
std::vector<std::string> TSpinTaste<FImpl, TField>::getOutput(void)
{
  std::vector<std::string> out = {getName()};

  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename TField>
void TSpinTaste<FImpl, TField>::setup(void)
{
  if (envHasType(TField, par().q)) {

    envCreateLat(TField, getName());

  } else if (envHasType(std::vector<TField>, par().q)) {

    auto &q = envGet(std::vector<TField>, par().q);

    envCreate(std::vector<TField>, getName(), 1, q.size(),
              envGetGrid(TField));

  } else {
    HADRONS_ERROR(Argument, "incompatible field type '" + env().getObjectType(par().q) )
  }

  // envTmpLat(LatticeComplex, "stagPhaseDiag");
  envCache(LatticeComplex, phName_,    1, envGetGrid(LatticeComplex));
}

// spin-taste operator ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename TField>
void TSpinTaste<FImpl, TField>::spinTasteOp(TField &out, const TField &in)
{

  out = Zero();

  // grid can't handle real * prop, so use complex
  auto &ph = envGet(LatticeComplex, phName_);

  if (!hasPhase_) {

    Lattice<iScalar<vInteger> > x(env().getGrid()); LatticeCoordinate(x, 0);
    Lattice<iScalar<vInteger> > y(env().getGrid()); LatticeCoordinate(y, 1);
    Lattice<iScalar<vInteger> > z(env().getGrid()); LatticeCoordinate(z, 2);

    // local taste non-singlet ops, including ``Hermiticity" phase,
    // see Tab. 11.2 in Degrand and Detar

    ph = 1.0;

    std::vector<Gamma::Algebra> gamDiag = strToVec<Gamma::Algebra>(par().gamma_spin);

    if (gamDiag.size() != 1) {
      HADRONS_ERROR(Argument, "Please only provide one spin-taste operation to '" + getName() )
    }

    LOG(Message) << "Using gamma: " << gamDiag[0] << std::endl;

    switch (gamDiag[0]) {

    case Gamma::Algebra::GammaX  :
      ph = where( mod(x, 2) == (Integer)0, ph, -ph);
      break;

    case Gamma::Algebra::GammaY  :
      ph = where( mod(y, 2) == (Integer)0, ph, -ph);
      break;

    case Gamma::Algebra::GammaZ  :
      ph = where( mod(z, 2) == (Integer)0, ph, -ph);
      break;

    case Gamma::Algebra::Gamma5  :
      break;

    default :
      HADRONS_ERROR(Argument, "Your choice of spin-taste operation for '" + getName() + "' is not currently supported'" )
    }
    hasPhase_ = true;
  }
  out = in * ph;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename TField>
void TSpinTaste<FImpl, TField>::execute(void)
{
  LOG(Message) << "Spin-taste acting on '" << par().q
               << "' using t0 '"          << par().t0
               << "' using gamma_spin '" << par().gamma_spin << "'."
               << std::endl;

  if (envHasType(TField, par().q))
  {
      auto  &src = envGet(TField, getName()); 
      auto  &q   = envGet(TField, par().q);

      LOG(Message) << "Using source '" << par().q << "'" << std::endl;
      spinTasteOp(src, q);
  }
  else
  {
      auto  &src = envGet(std::vector<TField>, getName()); 
      auto  &q   = envGet(std::vector<TField>, par().q);

      for (unsigned int i = 0; i < q.size(); ++i)
      {
          LOG(Message) << "Using element " << i << " of source vector '" 
                       << par().q << "'" << std::endl;
          spinTasteOp(src[i], q[i]);
      }
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_SpinTaste_hpp_
