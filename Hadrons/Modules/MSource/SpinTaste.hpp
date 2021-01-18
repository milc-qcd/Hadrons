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
				    std::string, t0,
                                    std::string, gamma_spin);
};

template <typename FImpl>
class TSpinTaste: public Module<SpinTastePar>
{
public:
    BASIC_TYPE_ALIASES(STAGIMPL,);
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
};

MODULE_REGISTER_TMP(SpinTaste, TSpinTaste<STAGIMPL>, MSource);

/******************************************************************************
 *                          TSpinTaste implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSpinTaste<FImpl>::TSpinTaste(const std::string name)
: Module<SpinTastePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSpinTaste<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSpinTaste<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSpinTaste<FImpl>::setup(void)
{
  envCreateLat(PropagatorField, getName());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSpinTaste<FImpl>::execute(void)
{
    LOG(Message) << "Spin-taste acting on '" << par().q
                 << "' using t0 '"          << par().t0
                 << "' using gamma_spin '" << par().gamma_spin << "'."
                 << std::endl;

    auto &q          = envGet(PropagatorField, par().q);
    int t0           = stoi(par().t0);
    auto &out        = envGet(PropagatorField, getName());
    
    out = Zero();
    out = q;  // For testing
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_SpinTaste_hpp_
