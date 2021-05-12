/*
 * Point.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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

#ifndef Hadrons_MSource_PointF_hpp_
#define Hadrons_MSource_PointF_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Point source
 ------------
 * src_x = delta_x,position
 
 * options:
 - position: space-separated integer sequence (e.g. "0 1 1 0")
 
 */

/******************************************************************************
 *                                  TPointF                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class PointFPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PointFPar,
                                    std::string, position);
};

template <typename FImpl>
class TPointF: public Module<PointFPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::FermionField FermionField;
    typedef typename FermionField::scalar_object FermionSite;
public:
    // constructor
    TPointF(const std::string name);
    // destructor
    virtual ~TPointF(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(StagPointF,   TPointF<STAGIMPL>,     MSource);

/******************************************************************************
 *                       TPointF template implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPointF<FImpl>::TPointF(const std::string name)
: Module<PointFPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPointF<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TPointF<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPointF<FImpl>::setup(void)
{
    auto nc = FImpl::Dimension;

    envCreate(std::vector<FermionField>, getName(), 1, nc, envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPointF<FImpl>::execute(void)
{
    LOG(Message) << "Creating point source at position [" << par().position
                << "]" << std::endl;

    FermionSite fSite;
    decltype(peekIndex<ColourIndex,FermionSite>(fSite,0)) id(1.);

    std::vector<int> position = strToVec<int>(par().position);
    auto             &src     = envGet(std::vector<FermionField>, getName());
    auto             nc       = FImpl::Dimension;

    if (position.size() != env().getNd())
    {
        HADRONS_ERROR(Size, "position has " + std::to_string(position.size())
                      + " components (must have " + std::to_string(env().getNd()) + ")");
    }

    for (int i=0;i<nc;i++) {
        src[i] = Zero();
        fSite = 0.;
        pokeIndex<ColourIndex,FermionSite>(fSite,id,i);
        pokeSite(fSite, src[i], position);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_Point_hpp_
