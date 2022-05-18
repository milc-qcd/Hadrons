/*
 * RandomPoint.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: Michael Lynch <michaellynch628@gmail.com>
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

#ifndef Hadrons_MSource_RandomPoint_hpp_
#define Hadrons_MSource_RandomPoint_hpp_

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
 *                                  TRandomPoint                                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class RandomPointPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RandomPointPar,
                                    int,  nSrc,
                                    bool, uniqueT);
};

template <typename FImpl>
class TRandomPoint: public Module<RandomPointPar>
{
public:
    BASIC_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TRandomPoint(const std::string name);
    // destructor
    virtual ~TRandomPoint(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(StagRandomPoint,   TRandomPoint<STAGIMPL>,     MSource);

/******************************************************************************
 *                       TRandomPoint template implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TRandomPoint<FImpl>::TRandomPoint(const std::string name)
: Module<RandomPointPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TRandomPoint<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TRandomPoint<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName()+"_shift"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandomPoint<FImpl>::setup(void)
{
    envCreate(std::vector<PropagatorField>, getName(),1,par().nSrc,env().getGrid());
    envCreate(std::vector<Integer>, getName()+"_shift", 1, par().nSrc, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandomPoint<FImpl>::execute(void)
{
    LOG(Message) << "Creating " << par().nSrc << " points at random positions on the lattice." << std::endl;

    int nt = env().getDim(Nd-1);

    std::vector<int> position(Nd,0);
    std::vector<int> times(nt);
    SitePropagator   id;
    RealD rnum;

    auto rng = env().getSerialRng();

    auto &src        = envGet(std::vector<PropagatorField>, getName());
    auto &time_shift = envGet(std::vector<Integer>,getName()+"_shift");

    if (par().uniqueT && par().nSrc > nt) {
        LOG(Error) << "Requested unique time indices, but requested sources (" << par().nSrc 
                    << ") > lattice time slices (" << nt << ")." << std::endl;
    }

    for (int i=0;i<nt;i++) {
        times[i] = i;
    }

    id = 1.;
    for (int i=0;i<par().nSrc;i++) {
        src[i] = Zero();

        for (int j=0;j<Nd;j++) {
            auto N = env().getDim(j);
            rng->fill(rnum,rng->_uniform);
            int idx = int(N*rnum);
            if (par().uniqueT && j == Nd-1) {
                N = times.size();
                int tidx = int(N*rnum);
                idx = times[tidx];
                times.erase(times.begin()+tidx);
            }
            position[j] = idx;
        }
        time_shift[i] = position[Nd-1];
        
        pokeSite(id, src[i], position);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_RandomPoint_hpp_
