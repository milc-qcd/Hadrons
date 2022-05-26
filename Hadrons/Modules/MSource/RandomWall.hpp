/*
 * RandomWall.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: fionnoh <fionnoh@gmail.com>
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

#ifndef Hadrons_MSource_RandomWall_hpp_
#define Hadrons_MSource_RandomWall_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Random Wall source
 -----------------------------
 
 * options:
 - tW:   source timeslice (integer)
 - size: number of sources (integer)
 
 */

/******************************************************************************
 *                         Random Wall                                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class RandomWallPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RandomWallPar,
                                    unsigned int, tStep,
                                    unsigned int, nSrc,
                                    std::string, reuset0);
};

template <typename FImpl>
class TRandomWall: public Module<RandomWallPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TRandomWall(const std::string name);
    // destructor
    virtual ~TRandomWall(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool reuset0_ = false;
};

// MODULE_REGISTER_TMP(RandomWall, TRandomWall<FIMPL>, MSource);
MODULE_REGISTER_TMP(StagRandomWall, TRandomWall<STAGIMPL>, MSource);

/******************************************************************************
 *                 TRandomWall implementation                                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TRandomWall<FImpl>::TRandomWall(const std::string name)
: Module<RandomWallPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TRandomWall<FImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TRandomWall<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName()+"_shift"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandomWall<FImpl>::setup(void)
{
    envTmp(TimeDilutedNoise<FImpl>, "noise", 1, envGetGrid(FermionField), par().nSrc);
    envTmp(PropagatorField, "shiftedField", 1, envGetGrid(PropagatorField));

    envCreate(std::vector<PropagatorField>, getName(), 1, 0, envGetGrid(PropagatorField));

    envCreate(std::vector<Integer>, getName()+"_shift", 1, 0, 0);

    if (!par().reuset0.empty()) {
        if (!(std::istringstream(par().reuset0) >> reuset0_)) {
            LOG(Error) << "parameter reuset0='" << par().reuset0 << "' must be 'true' or 'false'";
        }
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandomWall<FImpl>::execute(void)
{    
    envGetTmp(TimeDilutedNoise<FImpl>, noise);
    LOG(Message) << "Generating " << par().nSrc << " time-diluted, spin-color diagonal noise sources at every " << par().tStep << " time step(s)" << std::endl;
    noise.generateNoise(rng4d());

    auto &noisevec = envGet(std::vector<PropagatorField>,getName());
    auto &time_shift = envGet(std::vector<Integer>,getName()+"_shift");

    int nt    = envGetGrid(PropagatorField)->GlobalDimensions()[Tp];

    int tStep = par().tStep;
    int nSources = par().nSrc;

    int nsc   = noise.getNsc();
    int nSlices = nt/tStep;
    int nVecs = nSources*nSlices;


    time_shift.resize(nVecs,0);

    noisevec.resize(nVecs,envGetGrid(PropagatorField));

    envGetTmp(PropagatorField,shiftedField);

    for (int i=0;i<nSources;i++) {
        if (reuset0_) {
            shiftedField = noise.getProp(i*nt);
            noisevec[i*nSlices] = shiftedField;
        }
        for (int j=0;j<nSlices;j++) {
            int idx = i*nSlices+j;
            int offset = i*nt+j*tStep;
            if (!reuset0_) {
                noisevec[idx] = noise.getProp(offset);                
            } else {
                if (j != 0) {
                    noisevec[idx] = Cshift(noisevec[idx-1],Tp,tStep);
                }
            }
            time_shift[idx] = j*tStep;
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_RandomWall_hpp_
