/*
 * TimeDilutedSpinColorDiagonal.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
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
#ifndef Hadrons_MNoise_TimeDilutedSpinColorDiagonal_hpp_
#define Hadrons_MNoise_TimeDilutedSpinColorDiagonal_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *             Generate time diluted spin-color diagonal noise                *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNoise)

class TimeDilutedSpinColorDiagonalPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(TimeDilutedSpinColorDiagonalPar,
                                    unsigned int, nsrc);
};

template <typename FImpl>
class TTimeDilutedSpinColorDiagonal: public Module<TimeDilutedSpinColorDiagonalPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TTimeDilutedSpinColorDiagonal(const std::string name);
    // destructor
    virtual ~TTimeDilutedSpinColorDiagonal(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(TimeDilutedSpinColorDiagonal, TTimeDilutedSpinColorDiagonal<FIMPL>, MNoise);
MODULE_REGISTER_TMP(StagTimeDilutedSpinColorDiagonal, TTimeDilutedSpinColorDiagonal<STAGIMPL>, MNoise);
MODULE_REGISTER_TMP(ZTimeDilutedSpinColorDiagonal, TTimeDilutedSpinColorDiagonal<ZFIMPL>, MNoise);

/******************************************************************************
 *              TTimeDilutedSpinColorDiagonal implementation                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TTimeDilutedSpinColorDiagonal<FImpl>::TTimeDilutedSpinColorDiagonal(const std::string name)
: Module<TimeDilutedSpinColorDiagonalPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TTimeDilutedSpinColorDiagonal<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TTimeDilutedSpinColorDiagonal<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName()+"_vec", getName()+"_shift"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTimeDilutedSpinColorDiagonal<FImpl>::setup(void)
{
    envCreateDerived(SpinColorDiagonalNoise<FImpl>, 
                     TimeDilutedNoise<FImpl>,
                     getName(), 1, envGetGrid(FermionField), par().nsrc);

    envCreate(std::vector<FermionField>, getName() + "_vec", 1, 0, envGetGrid(FermionField));

    envCreate(std::vector<Integer>, getName()+"_shift", 1, 0, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TTimeDilutedSpinColorDiagonal<FImpl>::execute(void)
{

    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, getName());
    LOG(Message) << "Generating time-diluted, spin-color diagonal noise" << std::endl;
    noise.generateNoise(rng4d());

    auto &noisevec = envGet(std::vector<FermionField>,getName()+"_vec");

    int nferm = noise.fermSize();
    int nsc   = noise.getNsc();
    int nt    = envGetGrid(FermionField)->GlobalDimensions()[Tp];

    noisevec.resize(nferm,envGetGrid(FermionField));
    for (int i=0;i<nferm;i++) {
        noisevec[i] = noise.getFerm(i);
    }

    auto &time_shift = envGet(std::vector<Integer>,getName()+"_shift");

    time_shift.resize(nferm,0);

    for (int i = 0;i<time_shift.size();i++) {
        time_shift[i] = (i/nsc)%nt;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNoise_TimeDilutedSpinColorDiagonal_hpp_
