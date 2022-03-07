/*
 * ProjectedTimeDilutedSpinColorDiagonal.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MNoise_ProjectedTimeDilutedSpinColorDiagonal_hpp_
#define Hadrons_MNoise_ProjectedTimeDilutedSpinColorDiagonal_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *             Generate time diluted spin-color diagonal noise                *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MNoise)

class ProjectedTimeDilutedSpinColorDiagonalPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ProjectedTimeDilutedSpinColorDiagonalPar,
                                    unsigned int, nsrc,
                                    std::string,  basis);
};

template <typename FImpl>
class TProjectedTimeDilutedSpinColorDiagonal: public Module<ProjectedTimeDilutedSpinColorDiagonalPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TProjectedTimeDilutedSpinColorDiagonal(const std::string name);
    // destructor
    virtual ~TProjectedTimeDilutedSpinColorDiagonal(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ProjectedTimeDilutedSpinColorDiagonal, TProjectedTimeDilutedSpinColorDiagonal<FIMPL>, MNoise);
MODULE_REGISTER_TMP(StagProjectedTimeDilutedSpinColorDiagonal, TProjectedTimeDilutedSpinColorDiagonal<STAGIMPL>, MNoise);
MODULE_REGISTER_TMP(ZProjectedTimeDilutedSpinColorDiagonal, TProjectedTimeDilutedSpinColorDiagonal<ZFIMPL>, MNoise);

/******************************************************************************
 *              TProjectedTimeDilutedSpinColorDiagonal implementation                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TProjectedTimeDilutedSpinColorDiagonal<FImpl>::TProjectedTimeDilutedSpinColorDiagonal(const std::string name)
: Module<ProjectedTimeDilutedSpinColorDiagonalPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TProjectedTimeDilutedSpinColorDiagonal<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().basis};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TProjectedTimeDilutedSpinColorDiagonal<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TProjectedTimeDilutedSpinColorDiagonal<FImpl>::setup(void)
{
    auto &basis = envGet(std::vector<FermionField>,par().basis);
    envCreateDerived(SpinColorDiagonalNoise<FImpl>, 
                     ProjectedTimeDilutedNoise<FImpl>,
                     getName(), 1, envGetGrid(FermionField), par().nsrc, basis);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TProjectedTimeDilutedSpinColorDiagonal<FImpl>::execute(void)
{
    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, getName());
    auto nt    = noise.getGrid()->GlobalDimensions()[Tp];

    LOG(Message) << "Generating time-diluted, spin-color diagonal noise" << std::endl;
    noise.generateNoise(rng4d());
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MNoise_ProjectedTimeDilutedSpinColorDiagonal_hpp_
