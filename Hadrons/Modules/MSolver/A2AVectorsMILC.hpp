/*
 * A2AVectorsMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
 * Author: fionnoh <fionnoh@gmail.com>
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
#ifndef Hadrons_MSolver_A2AVectorsMILC_hpp_
#define Hadrons_MSolver_A2AVectorsMILC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Create all-to-all V & W vectors                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2AVectorsMILCPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsMILCPar,
                                  std::string, noise,
                                  std::string, action,
                                  std::string, eigenPack,
                                  std::string, solver,
                                  std::string, highOutput,
                                  bool,        highMultiFile);
};

template <typename FImpl, typename Pack>
class TA2AVectorsMILC : public Module<A2AVectorsMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef A2AVectorsSchur<FImpl> A2A;

public:
    // constructor
    TA2AVectorsMILC(const std::string name);
    // destructor
    virtual ~TA2AVectorsMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);

    // setup
    virtual void setup(void);

    // execute
    virtual void execute(void);
private:
    std::string  solverName_;
    unsigned int Nh_{0};
    bool hasEpack_{false};
};

MODULE_REGISTER_TMP(A2AVectorsMILC, 
    ARG(TA2AVectorsMILC<STAGIMPL, BaseFermionEigenPack<STAGIMPL,Complex>>), MSolver);

/******************************************************************************
 *                       TA2AVectorsMILC implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TA2AVectorsMILC<FImpl, Pack>::TA2AVectorsMILC(const std::string name)
: Module<A2AVectorsMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectorsMILC<FImpl, Pack>::getInput(void)
{
    std::vector<std::string> in {par().action,par().solver, par().noise};

    if (!par().eigenPack.empty()) {
        in.push_back(par().eigenPack);
        in.push_back(par().eigenPack+"_evec");
        in.push_back(par().eigenPack+"_eval");
    }
    
    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectorsMILC<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {};

    if (!par().noise.empty()) {
        out.push_back(getName() + "_w");
        out.push_back(getName() + "_v");
    }

    return out;
}

/******************************************************************************
 *              TA2AVectorsMILC setup                                         *
 ******************************************************************************/
template <typename FImpl, typename Pack>
void TA2AVectorsMILC<FImpl, Pack>::setup(void)
{
    auto        &action     = envGet(FMat, par().action);
    auto        &solver     = envGet(Solver, par().solver);
    int         Ls          = env().getObjectLs(par().action);

    envTmp(A2A, "a2a", 1, action, solver);

    hasEpack_ = !par().eigenPack.empty();

    if (Ls > 1) {
       HADRONS_ERROR(Argument, "Ls > 1 not implemented");
    }


    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);

    Nh_ = noise.fermSize();

    envCreate(std::vector<FermionField>, getName() + "_v", 1, 
              Nh_, envGetGrid(FermionField));
    envCreate(std::vector<FermionField>, getName() + "_w", 1, 
              Nh_, envGetGrid(FermionField));
}

/******************************************************************************
 *              TA2AVectorsMILC execution                                     *
 ******************************************************************************/
template <typename FImpl, typename Pack>
void TA2AVectorsMILC<FImpl, Pack>::execute(void)
{
    auto        &action    = envGet(FMat, par().action);
    auto        &solver    = envGet(Solver, par().solver);

    envGetTmp(A2A, a2a);

    if (hasEpack_)
    {
       LOG(Message) << "Computing all-to-all vectors "
                    << " using eigenpack '" << par().eigenPack << "' (low modes) and noise '"
                    << par().noise << "' (" << Nh_ 
                    << " noise vectors)" << std::endl;
    } else {
       LOG(Message) << "Computing all-to-all vectors "
                    << " using noise '" << par().noise << "' (" << Nh_ 
                    << " noise vectors)" << std::endl;
    }

    typename std::vector<Real>::iterator it_eval;

    auto &v     = envGet(std::vector<FermionField>, getName() + "_v");
    auto &w     = envGet(std::vector<FermionField>, getName() + "_w");
    auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);

    int nsrc = noise.size()/FImpl::Dimension;  

    // Normalization for the noise sources
    RealD norm = 1.0/::sqrt(Real(nsrc));

    std::cout << "Normalizing stochastic vectors by 1/sqrt(" << nsrc << ")" << std::endl;

    startTimer("W high mode");
    for (int ih = 0; ih < Nh_; ih++) {
        w[ih] = noise.getFerm(ih);
        w[ih] = norm*w[ih];
    }
    if (hasEpack_) {
        LOG(Message) << "Projecting low contribution from stochastic high mode sources" << std::endl;

        auto &epack = envGet(Pack, par().eigenPack);
        a2a.removeLowModeProj(w,epack.evec,epack.eval);
    }
    stopTimer("W high mode");

    for (int ih = 0; ih < Nh_; ih++)
    {
       startTimer("V high mode");
       LOG(Message) << "V vector (solve) i = " << ih
                    << " (" << ((hasEpack_) ? "high " : "") 
                    << "stochastic mode)" << std::endl;

       a2a.makeHighModeV(v[ih],w[ih]);

       stopTimer("V high mode");
    }

    if (!par().highOutput.empty())
    {
       startTimer("V I/O");
       A2AVectorsIo::write(par().highOutput + "_v", v, par().highMultiFile, vm().getTrajectory());
       stopTimer("V I/O");
       startTimer("W I/O");
       A2AVectorsIo::write(par().highOutput + "_w", w, par().highMultiFile, vm().getTrajectory());
       stopTimer("W I/O");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectorsMILC_hpp_
