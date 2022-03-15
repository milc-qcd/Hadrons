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
                                  std::string, lowModes,
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
    unsigned int Nl_{0},Nh_{0};
    bool hasEpack_{false}, usesMultiRHS{false};
};

MODULE_REGISTER_TMP(A2AVectorsMILC, 
    ARG(TA2AVectorsMILC<STAGIMPL, BaseFermionEigenPack<STAGIMPL>>), MSolver);

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
    std::vector<std::string> in {par().action,par().solver};

    if (!par().eigenPack.empty()) {

        in.push_back(par().eigenPack);

        if (IsStaggeredImpl<FImpl>()) {
           in.push_back(par().eigenPack+"_mass");
           in.push_back(par().eigenPack+"_evenEigen");
        }
    }
    
    if (!par().noise.empty())
        in.push_back(par().noise);

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

    if (!par().eigenPack.empty()) {
        out.push_back(getName() + "_evec");
        out.push_back(getName() + "_eval");
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
    int         actionLs    = env().getObjectLs(par().action);
    int         solverLs    = env().getObjectLs(par().solver);

    envTmp(A2A, "a2a", 1, action, solver);

    hasEpack_ = !par().eigenPack.empty();
    if (hasEpack_) {

        auto &epack = envGet(Pack, par().eigenPack);
        Nl_ = epack.evec.size()*(IsStaggeredImpl<FImpl>()?2:1);
        envCreate(std::vector<FermionField>, getName() + "_evec", actionLs, Nl_, envGetGrid(FermionField, actionLs));
        envCreate(std::vector<ComplexD>, getName() + "_eval", actionLs, Nl_, 0);
    }

    if (actionLs > 1) {
       HADRONS_ERROR(Argument, "Ls > 1 not implemented");
    }


    if (!par().noise.empty()) {

        auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);

        Nh_ = noise.fermSize();

        if (solverLs == Nh_) {
            usesMultiRHS = true;
            envTmpLat(FermionField,"multiRHSource",solverLs);
            envTmpLat(FermionField,"multiRHSolve",solverLs);
        }

        envCreate(std::vector<FermionField>, getName() + "_v", 1, 
                  Nh_, envGetGrid(FermionField));
        envCreate(std::vector<FermionField>, getName() + "_w", 1, 
                  Nh_, envGetGrid(FermionField));
    }
}

/******************************************************************************
 *              TA2AVectorsMILC execution                                     *
 ******************************************************************************/
template <typename FImpl, typename Pack>
void TA2AVectorsMILC<FImpl, Pack>::execute(void)
{
    auto        &action    = envGet(FMat, par().action);
    auto        &solver    = envGet(Solver, par().solver);
    int         solverLs   = env().getObjectLs(par().solver);

    Real        mass;
    envGetTmp(A2A, a2a);

    if (Nl_ > 0)
    {
       if (Nh_ > 0)
       {
           LOG(Message) << "Computing all-to-all vectors "
                        << " using eigenpack '" << par().eigenPack << "' ("
                        << Nl_ << " low modes) and noise '"
                        << par().noise << "' (" << Nh_ 
                        << " noise vectors)" << std::endl;
       } else {
           LOG(Message) << "Computing all-to-all vectors "
                        << " using eigenpack '" << par().eigenPack << "' ("
                        << Nl_ << " low modes) and zero noise vectors" << std::endl;

       }
       LOG(Message) << "Eigenpack with conjugate pair evecs and corresponding evals available in '" 
                    << getName() << "_evec,' and " << getName() << "_eval,' respectively"<< std::endl;

    } else {

       LOG(Message) << "Computing all-to-all vectors "
                    << " using noise '" << par().noise << "' (" << Nh_ 
                    << " noise vectors)" << std::endl;
    }

    typename std::vector<FermionField>::iterator it_evec, it_lowModeEvec;
    typename std::vector<Real>::iterator it_eval;
    typename std::vector<ComplexD>::iterator it_lowModeEval;

    if (Nl_ > 0) {

       auto &lowModeVecs = envGet(std::vector<FermionField>, getName() + "_evec");
       auto &lowModeVals = envGet(std::vector<ComplexD>, getName() + "_eval");
       auto &epack  = envGet(Pack, par().eigenPack);
       it_evec = epack.evec.begin();
       it_lowModeEvec = lowModeVecs.begin();
       it_lowModeEval = lowModeVals.begin();

       mass = (envGet(std::vector<Real>, par().eigenPack+"_mass"))[0];

       // Low modes
       for (auto it_eval = epack.eval.begin(); it_eval < epack.eval.end(); it_eval++)
       {
           int il = it_eval-epack.eval.begin();

           auto cbEven = (envGet(std::vector<bool>, par().eigenPack+"_evenEigen"))[0];
           startTimer("low mode pair");
           LOG(Message) << "Generating eigenvector pairs for i = " << 2*il << " and " << 2*il+1 << " (low mode)" << std::endl;

           a2a.makeLowModePairs(it_lowModeEvec, it_lowModeEval, it_evec, mass, *it_eval, cbEven);

           stopTimer("low mode pair");

           it_lowModeEvec+=2;
           it_lowModeEval+=2;
           it_evec++;
       }
    }

    // High modes
    if (Nh_ > 0) {

        auto        &v         = envGet(std::vector<FermionField>, getName() + "_v");
        auto        &w         = envGet(std::vector<FermionField>, getName() + "_w");
        auto &noise = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);

        int nsrc = noise.size()/FImpl::Dimension;  

        // Normalization for the noise sources
        RealD norm = 1.0/::sqrt(Real(nsrc));

        std::cout << "Normalizing stochastic vectors by 1/sqrt(" << nsrc << ")" << std::endl;

        FermionField *multiRHSource, *multiRHSolve;
        if (usesMultiRHS) {
            multiRHSource = env().template getObject<FermionField>(getName() + "_tmp_multiRHSource");
            multiRHSolve = env().template getObject<FermionField>(getName() + "_tmp_multiRHSolve");
        }

        for (int ih = 0; ih < Nh_; ih++)
        {
           startTimer("W high mode");
           LOG(Message) << "W vector i = " << Nl_ + ih
                        << " (" << ((hasEpack_) ? "high " : "") 
                        << "stochastic mode)" << std::endl;
            if (hasEpack_ || !par().lowModes.empty()) {
                FermionField *lowModeVecs;
                if (hasEpack_)
                    lowModeVecs = &(envGet(std::vector<FermionField>, getName() + "_evec"));
                else
                    lowModeVecs = &(envGet(std::vector<FermionField>, par().lowModes + "_evec"));
                a2a.makeHighModeW(w[ih], noise.getFerm(ih),lowModeVecs,lowModeVecs.size());
            } else {
                a2a.makeHighModeW(w[ih], noise.getFerm(ih));
            }

            w[ih] = norm*w[ih];
            if (usesMultiRHS) {
                InsertSlice(w[ih],*multiRHSource,ih,0);
            }

           stopTimer("W high mode");
           if (!usesMultiRHS) {
               startTimer("V high mode");
               LOG(Message) << "V vector i = " << Nl_ + ih
                            << " (" << ((hasEpack_) ? "high " : "") 
                            << "stochastic mode)" << std::endl;

               a2a.makeHighModeV(v[ih],w[ih]);

               stopTimer("V high mode");
           }
        }

        if (usesMultiRHS) {

            a2a.makeHighModeV(*multiRHSolve,*multiRHSource);
            for (int ih = 0; ih < Nh_; ih++) {
                ExtractSlice(v[ih],*multiRHSolve,ih,0);
            }
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
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectorsMILC_hpp_
