/*
 * A2AVectors.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
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
    unsigned int Nl_{0};
    bool hasEpack_{false};
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
    std::vector<std::string> out = {getName() + "_v", getName() + "_w", getName()+"_lowModes_evec", getName()+"_lowModes_eval"};

    return out;
}

/******************************************************************************
 *              TA2AVectorsMILC setup                                         *
 ******************************************************************************/
template <typename FImpl, typename Pack>
void TA2AVectorsMILC<FImpl, Pack>::setup(void)
{
    auto        &noise      = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto        &action     = envGet(FMat, par().action);
    auto        &solver     = envGet(Solver, par().solver);
    int         Ls          = env().getObjectLs(par().action);

    envTmp(A2A, "a2a", 1, action, solver);
    hasEpack_ = !par().eigenPack.empty();
    if (hasEpack_) {

        auto &epack = envGet(Pack, par().eigenPack);
        Nl_ = epack.evec.size()*(IsStaggeredImpl<FImpl>()?2:1);
        envCreate(std::vector<FermionField>, getName() + "_lowModes_evec", Ls, Nl_, envGetGrid(FermionField, Ls));
        envCreate(std::vector<ComplexD>, getName() + "_lowModes_eval", Ls, Nl_, 0);
    }

    envCreate(std::vector<FermionField>, getName() + "_v", 1, 
              noise.fermSize(), envGetGrid(FermionField));
    envCreate(std::vector<FermionField>, getName() + "_w", 1, 
              noise.fermSize(), envGetGrid(FermionField));

    if (Ls > 1) {
       HADRONS_ERROR(Argument, "Ls > 1 not implemented");
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
       auto        &noise     = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
       auto        &v         = envGet(std::vector<FermionField>, getName() + "_v");
       auto        &w         = envGet(std::vector<FermionField>, getName() + "_w");
//       int         Ls         = env().getObjectLs(par().action);
       Real        mass;
       envGetTmp(A2A, a2a);

       if (Nl_ > 0)
       {
           LOG(Message) << "Computing all-to-all vectors "
                        << " using eigenpack '" << par().eigenPack << "' ("
                        << Nl_ << " low modes) and noise '"
                        << par().noise << "' (" << noise.fermSize() 
                        << " noise vectors)" << std::endl;

           LOG(Message) << "Eigenpack with conjugate pair evecs and corresponding evals available in '" 
                        << getName() << "_lowModes'" << std::endl;
       } else {

           LOG(Message) << "Computing all-to-all vectors "
                        << " using noise '" << par().noise << "' (" << noise.fermSize() 
                        << " noise vectors)" << std::endl;
       }

       typename std::vector<FermionField>::iterator it_evec, it_lowModeEvec;
       typename std::vector<Real>::iterator it_eval;
       typename std::vector<ComplexD>::iterator it_lowModeEval;

       if (Nl_ > 0) {

           auto &lowModeVecs = envGet(std::vector<FermionField>, getName() + "_lowModes_evec");
           auto &lowModeVals = envGet(std::vector<ComplexD>, getName() + "_lowModes_eval");
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
       int nsrc = noise.size();  

       if (nsrc > 0) {
           // Normalization for the noise sources
           RealD norm = 1.0/::sqrt(Real(nsrc));
           
           std::cout << "Normalizing stochastic vectors by 1/sqrt(" << nsrc << ")" << std::endl;

           for (unsigned int ih = 0; ih < noise.fermSize(); ih++)
           {
               startTimer("W high mode");
               LOG(Message) << "W vector i = " << Nl_ + ih
                            << " (" << ((hasEpack_) ? "high " : "") 
                            << "stochastic mode)" << std::endl;
                if (hasEpack_) {
                    auto &lowModeVecs = envGet(std::vector<FermionField>, getName() + "_lowModes_evec");
                    a2a.makeHighModeW(w[ih], noise.getFerm(ih),lowModeVecs,lowModeVecs.size());
                } else {
                    a2a.makeHighModeW(w[ih], noise.getFerm(ih));
                }

                w[ih] = norm*w[ih];

               stopTimer("W high mode");
               startTimer("V high mode");
               LOG(Message) << "V vector i = " << Nl_ + ih
                            << " (" << ((hasEpack_) ? "high " : "") 
                            << "stochastic mode)" << std::endl;

               a2a.makeHighModeV(v[ih],w[ih]);

               stopTimer("V high mode");
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

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectorsMILC_hpp_
