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
#ifndef Hadrons_MSolver_A2AVectors_hpp_
#define Hadrons_MSolver_A2AVectors_hpp_

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

class A2AVectorsPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsPar,
                                  std::string, noise,
                                  std::string, action,
                                  double     , mass,
                                  std::string, eigenPack,
                                  std::string, solver,
                                  std::string, output,
                                  bool,        evenEigen,
                                  bool,        multiFile);
};

template <typename FImpl, typename Pack>
class TA2AVectors : public Module<A2AVectorsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef A2AVectorsSchur<FImpl> A2A;

public:
    // constructor
    TA2AVectors(const std::string name);
    // destructor
    virtual ~TA2AVectors(void) {};
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
};

MODULE_REGISTER_TMP(A2AVectors, 
    ARG(TA2AVectors<FIMPL, BaseFermionEigenPack<FIMPL>>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectors, 
    ARG(TA2AVectors<ZFIMPL, BaseFermionEigenPack<ZFIMPL>>), MSolver);
MODULE_REGISTER_TMP(StagA2AVectors, 
    ARG(TA2AVectors<STAGIMPL, BaseFermionEigenPack<STAGIMPL>>), MSolver);

/******************************************************************************
 *                       TA2AVectors implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TA2AVectors<FImpl, Pack>::TA2AVectors(const std::string name)
: Module<A2AVectorsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectors<FImpl, Pack>::getInput(void)
{
    std::string              sub_string;
    std::vector<std::string> in;

    if (!par().eigenPack.empty())
    {
        in.push_back(par().eigenPack);
        sub_string = (!par().eigenPack.empty()) ? "_subtract" : "";
    }
    in.push_back(par().solver + sub_string);
    in.push_back(par().noise);

    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TA2AVectors<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "_v", getName() + "_w"};

    return out;
}

/******************************************************************************
 *              TA2AVectors setup   (NON-STAGGERED)                           *
 ******************************************************************************/
template <typename FImpl, typename Pack>
void TA2AVectors<FImpl, Pack>::setup(void)
{
    bool        hasLowModes = (!par().eigenPack.empty());
    std::string sub_string  = (hasLowModes) ? "_subtract" : "";
    auto        &noise      = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto        &action     = envGet(FMat, par().action);
    auto        &solver     = envGet(Solver, par().solver + sub_string);
    int         Ls          = env().getObjectLs(par().action);

    envTmp(A2A, "a2a", 1, action, solver);
    envGetTmp(A2A, a2a);

    if (hasLowModes)
    {
        auto &epack = envGet(Pack, par().eigenPack);
        Nl_ = epack.evec.size()*(a2a.isStaggered()?2:1);
    }
    envCreate(std::vector<FermionField>, getName() + "_v", 1, 
              Nl_ + noise.fermSize(), envGetGrid(FermionField));
    envCreate(std::vector<FermionField>, getName() + "_w", 1, 
              Nl_ + noise.fermSize(), envGetGrid(FermionField));
    if (Ls > 1)
    {
        if (a2a.isStaggered()) {
            envTmp(std::vector<FermionField>, "f5", Ls, 2, envGetGrid(FermionField,Ls));
            envTmp(std::vector<FermionField>, "f5_2", Ls, 2, envGetGrid(FermionField,Ls));
        } else 
            envTmpLat(FermionField, "f5", Ls);

    }
}

/******************************************************************************
 *              TA2AVectors execution   (NON-STAGGERED)                       *
 ******************************************************************************/
template <typename FImpl, typename Pack>
void TA2AVectors<FImpl, Pack>::execute(void)
{
    std::string sub_string = (Nl_ > 0) ? "_subtract" : "";
       auto        &action    = envGet(FMat, par().action);
       auto        &solver    = envGet(Solver, par().solver + sub_string);
       auto        &noise     = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
       auto        &v         = envGet(std::vector<FermionField>, getName() + "_v");
       auto        &w         = envGet(std::vector<FermionField>, getName() + "_w");
       int         Ls         = env().getObjectLs(par().action);

       envGetTmp(A2A, a2a);

       if (Nl_ > 0)
       {
           LOG(Message) << "Computing all-to-all vectors "
                        << " using eigenpack '" << par().eigenPack << "' ("
                        << Nl_ << " low modes) and noise '"
                        << par().noise << "' (" << noise.fermSize() 
                        << " noise vectors)" << std::endl;
       }
       else
       {
           LOG(Message) << "Computing all-to-all vectors "
                        << " using noise '" << par().noise << "' (" << noise.fermSize() 
                        << " noise vectors)" << std::endl;
       }

       typename std::vector<FermionField>::iterator it_w, it_v, it_evec;
       typename std::vector<Real>::iterator it_eval;

       if (!par().eigenPack.empty()) {

           auto &epack  = envGet(Pack, par().eigenPack);
           it_w = w.begin();
           it_v = v.begin();
           it_evec = epack.evec.begin();

           // Low modes
           for (auto it_eval = epack.eval.begin(); it_eval < epack.eval.end(); it_eval++)
           {
               int il = it_eval-epack.eval.begin();

               if(a2a.isStaggered()) {
                   startTimer("low mode pair");
                   LOG(Message) << "V,W vector pairs for i = " << il << " and " << il+1 << " (low mode)" << std::endl;
                   if (Ls == 1)
                       a2a.makeLowModePairs(it_v, it_w, it_evec, par().mass, *it_eval, par().evenEigen == true);
                   else {
                       envGetTmp(std::vector<FermionField>, f5);
                       envGetTmp(std::vector<FermionField>, f5_2);
                       typename std::vector<FermionField>::iterator it_f5 = f5.begin(), it_f5_2 = f5_2.begin();
                       a2a.makeLowModePairs5D(it_v, it_f5, it_w, it_f5_2, it_evec, par().mass, *it_eval, par().evenEigen == true); 
                   }
                   stopTimer("low mode pair");
               } else {
                   if (Ls == 1) {
                       startTimer("V low mode");
                       LOG(Message) << "V vector i = " << il << " (low mode)" << std::endl;
           
                       a2a.makeLowModeV(*it_v, *it_evec, *it_eval);
                   
                       stopTimer("V low mode");
                       startTimer("W low mode");
                       LOG(Message) << "W vector i = " << il << " (low mode)" << std::endl;

                       a2a.makeLowModeW(*it_w, *it_evec, *it_eval); 
                   } else {
                       startTimer("V low mode");
                       LOG(Message) << "V vector i = " << il << " (low mode)" << std::endl;

                       envGetTmp(FermionField, f5);
                       a2a.makeLowModeV5D(*it_v, f5, *it_evec, *it_eval);

                       stopTimer("V low mode");
                       startTimer("W low mode");
                       LOG(Message) << "W vector i = " << il << " (low mode)" << std::endl;

                       a2a.makeLowModeW5D(*it_w, f5, *it_evec, *it_eval);
                   }
                   stopTimer("W low mode");
               }

               if (a2a.isStaggered()) {
                   it_w+=2;
                   it_v+=2;
               } else {
                   it_w++;
                   it_v++;
               }

               it_evec++;
           }
       }

       // High modes
       for (unsigned int ih = 0; ih < noise.fermSize(); ih++)
       {
           startTimer("V high mode");
           LOG(Message) << "V vector i = " << Nl_ + ih
                        << " (" << ((Nl_ > 0) ? "high " : "") 
                        << "stochastic mode)" << std::endl;
           if (Ls == 1)
           {
               a2a.makeHighModeV(v[Nl_ + ih], noise.getFerm(ih));
           }
           else
           {
               envGetTmp(FermionField, f5);
               a2a.makeHighModeV5D(v[Nl_ + ih], f5, noise.getFerm(ih));
           }
           stopTimer("V high mode");
           startTimer("W high mode");
           LOG(Message) << "W vector i = " << Nl_ + ih
                        << " (" << ((Nl_ > 0) ? "high " : "") 
                        << "stochastic mode)" << std::endl;
           if (Ls == 1)
           {
               a2a.makeHighModeW(w[Nl_ + ih], noise.getFerm(ih));
           }
           else
           {
               envGetTmp(FermionField, f5);
               a2a.makeHighModeW5D(w[Nl_ + ih], f5, noise.getFerm(ih));
           }
           stopTimer("W high mode");
       }

       // I/O if necessary
       if (!par().output.empty())
       {
           startTimer("V I/O");
           A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
           stopTimer("V I/O");
           startTimer("W I/O");
           A2AVectorsIo::write(par().output + "_w", w, par().multiFile, vm().getTrajectory());
           stopTimer("W I/O");
       }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectors_hpp_
