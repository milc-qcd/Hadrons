/*
 * RBPrecCGMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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

#ifndef Hadrons_MSolver_RBPrecCGMILC_hpp_
#define Hadrons_MSolver_RBPrecCGMILC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MSolver/GuesserMILC.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     Schur red-black preconditioned CG                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class RBPrecCGMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RBPrecCGMILCPar ,
                                    std::string , action,
                                    unsigned int, maxIteration,
                                    double      , residual,
                                    std::string , eigenPack,
                                    std::string , mustConverge,
                                    std::string , evenEigen);
};

template <typename FImpl, int nBasis = HADRONS_DEFAULT_LANCZOS_NBASIS>
class TRBPrecCGMILC: public Module<RBPrecCGMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    HADRONS_DEFINE_SCHUR_SOLVE(schurSolve_t,FImpl);
    HADRONS_DEFINE_SCHUR_OP(schurOp_t,FImpl);
public:
    // constructor
    TRBPrecCGMILC(const std::string name);
    // destructor
    virtual ~TRBPrecCGMILC(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
public:
    bool mustConverge_ = true;
    bool evenEigen_    = false;
};

MODULE_REGISTER_TMP(StagRBPrecCG, ARG(TRBPrecCGMILC<STAGIMPL>), MSolver);

/******************************************************************************
 *                      TRBPrecCGMILC template implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
TRBPrecCGMILC<FImpl, nBasis>::TRBPrecCGMILC(const std::string name)
: Module(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
std::vector<std::string> TRBPrecCGMILC<FImpl, nBasis>::getInput(void)
{
    std::vector<std::string> in = {par().action};
    
    if (!par().eigenPack.empty())
    {
        in.push_back(par().eigenPack);
    }
    
    return in;
}

template <typename FImpl, int nBasis>
DependencyMap TRBPrecCGMILC<FImpl, nBasis>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().action, getName()});
    dep.insert({par().action, getName() + "_subtract"});

    if (!par().eigenPack.empty())
    {
        dep.insert({par().eigenPack, getName(),             });
        dep.insert({par().eigenPack, getName() + "_subtract"});

        if (env().hasObject(par().eigenPack + "_eval")) {
            dep.insert({par().eigenPack+"_eval", getName(),             });
            dep.insert({par().eigenPack+"_eval", getName() + "_subtract"});
        }
    }

    return dep;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TRBPrecCGMILC<FImpl, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_subtract"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TRBPrecCGMILC<FImpl, nBasis>::setup(void)
{
    if (par().maxIteration == 0)
    {
        HADRONS_ERROR(Argument, "zero maximum iteration");
    }

    LOG(Message) << "setting up Schur red-black preconditioned CG for"
                 << " action '" << par().action << "' with residual "
                 << par().residual << ", maximum iteration " 
                 << par().maxIteration << std::endl;

    auto Ls        = env().getObjectLs(par().action);
    auto &mat      = envGet(FMat, par().action);
    auto guesserPt = makeGuesser<FImpl, nBasis>(par().eigenPack);

    if (!par().mustConverge.empty()) {
        if (!(std::istringstream(par().mustConverge) >> std::boolalpha >> mustConverge_)) {
            HADRONS_ERROR(Logic,"parameter mustConverge='" + par().mustConverge + "' must be 'true' or 'false'");
        }
    }

    if (!par().evenEigen.empty()) {
        if (!(std::istringstream(par().evenEigen) >> std::boolalpha >> evenEigen_)) {
            HADRONS_ERROR(Logic,"parameter evenEigen='" + par().evenEigen + "' must be 'true' or 'false'");
        }
    }

    auto makeSolver = [&mat, guesserPt, this](bool subGuess) {
        return [&mat, guesserPt, subGuess, this](FermionField &sol,
                                     const FermionField &source) {

            if (subGuess) {
                HADRONS_ERROR(Implementation,"Guess subtraction not supported for solver '"+getName()+"'.");
            }

            GridBase *g = envGetRbGrid(FermionField);
            FermionField tmp(envGetGrid(FermionField));

            MdagMLinearOperator<FMat, FermionField> hermOp(mat);
            ConjugateGradient<FermionField> cg(par().residual,
                                               par().maxIteration,
                                               this->mustConverge_);
            schurSolve_t<FermionField> schurSolver(cg,false,false,(this->evenEigen_?Even:Odd));
            // schurSolver.subtractGuess(subGuess);

            // Perform initial solve
            schurSolver(mat, source, sol, *guesserPt);

            LOG(Message) << "Improving residual of full field." << std::endl;

            // Improve residual of full field to meet desired result
            mat.Mdag(source, tmp);
            cg(hermOp, tmp, sol);
        };
    };
    auto solver = makeSolver(false);
    envCreate(Solver, getName(), Ls, solver, mat);
    auto solver_subtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls, solver_subtract, mat);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TRBPrecCGMILC<FImpl, nBasis>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_RBPrecCGMILC_hpp_
