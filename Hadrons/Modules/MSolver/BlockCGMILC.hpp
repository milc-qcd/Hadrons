/*
 * BlockCGMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MSolver_BlockCGMILC_hpp_
#define Hadrons_MSolver_BlockCGMILC_hpp_

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

class BlockCGMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(BlockCGMILCPar ,
                                    std::string , action5D,
                                    unsigned int, maxIteration,
                                    double      , residual,
                                    std::string , eigenPack);
};

template <typename FImpl, int nBasis = HADRONS_DEFAULT_LANCZOS_NBASIS>
class TBlockCGMILC: public Module<BlockCGMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    HADRONS_DEFINE_SCHUR_SOLVE(schurSolve_t,FImpl);
public:
    // constructor
    TBlockCGMILC(const std::string name);
    // destructor
    virtual ~TBlockCGMILC(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getReference(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(BlockCG, ARG(TBlockCGMILC<STAGIMPL>), MSolver);

/******************************************************************************
 *                      TBlockCGMILC template implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
TBlockCGMILC<FImpl, nBasis>::TBlockCGMILC(const std::string name)
: Module(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
std::vector<std::string> TBlockCGMILC<FImpl, nBasis>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TBlockCGMILC<FImpl, nBasis>::getReference(void)
{
    std::vector<std::string> ref = {par().action5D};
    
    if (!par().eigenPack.empty())
    {
        ref.push_back(par().eigenPack);
        ref.push_back(par().eigenPack+"_evenEigen");
    }

    return ref;
}

template <typename FImpl, int nBasis>
std::vector<std::string> TBlockCGMILC<FImpl, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_subtract"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TBlockCGMILC<FImpl, nBasis>::setup(void)
{
    if (par().maxIteration == 0)
    {
        HADRONS_ERROR(Argument, "zero maximum iteration");
    }

    LOG(Message) << "setting up Schur red-black preconditioned CG for"
                 << " action '" << par().action5D << "' with residual "
                 << par().residual << ", maximum iteration " 
                 << par().maxIteration << std::endl;

    auto Ls        = env().getObjectLs(par().action5D);
    auto &mat      = envGet(FMat, par().action5D);

    if (Ls == 1) {
        HADRONS_ERROR(Argument, "Action must be 5D with Ls equal to number of right hand solves")
    }

    auto guesserPt = makeGuesser<FImpl, nBasis>(par().eigenPack);

    int checkerboard = Odd;
    if (!par().eigenPack.empty())
        checkerboard = (envGet(std::vector<bool>, par().eigenPack+"_evenEigen"))[0] ? Even : Odd;

    int blockDim = 0;
    auto makeSolver = [&mat, guesserPt, checkerboard,blockDim, this](bool subGuess) {
        return [&mat, guesserPt, subGuess, checkerboard,blockDim, this](FermionField &sol,
                                     const FermionField &source) {
            BlockConjugateGradient<FermionField> bcg(CGmultiRHS,blockDim,(RealD)par().residual,
                                               (Integer)par().maxIteration);
            schurSolve_t<FermionField> schurSolver(bcg,false,false,checkerboard);
            schurSolver.subtractGuess(subGuess);
            schurSolver(mat, source, sol, *guesserPt);
        };
    };
    auto solver = makeSolver(false);
    envCreate(Solver, getName(), Ls, solver, mat);
    auto solver_subtract = makeSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls, solver_subtract, mat);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TBlockCGMILC<FImpl, nBasis>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_BlockCGMILC_hpp_
