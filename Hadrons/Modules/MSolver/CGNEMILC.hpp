#ifndef Hadrons_MSolver_CGNEMILC_hpp_
#define Hadrons_MSolver_CGNEMILC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/Modules/MSolver/GuesserMILC.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         CGNEMILC                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class CGNEMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CGNEMILCPar,
                                    std::string , action,
                                    unsigned int, maxIteration,
                                    double      , residual,
                                    std::string , eigenPack,
                                    std::string , mustConverge);
};

template <typename FImpl, int nBasis = HADRONS_DEFAULT_LANCZOS_NBASIS>
class TCGNEMILC: public Module<CGNEMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TCGNEMILC(const std::string name);
    // destructor
    virtual ~TCGNEMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
public:
    bool mustConverge_ = true;
};

MODULE_REGISTER_TMP(StagCGNE, TCGNEMILC<STAGIMPL>, MSolver);

/******************************************************************************
 *                 TCGNEMILC implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
TCGNEMILC<FImpl, nBasis>::TCGNEMILC(const std::string name)
: Module<CGNEMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
std::vector<std::string> TCGNEMILC<FImpl, nBasis>::getInput(void)
{
    std::vector<std::string> in = {par().action};
    
    if (!par().eigenPack.empty())
    {
        in.push_back(par().eigenPack);
    }
    
    return in;
}

template <typename FImpl, int nBasis>
DependencyMap TCGNEMILC<FImpl, nBasis>::getObjectDependencies(void)
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
std::vector<std::string> TCGNEMILC<FImpl, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {getName(),getName()+"_subtract"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TCGNEMILC<FImpl, nBasis>::setup(void)
{
    if (par().maxIteration == 0)
    {
        HADRONS_ERROR(Argument, "zero maximum iteration");
    }

    LOG(Message) << "setting up normal equation CG for"
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

    auto makeSolver = [&mat, guesserPt, this](bool subGuess) 
    {
        return [&mat, guesserPt, subGuess, this](FermionField &sol,
                                                 const FermionField &source) 
        {
            GridBase                                *g = sol.Grid();
            FermionField                            guess(g), tmp(g);
            MdagMLinearOperator<FMat, FermionField> hermOp(mat);
            ConjugateGradient<FermionField>         cg(par().residual,
                                                       par().maxIteration,
                                                       this->mustConverge_);

            guess = sol;
            mat.Mdag(source, tmp);
            (*guesserPt)(tmp, sol);
            cg(hermOp, tmp, sol);
            if (subGuess)
            {
                sol -= guess;
            }
        };
    };

    auto makeGuessSolver = [&mat, this](bool subGuess) 
    {
        return [&mat, subGuess, this](FermionField &sol,
                                                 const FermionField &source, const FermionField &guess) 
        {
            GridBase                                *g = sol.Grid();
            FermionField                            tmp(g);
            MdagMLinearOperator<FMat, FermionField> hermOp(mat);
            ConjugateGradient<FermionField>         cg(par().residual,
                                                       par().maxIteration,
                                                       this->mustConverge_);

            sol = guess;
            mat.Mdag(source, tmp);

            cg(hermOp, tmp, sol);
            if (subGuess)
            {
                sol -= guess;
            }
        };
    };

    auto solver = makeSolver(false);
    auto guessSolver = makeGuessSolver(false);
    envCreate(Solver, getName(), Ls, solver, guessSolver, mat);
    auto solver_subtract = makeSolver(true);
    auto guessSolver_subtract = makeGuessSolver(true);
    envCreate(Solver, getName() + "_subtract", Ls, solver_subtract, guessSolver_subtract, mat);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int nBasis>
void TCGNEMILC<FImpl, nBasis>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_CGNEMILC_hpp_
