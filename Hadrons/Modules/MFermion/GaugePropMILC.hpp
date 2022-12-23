/*
 * GaugePropMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Guido Cossu <guido.cossu@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Nils Asmussen <n.asmussen@soton.ac.uk>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: pretidav <david.preti@csic.es>
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

#ifndef Hadrons_MFermion_GaugePropMILC_hpp_
#define Hadrons_MFermion_GaugePropMILC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                GaugePropMILC                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;

class GaugePropMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugePropMILCPar,
                                    std::string, source,
                                    std::string, gammas,
                                    std::string, solver,
                                    std::string, guess);
};

template <typename FImpl>
class TGaugePropMILC: public Module<GaugePropMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);

public:
    // constructor
    TGaugePropMILC(const std::string name);
    // destructor
    virtual ~TGaugePropMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    template <typename TField>
    void setupHelper(void);
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void solvePropagator(FermionField &prop, const FermionField &src, const FermionField* guess = nullptr);
    void solvePropagator(PropagatorField &prop, const PropagatorField &src, const PropagatorField* guess = nullptr);

    template <typename TField>
    void solvePropagator(std::map<StagGamma,TField> &prop, const TField &src);
    template <typename TField>
    void solvePropagator(std::vector<TField> &prop, const std::vector<TField> &src);
    template <typename TField>
    void solvePropagator(std::map<StagGamma,std::vector<TField>> &prop, const std::vector<TField> &src);
private:
    bool hasGammas_;
};

MODULE_REGISTER_TMP(StagGaugeProp, TGaugePropMILC<STAGIMPL>, MFermion);

/******************************************************************************
 *                      TGaugePropMILC implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TGaugePropMILC<FImpl>::TGaugePropMILC(const std::string name)
: Module<GaugePropMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TGaugePropMILC<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().solver};
    hasGammas_ = !par().gammas.empty();

    if (!par().guess.empty()) {
        in.push_back(par().guess);
    }
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TGaugePropMILC<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
template <typename TField>
void TGaugePropMILC<FImpl>::setupHelper() {
    envTmpLat(TField, "field");
    envTmpLat(FermionField, "fermIn");
    envTmpLat(FermionField, "fermOut");
    envTmpLat(FermionField, "fermGuess");

    envGetTmp(FermionField, fermIn);
    envGetTmp(FermionField, fermOut);
    envGetTmp(FermionField, fermGuess);
    fermIn = Zero();
    fermOut = Zero();
    fermGuess = Zero();

    if (envHasType(TField,par().source)) {
        if (hasGammas_) {
            envGetTmp(std::vector<GammaPair>,gammaList);

            std::map<Gamma::Algebra,TField> dummy;
            envCreate(ARG(std::map<Gamma::Algebra,TField>), getName(), 1, dummy);
            auto &sol = envGet(ARG(std::map<Gamma::Algebra,TField>), getName());
            for (auto & gamma:gammaList) {
                sol.insert({gamma,envGetGrid(TField)});
                sol.at(gamma) = Zero();
            }
        } else {
            envCreateLat(TField, getName());
        }
    } else {
        int srcSize = 0;
        if (envHasType(ARG(std::map<Gamma::Algebra,std::vector<TField> >),par().source)) {
            auto &src = envGet(ARG(std::map<Gamma::Algebra,std::vector<TField> >), par().source);
            srcSize = src.at(Gamma::Algebra::Gamma5).size();
        } else {
            auto &src = envGet(std::vector<TField>, par().source);
            srcSize = src.size();
        }

        if (hasGammas_) {
            envGetTmp(std::vector<Gamma::Algebra>,gammaList);
            std::map<Gamma::Algebra,std::vector<TField>> dummy;
            envCreate(ARG(std::map<Gamma::Algebra,std::vector<TField>>), getName(), 1, dummy);
            auto &sol = envGet(ARG(std::map<Gamma::Algebra,std::vector<TField>>), getName());
            for (auto & gamma:gammaList) {
                sol.insert({gamma,std::vector<TField>(srcSize,envGetGrid(TField))});
                for (auto & s:sol.at(gamma)) {
                    s = Zero();
                }
            }
        } else {
            envCreate(std::vector<TField>, getName(), 1, srcSize,
                      envGetGrid(TField));
        }
    }
}

template <typename FImpl>
void TGaugePropMILC<FImpl>::setup(void)
{
    envTmp(StagGamma,"spinTaste",1,0,0);
   envTmp(std::vector<GammaPair>,"gammaList",1,0);

    envGetTmp(std::vector<GammaPair>,gammaList);
    gammaList.clear();

    if (hasGammas_)  {
        gammaList = strToVec<GammaPair>(par().gammas);
    }
    
    if (envHasType(PropagatorField,par().source) || envHasType(std::vector<PropagatorField>,par().source) || envHasType(ARG(std::map<StagGamma,std::vector<PropagatorField> >),par().source)) {
        setupHelper<PropagatorField>();
    } else if (envHasType(FermionField,par().source) || envHasType(std::vector<FermionField>,par().source)|| envHasType(ARG(std::map<StagGamma,std::vector<FermionField> >),par().source)) {
        setupHelper<FermionField>();
    } else {
        HADRONS_ERROR(Logic,"Type of source '" + par().source + "' not recognized.");
    }

}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaugePropMILC<FImpl>::solvePropagator(FermionField &sol, const FermionField &src, const FermionField *guess)
{
    auto &solver  = envGet(Solver, par().solver);
    
    if (guess != nullptr) {
        solver(sol, src,*guess);
    } else {
        solver(sol, src);
    }
}

template <typename FImpl>
void TGaugePropMILC<FImpl>::solvePropagator(PropagatorField &sol, 
                                            const PropagatorField &src, const PropagatorField *guess)
{
    auto &solver  = envGet(Solver, par().solver);
    
    envGetTmp(FermionField, fermIn);
    envGetTmp(FermionField, fermOut);
    envGetTmp(FermionField, fermGuess);

    for (unsigned int c = 0; c < FImpl::Dimension; ++c)
    {
        PropToFerm<FImpl>(fermIn, src, c);
        if (guess != nullptr) {
            PropToFerm<FImpl>(fermGuess,*guess,c);
            solver(fermOut, fermIn,fermGuess);
        } else {
            solver(fermOut, fermIn);
        }
        FermToProp<FImpl>(sol, fermOut, c);
    }
}

template <typename FImpl>
template<typename TField>
void TGaugePropMILC<FImpl>::solvePropagator(std::vector<TField> &sol, const std::vector<TField> &src)
{
    for (int i = 0;i<src.size();i++) {
        LOG(Message) << "Solving element " << i << " of '" << par().source << "'" << std::endl;
        solvePropagator(sol[i],src[i]);
    }
}

template <typename FImpl>
template<typename TField>
void TGaugePropMILC<FImpl>::solvePropagator(std::map<StagGamma,TField> &sol, const TField &src)
{
    envGetTmp(std::vector<GammaPair>,gammaList);
    envGetTmp(TField,field);

    std::map<StagGamma,TField> *guess;
    if (!par().guess.empty()) {
        if (!envHasType(ARG(std::map<StagGamma,TField>),par().guess)) {
            HADRONS_ERROR(Argument, "guess parameter '" + par().guess + "' must have same data structure as source, '"+par().source+"'");
        }
        guess = env().getObject<std::map<StagGamma,TField>>(par().guess);
    }


    envGetTmp(StagGamma,spinTaste);

    for (auto &gamma:gammaList) {
        spinTaste.g_spin = gamma.first;
        spinTaste.g_taste = gamma.second;
        LOG(Message) << "Solve for '" << par().source << "' with '(" <<
            Gamma::name[gamma.first] << ", " << Gamma::name[gamma.second] << ")'" << std::endl;


        field = spinTaste*src;

        if (!par().guess.empty()) {
            const TField& guessField = guess->at(spinTaste);
            solvePropagator(sol.at(spinTaste),field,&guessField);
        } else {
            solvePropagator(sol.at(spinTaste),field);
        }

    }
}

template <typename FImpl>
template<typename TField>
void TGaugePropMILC<FImpl>::solvePropagator(std::map<StagGamma,std::vector<TField>> &sol, const std::vector<TField> &src)
{
    envGetTmp(std::vector<GammaPair>,gammaList);
    envGetTmp(TField,field);

    std::map<Gamma::Algebra,std::vector<TField> > *guess;
    if (!par().guess.empty()) {
        if (!envHasType(ARG(std::map<Gamma::Algebra,std::vector<TField> >),par().guess)) {
            HADRONS_ERROR(Argument, "guess parameter '" + par().guess + "' must have same data structure as source, '"+par().source+"'");
        }
        guess = env().getObject<std::map<Gamma::Algebra,std::vector<TField> >>(par().guess);
    }

    envGetTmp(StagGamma,spinTaste);

    for (auto &gamma:gammaList) {
        spinTaste.g_spin = gamma.first;
        spinTaste.g_taste = gamma.second;
        LOG(Message) << "Solve for '" << par().source << "' with '(" << 
        Gamma::name[gamma.first] << ", " << Gamma::name[gamma.second] << ")'" << std::endl;

        for (int i = 0;i<src.size();i++) {
            LOG(Message) << "Solving element " << i << " of '" << par().source << "'" << std::endl;
            const TField& srctmp = src[i];

            field = spinTaste*srctmp;

            if (!par().guess.empty()) {
                const TField& guessField = guess->at(spinTaste)[i];
                solvePropagator(sol.at(spinTaste)[i],field,&guessField);
            } else {
                solvePropagator(sol.at(spinTaste)[i],field);
            }
        }
    }
}

template <typename FImpl>
void TGaugePropMILC<FImpl>::execute(void)
{
    LOG(Message) << "Computing quark propagator '" << getName() << "'"
                 << std::endl;
    
    if (envHasType(PropagatorField,par().source)) {
        auto &src = envGet(PropagatorField,par().source);
        if (hasGammas_) {
            auto &sol = envGet(ARG(std::map<StagGamma,PropagatorField>),getName());
            solvePropagator(sol,src);
        } else {
            auto &sol = envGet(PropagatorField,getName());
            solvePropagator(sol,src);
        }
    } else if (envHasType(std::vector<PropagatorField>, par().source) || envHasType(ARG(std::map<StagGamma,std::vector<PropagatorField> >), par().source)) {
        std::vector<PropagatorField> *src;
        if (envHasType(std::vector<PropagatorField>, par().source)) {
            auto &srctmp = envGet(std::vector<PropagatorField>, par().source);
            src = &srctmp;
        } else {
            auto &srctmp = envGet(ARG(std::map<StagGamma,std::vector<PropagatorField> >), par().source);
            src = &(srctmp.at(StagGamma(Gamma::Algebra::Gamma5,Gamma::Algebra::Gamma5)));
        }
        if (hasGammas_) {
            auto &sol = envGet(ARG(std::map<StagGamma,std::vector<PropagatorField> >),getName());
            solvePropagator(sol,*src);
        } else {
            auto &sol = envGet(std::vector<PropagatorField>,getName());
            solvePropagator(sol,*src);
        }
    } else if (envHasType(FermionField,par().source)) {
        auto &src = envGet(FermionField,par().source);
        if (hasGammas_) {
            auto &sol = envGet(ARG(std::map<StagGamma,FermionField>),getName());
            solvePropagator(sol,src);
        } else {
            auto &sol = envGet(FermionField,getName());
            solvePropagator(sol,src);
        }
    } else {
        std::vector<FermionField> *src;
        if (envHasType(std::vector<FermionField>, par().source)) {
            auto &srctmp = envGet(std::vector<FermionField>, par().source);
            src = &srctmp;
        } else {
            auto &srctmp = envGet(ARG(std::map<StagGamma,std::vector<FermionField> >), par().source);
            src = &(srctmp.at(StagGamma(Gamma::Algebra::Gamma5,Gamma::Algebra::Gamma5)));
        }
        if (hasGammas_) {
            auto &sol = envGet(ARG(std::map<StagGamma,std::vector<FermionField>>),getName());
            solvePropagator(sol,*src);
        } else {
            auto &sol = envGet(std::vector<FermionField>,getName());
            solvePropagator(sol,*src);
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_GaugePropMILC_hpp_
