/*
 * StagGaugeProp.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MFermion_StagGaugeProp_hpp_
#define Hadrons_MFermion_StagGaugeProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                StagGaugeProp                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class StagGaugePropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagGaugePropPar,
                                    std::string, source,
                                    std::string, gammas,
                                    std::string, gammaFunc,
                                    std::string, solver);
};

template <typename FImpl>
class TStagGaugeProp: public Module<StagGaugePropPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);

    typedef std::function<LatticeComplex (Gamma::Algebra gamma)> GammaFn;
public:
    // constructor
    TStagGaugeProp(const std::string name);
    // destructor
    virtual ~TStagGaugeProp(void) {};
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
    void solvePropagator(FermionField &prop, const FermionField &src);
    void solvePropagator(PropagatorField &prop, const PropagatorField &src);

    template <typename TField>
    void solvePropagator(std::map<Gamma::Algebra,TField> &prop, const TField &src);
    template <typename TField>
    void solvePropagator(std::vector<TField> &prop, const std::vector<TField> &src);
    template <typename TField>
    void solvePropagator(std::map<Gamma::Algebra,std::vector<TField>> &prop, const std::vector<TField> &src);
private:
    bool hasGammas_;
};

MODULE_REGISTER_TMP(StagGaugeProp, TStagGaugeProp<STAGIMPL>, MFermion);

/******************************************************************************
 *                      TStagGaugeProp implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStagGaugeProp<FImpl>::TStagGaugeProp(const std::string name)
: Module<StagGaugePropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStagGaugeProp<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().solver};
    hasGammas_ = !par().gammas.empty();

    if (hasGammas_) {
        in.push_back(par().gammaFunc);
    }
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TStagGaugeProp<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
template <typename TField>
void TStagGaugeProp<FImpl>::setupHelper() {
    envTmpLat(TField, "field");
    envTmpLat(FermionField, "ferm1");
    envTmpLat(FermionField, "ferm2");
    envTmpLat(LatticeComplex,"stagPhase");

    envGetTmp(FermionField, ferm1);
    envGetTmp(FermionField, ferm2);
    ferm1 = Zero();
    ferm2 = Zero();

    if (envHasType(TField,par().source)) {
        if (hasGammas_) {
            envGetTmp(std::vector<Gamma::Algebra>,gammaList);

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
void TStagGaugeProp<FImpl>::setup(void)
{
    envTmp(std::vector<Gamma::Algebra>,"gammaList",1,0);

    envGetTmp(std::vector<Gamma::Algebra>,gammaList);
    gammaList.clear();

    if (hasGammas_)  {
        gammaList = strToVec<Gamma::Algebra>(par().gammas);
    }
    
    if (envHasType(PropagatorField,par().source) || envHasType(std::vector<PropagatorField>,par().source) || envHasType(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField> >),par().source)) {
        setupHelper<PropagatorField>();
    } else if (envHasType(FermionField,par().source) || envHasType(std::vector<FermionField>,par().source)|| envHasType(ARG(std::map<Gamma::Algebra,std::vector<FermionField> >),par().source)) {
        setupHelper<FermionField>();
    } else {
        HADRONS_ERROR(Logic,"Type of source '" + par().source + "' not recognized.");
    }

}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagGaugeProp<FImpl>::solvePropagator(FermionField &sol, const FermionField &src)
{
    auto &solver  = envGet(Solver, par().solver);
    
    solver(sol, src);
}

template <typename FImpl>
void TStagGaugeProp<FImpl>::solvePropagator(PropagatorField &sol, 
                                            const PropagatorField &src)
{
    auto &solver  = envGet(Solver, par().solver);
    
    envGetTmp(FermionField, ferm1);
    envGetTmp(FermionField, ferm2);

    for (unsigned int c = 0; c < FImpl::Dimension; ++c)
    {
        PropToFerm<FImpl>(ferm1, src, c);
        solver(ferm2, ferm1);
        FermToProp<FImpl>(sol, ferm2, c);
    }
}

template <typename FImpl>
template<typename TField>
void TStagGaugeProp<FImpl>::solvePropagator(std::vector<TField> &sol, const std::vector<TField> &src)
{
    for (int i = 0;i<src.size();i++) {
        LOG(Message) << "Solving element " << i << " of '" << par().source << "'" << std::endl;
        solvePropagator(sol[i],src[i]);
    }
}

template <typename FImpl>
template<typename TField>
void TStagGaugeProp<FImpl>::solvePropagator(std::map<Gamma::Algebra,TField> &sol, const TField &src)
{
    envGetTmp(std::vector<Gamma::Algebra>,gammaList);
    envGetTmp(TField,field);
    envGetTmp(LatticeComplex,stagPhase);
    auto &func = envGet(GammaFn, par().gammaFunc);

    for (auto &gamma:gammaList) {
        std::string gammaStr = Gamma::name[gamma];
        LOG(Message) << "Solve for '" << par().source << "' with '" << gammaStr << "'" << std::endl;
        stagPhase = func(gamma);
        field = stagPhase*src;
        solvePropagator(sol.at(gamma),field);

    }
}

template <typename FImpl>
template<typename TField>
void TStagGaugeProp<FImpl>::solvePropagator(std::map<Gamma::Algebra,std::vector<TField>> &sol, const std::vector<TField> &src)
{
    envGetTmp(std::vector<Gamma::Algebra>,gammaList);
    envGetTmp(TField,field);
    envGetTmp(LatticeComplex,stagPhase);
    auto &func = envGet(GammaFn, par().gammaFunc);

    for (auto &gamma:gammaList) {
        std::string gammaStr = Gamma::name[gamma];
        LOG(Message) << "Solve for '" << gammaStr << "'" << std::endl;
    
        stagPhase = func(gamma);
        for (int i = 0;i<src.size();i++) {
            LOG(Message) << "Solving element " << i << " of '" << par().source << "'" << std::endl;
            const TField& srctmp = src[i];
            field = stagPhase*srctmp;
            solvePropagator(sol.at(gamma)[i],field);
        }
    }
}

template <typename FImpl>
void TStagGaugeProp<FImpl>::execute(void)
{
    LOG(Message) << "Computing quark propagator '" << getName() << "'"
                 << std::endl;
    
    if (envHasType(PropagatorField,par().source)) {
        auto &src = envGet(PropagatorField,par().source);
        if (hasGammas_) {
            auto &sol = envGet(ARG(std::map<Gamma::Algebra,PropagatorField>),getName());
            solvePropagator(sol,src);
        } else {
            auto &sol = envGet(PropagatorField,getName());
            solvePropagator(sol,src);
        }
    } else if (envHasType(std::vector<PropagatorField>, par().source) || envHasType(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField> >), par().source)) {
        std::vector<PropagatorField> *src;
        if (envHasType(std::vector<PropagatorField>, par().source)) {
            auto &srctmp = envGet(std::vector<PropagatorField>, par().source);
            src = &srctmp;
        } else {
            auto &srctmp = envGet(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField> >), par().source);
            src = &(srctmp.at(Gamma::Algebra::Gamma5));
        }
        if (hasGammas_) {
            auto &sol = envGet(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField> >),getName());
            solvePropagator(sol,*src);
        } else {
            auto &sol = envGet(std::vector<PropagatorField>,getName());
            solvePropagator(sol,*src);
        }
    } else if (envHasType(FermionField,par().source)) {
        auto &src = envGet(FermionField,par().source);
        if (hasGammas_) {
            auto &sol = envGet(ARG(std::map<Gamma::Algebra,FermionField>),getName());
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
            auto &srctmp = envGet(ARG(std::map<Gamma::Algebra,std::vector<FermionField> >), par().source);
            src = &(srctmp.at(Gamma::Algebra::Gamma5));
        }
        if (hasGammas_) {
            auto &sol = envGet(ARG(std::map<Gamma::Algebra,std::vector<FermionField>>),getName());
            solvePropagator(sol,*src);
        } else {
            auto &sol = envGet(std::vector<FermionField>,getName());
            solvePropagator(sol,*src);
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_StagGaugeProp_hpp_
