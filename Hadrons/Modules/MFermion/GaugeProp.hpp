/*
 * GaugeProp.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Guido Cossu <guido.cossu@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Nils Asmussen <n.asmussen@soton.ac.uk>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: pretidav <david.preti@csic.es>
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

#ifndef Hadrons_MFermion_GaugeProp_hpp_
#define Hadrons_MFermion_GaugeProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                                GaugeProp                                   *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class GaugePropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(GaugePropPar,
                                    std::string, source,
                                    std::string, solver);
};

template <typename FImpl>
class TGaugeProp: public Module<GaugePropPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TGaugeProp(const std::string name);
    // destructor
    virtual ~TGaugeProp(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    template <typename TField>
    void setupHelper(void);
    virtual void setup(void);
    // execution
    template <typename TField>
    void executeHelper(void);
    virtual void execute(void);
private:
    void solvePropagator(FermionField &result, FermionField &propPhysical,
                         const FermionField &source);
    template <typename TFImpl = FImpl, IfNotStag<TFImpl> = 0>
    void solvePropagator(PropagatorField &result, PropagatorField &propPhysical,
                         const PropagatorField &source);
    template <typename TFImpl = FImpl, IfStag<TFImpl> = 0>
    void solvePropagator(PropagatorField &result, PropagatorField &propPhysical,
                         const PropagatorField &source);
private:
    unsigned int Ls_;
    Solver       *solver_{nullptr};
};

MODULE_REGISTER_TMP(GaugeProp, TGaugeProp<FIMPL>, MFermion);
MODULE_REGISTER_TMP(ZGaugeProp, TGaugeProp<ZFIMPL>, MFermion);

/******************************************************************************
 *                      TGaugeProp implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TGaugeProp<FImpl>::TGaugeProp(const std::string name)
: Module<GaugePropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TGaugeProp<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().source, par().solver};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TGaugeProp<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_5d"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
template <typename TField>
void TGaugeProp<FImpl>::setupHelper() {
    envTmpLat(FermionField, "tmp");
    if (Ls_ > 1)
    {
        envTmpLat(FermionField, "source", Ls_);
        envTmpLat(FermionField, "sol", Ls_);
    }
    else
    {
        envTmpLat(FermionField, "source");
        envTmpLat(FermionField, "sol");
    }

    if (envHasType(TField,par().source)) {
        envCreateLat(TField, getName());
        if (Ls_ > 1)
        {
            envCreateLat(TField, getName() + "_5d", Ls_);
        }
    } else {
        auto &src = envGet(std::vector<TField>, par().source);

        envCreate(std::vector<TField>, getName(), 1, src.size(),
                  envGetGrid(TField));
        if (Ls_ > 1)
        {
            envCreate(std::vector<TField>, getName() + "_5d", Ls_,
                      src.size(), envGetGrid(TField, Ls_));
        }
    }
}

template <typename FImpl>
void TGaugeProp<FImpl>::setup(void)
{
    Ls_ = env().getObjectLs(par().solver);
    
    if (envHasType(PropagatorField,par().source) || envHasType(std::vector<PropagatorField>,par().source))
        setupHelper<PropagatorField>();
    else if (envHasType(FermionField,par().source) || envHasType(std::vector<FermionField>,par().source)) 
        setupHelper<FermionField>();
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TGaugeProp<FImpl>::solvePropagator(FermionField &prop, 
                                        FermionField &propPhysical,
                                        const FermionField &fullSrc)
{
    auto &solver  = envGet(Solver, par().solver);
    auto &mat     = solver.getFMat();
    
    envGetTmp(FermionField, source);
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, tmp);
    
    LOG(Message) << "Inverting using solver '" << par().solver << "'" 
                 << std::endl;
    // source conversion for 4D sources
    LOG(Message) << "Import source" << std::endl;
    if (!env().isObject5d(par().source))
    {
        if (Ls_ == 1)
        {
           source = fullSrc;
        }
        else
        {
            tmp = fullSrc;
            mat.ImportPhysicalFermionSource(tmp, source);
        }
    }
    // source conversion for 5D sources
    else
    {
        if (Ls_ != env().getObjectLs(par().source))
        {
            HADRONS_ERROR(Size, "Ls mismatch between quark action and source");
        }
        else
        {
            source = fullSrc;
        }
    }
    sol = Zero();
    LOG(Message) << "Solve" << std::endl;
    solver(sol, source);
    LOG(Message) << "Export solution" << std::endl;
    prop = sol;
    // create 4D propagators from 5D one if necessary
    if (Ls_ > 1)
    {
        mat.ExportPhysicalFermionSolution(sol, tmp);
        propPhysical = tmp;
    }
}

template <typename FImpl>
template <typename TFImpl, IfNotStag<TFImpl> >
void TGaugeProp<FImpl>::solvePropagator(PropagatorField &prop, 
                                        PropagatorField &propPhysical,
                                        const PropagatorField &fullSrc)
{
    auto &solver  = envGet(Solver, par().solver);
    auto &mat     = solver.getFMat();
    
    envGetTmp(FermionField, source);
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, tmp);
    
    LOG(Message) << "Inverting using solver '" << par().solver << "'" 
                 << std::endl;
    for (unsigned int s = 0; s < Ns; ++s)
    for (unsigned int c = 0; c < FImpl::Dimension; ++c)
    {
        LOG(Message) << "Inversion for spin= " << s << ", color= " << c
                     << std::endl;
        // source conversion for 4D sources
        LOG(Message) << "Import source" << std::endl;
        if (!env().isObject5d(par().source))
        {
            if (Ls_ == 1)
            {
               PropToFerm<FImpl>(source, fullSrc, s, c);
            }
            else
            {
                PropToFerm<FImpl>(tmp, fullSrc, s, c);
                mat.ImportPhysicalFermionSource(tmp, source);
            }
        }
        // source conversion for 5D sources
        else
        {
            if (Ls_ != env().getObjectLs(par().source))
            {
                HADRONS_ERROR(Size, "Ls mismatch between quark action and source");
            }
            else
            {
                PropToFerm<FImpl>(source, fullSrc, s, c);
            }
        }
        sol = Zero();
        LOG(Message) << "Solve" << std::endl;
        solver(sol, source);
        LOG(Message) << "Export solution" << std::endl;
        FermToProp<FImpl>(prop, sol, s, c);
        // create 4D propagators from 5D one if necessary
        if (Ls_ > 1)
        {
            mat.ExportPhysicalFermionSolution(sol, tmp);
            FermToProp<FImpl>(propPhysical, tmp, s, c);
        }
    }
}

template <typename FImpl>
template <typename TFImpl, IfStag<TFImpl> >
void TGaugeProp<FImpl>::solvePropagator(PropagatorField &prop, 
                                        PropagatorField &propPhysical,
                                        const PropagatorField &fullSrc)
{
    auto &solver  = envGet(Solver, par().solver);
    auto &mat     = solver.getFMat();
    
    envGetTmp(FermionField, source);
    envGetTmp(FermionField, sol);
    envGetTmp(FermionField, tmp);
    
    LOG(Message) << "Inverting using solver '" << par().solver << "'" 
                 << std::endl;
    for (unsigned int c = 0; c < FImpl::Dimension; ++c)
    {
        LOG(Message) << "Inversion for color= " << c
                     << std::endl;
        // source conversion for 4D sources
        LOG(Message) << "Import source" << std::endl;
        if (!env().isObject5d(par().source))
        {
            if (Ls_ == 1)
            {
               PropToFerm<FImpl>(source, fullSrc, c);
            }
            else
            {
                PropToFerm<FImpl>(tmp, fullSrc, c);
                mat.ImportPhysicalFermionSource(tmp, source);
            }
        }
        // source conversion for 5D sources
        else
        {
            if (Ls_ != env().getObjectLs(par().source))
            {
                HADRONS_ERROR(Size, "Ls mismatch between quark action and source");
            }
            else
            {
                PropToFerm<FImpl>(source, fullSrc, c);
            }
        }
        sol = Zero();
        LOG(Message) << "Solve" << std::endl;
        solver(sol, source);
        LOG(Message) << "Export solution" << std::endl;
        FermToProp<FImpl>(prop, sol, c);
        // create 4D propagators from 5D one if necessary
        if (Ls_ > 1)
        {
            mat.ExportPhysicalFermionSolution(sol, tmp);
            FermToProp<FImpl>(propPhysical, tmp, c);
        }
    }
}

template <typename FImpl>
template<typename TField>
void TGaugeProp<FImpl>::executeHelper(void) {

    std::string propName = (Ls_ == 1) ? getName() : (getName() + "_5d");

    if (envHasType(TField, par().source))
    {
        auto &prop         = envGet(TField, propName);
        auto &propPhysical = envGet(TField, getName());
        auto &fullSrc      = envGet(TField, par().source);

        LOG(Message) << "Using source '" << par().source << "'" << std::endl;
        solvePropagator(prop, propPhysical, fullSrc);
    }
    else
    {
        auto &prop         = envGet(std::vector<TField>, propName);
        auto &propPhysical = envGet(std::vector<TField>, getName());
        auto &fullSrc      = envGet(std::vector<TField>, par().source);

        for (unsigned int i = 0; i < fullSrc.size(); ++i)
        {
            LOG(Message) << "Using element " << i << " of source vector '" 
                         << par().source << "'" << std::endl;
            solvePropagator(prop[i], propPhysical[i], fullSrc[i]);
        }
    }

}

template <typename FImpl>
void TGaugeProp<FImpl>::execute(void)
{
    LOG(Message) << "Computing quark propagator '" << getName() << "'"
                 << std::endl;
    
if (envHasType(PropagatorField,par().source) || envHasType(std::vector<PropagatorField>,par().source))
        executeHelper<PropagatorField>();
    else if (envHasType(FermionField,par().source) || envHasType(std::vector<FermionField>,par().source)) 
        executeHelper<FermionField>();
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_GaugeProp_hpp_
