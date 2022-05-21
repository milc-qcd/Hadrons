/*
 * ModifyEigenPackMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MUtilities_ModifyEigenPackMILC_hpp_
#define Hadrons_MUtilities_ModifyEigenPackMILC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Load eigen vectors/values package                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class ModifyEigenPackMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ModifyEigenPackMILCPar,
                                    std::string,  eigenPack,
                                    std::string,  checkerSwapAction,
                                    bool,         evenEigen,
                                    bool,         normalizeCheckerboard,
                                    double,       mass);
    ModifyEigenPackMILCPar() {
        mass = 0.0;
        normalizeCheckerboard = true;
        evenEigen = false;
    };
};

template <typename FImpl, typename Pack>
class TModifyEigenPackMILC: public Module<ModifyEigenPackMILCPar>
{
public:
    typedef typename Pack::Field   Field;
    typedef BaseEigenPack<Field>   BasePack;

    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TModifyEigenPackMILC(const std::string name);
    // destructor
    virtual ~TModifyEigenPackMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ModifyEigenPackMILC, ARG(TModifyEigenPackMILC<STAGIMPL,BaseFermionEigenPack<STAGIMPL> >), MUtilities);

/******************************************************************************
 *                    TModifyEigenPackMILC implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TModifyEigenPackMILC<FImpl,Pack>::TModifyEigenPackMILC(const std::string name)
: Module<ModifyEigenPackMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TModifyEigenPackMILC<FImpl,Pack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};

    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TModifyEigenPackMILC<FImpl,Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_eval", getName() + "_evalM"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TModifyEigenPackMILC<FImpl,Pack>::setup(void)
{
    auto Ls = env().getObjectLs(par().eigenPack);

    auto &epack = envGet(BasePack, par().eigenPack);

    envCreate(std::vector<Field>,getName(), Ls, 0, envGetRbGrid(Field, Ls));
    envCreate(std::vector<RealD>,getName() + "_eval", Ls, 0);
    envCreate(std::vector<ComplexD>,getName() + "_evalM", Ls, 0);
    envTmp(FermionField, "tempRb", 1, envGetRbGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TModifyEigenPackMILC<FImpl,Pack>::execute(void)
{
    int  Ls = env().getObjectLs(par().eigenPack);
    auto &epack = envGet(BasePack, par().eigenPack);

    auto &evec = envGet(std::vector<Field>,getName());
    auto &eval = envGet(std::vector<RealD>,getName() + "_eval");
    auto &evalM = envGet(std::vector<ComplexD>,getName() + "_evalM");

    bool normalizeCheckerboard = par().normalizeCheckerboard;
    bool evenEigen             = par().evenEigen;
    double mass                = par().mass;

    eval.insert(eval.end(),epack.eval.begin(),epack.eval.end());
    evalM.resize(eval.size(),0.0);

    LOG(Message) << "Eigenvalues of the Dirac operator, i.e. mass + i*lambda_D. are stored in '" << getName()+"_evalM" << std::endl;

    Real m = 2*mass;

    for (int i=0;i<eval.size();i++) {
        evalM[i] = ComplexD(m,sqrt(eval[i]));
    }

    if (mass > 0.0) {
        m = ::pow(m,2);
        for (auto &lam:eval) {
            lam += m;
        }        
        LOG(Message) << "Shifted eigenvalues by mass (including MILC factor of 2) = " << m << std::endl;
    }

    int cb = (evenEigen ? Even : Odd);
    if (!par().checkerSwapAction.empty()) {
        LOG(Message) << "Swapping checkerboard from " << (evenEigen?"Odd to ":"Even to ") << (evenEigen?"Even":"Odd")<< std::endl;

        evec.resize(epack.evec.size(),envGetRbGrid(Field, Ls));

        auto &action = envGet(FMat, par().checkerSwapAction);
        int cbNeg = (!evenEigen ? Even : Odd);

        envGetTmp(FermionField,tempRb);
        tempRb = Zero();
        tempRb.Checkerboard() = cb;

        for (int i = 0; i < epack.evec.size();i++) {

            epack.evec[i].Checkerboard() = cbNeg;
            action.Meooe(epack.evec[i],tempRb);

            if (mass == 0.0) {
                evec[i] = (1.0/evalM[i])*tempRb;
            } else {
                evec[i] = (1.0/(evalM[i]-evalM[i].real()))*tempRb;
            }
        }
    } else {
        evec.insert(evec.end(),epack.evec.begin(),epack.evec.end());
    }

    ComplexD norm(1.0,0.0);
    if (normalizeCheckerboard) {
        norm *= 1.0/sqrt(2.0);
        LOG(Message) << "Normalizing eigenvectors by 1/sqrt(2)" << std::endl;
    }
    LOG(Message) << "Setting eigenvector checkerboard to " << (evenEigen ? "'Even'" : "'Odd'" ) << std::endl;
    for (auto &e:evec) {
        e *= norm;
        e.Checkerboard() = cb;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_ModifyEigenPackMILC_hpp_
