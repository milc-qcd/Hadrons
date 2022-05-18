/*
 * LMAProp.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MFermion_LMAProp_hpp_
#define Hadrons_MFermion_LMAProp_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                       Create all-to-all V & W vectors                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MFermion)

class LMAPropPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(LMAPropPar,
                                  std::string, source,
                                  std::string, action,
                                  std::string, gammas,
                                  std::string, gammaFunc,
                                  std::string, lowModes);
};

template <typename FImpl>
class TLMAProp : public Module<LMAPropPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef std::function<LatticeComplex (Gamma::Algebra gamma)> GammaFn;
private:
    HADRONS_DEFINE_setProp_setFerm(FImpl);
public:
    // constructor
    TLMAProp(const std::string name);
    // destructor
    virtual ~TLMAProp(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);

    // setup
    virtual void setup(void);

    // execute
    virtual void execute(void);
private:
    bool hasGammas_;
};

MODULE_REGISTER_TMP(StagLMAProp, TLMAProp<STAGIMPL>, MFermion);

/******************************************************************************
 *                       TLMAProp implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLMAProp<FImpl>::TLMAProp(const std::string name)
: Module<LMAPropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLMAProp<FImpl>::getInput(void)
{
    std::vector<std::string> in {par().action, par().source};

    if (!par().lowModes.empty()) {
        in.push_back(par().lowModes);
        in.push_back(par().lowModes+"_evalM");
    }
    
    hasGammas_ = !par().gammas.empty();

    if (hasGammas_) {
        in.push_back(par().gammaFunc);
    }

    return in;
}

template <typename FImpl>
std::vector<std::string> TLMAProp<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

/******************************************************************************
 *              TLMAProp setup                                         *
 ******************************************************************************/
template <typename FImpl>
void TLMAProp<FImpl>::setup(void)
{
    auto        &action     = envGet(FMat, par().action);
    int         Ls          = env().getObjectLs(par().action);

    if (Ls > 1) {
       HADRONS_ERROR(Argument, "Ls > 1 not implemented");
    }

    auto &source = envGet(std::vector<PropagatorField>, par().source);

    envTmpLat(FermionField, "evecM");
    envTmpLat(FermionField, "evecMDag");
    envTmpLat(PropagatorField, "prop");
    envTmpLat(PropagatorField, "propDag");
    envTmpLat(FermionField, "ferm");
    envTmpLat(LatticeComplex,"stagPhase");
    envTmpLat(LatticeComplex,"flipOddcb");
    envTmp(FermionField, "tempRb", 1, envGetRbGrid(FermionField));
    envTmp(LatticeComplex, "tempRbScalar", 1, envGetRbGrid(LatticeComplex));
    envTmp(FermionField, "evecNeg", 1, envGetRbGrid(FermionField));

    envGetTmp(FermionField,evecM);
    envGetTmp(FermionField,evecMDag);
    envGetTmp(FermionField,tempRb);
    envGetTmp(LatticeComplex,tempRbScalar);
    envGetTmp(FermionField,evecNeg);
    envGetTmp(FermionField, ferm);
    envGetTmp(PropagatorField, prop);
    envGetTmp(PropagatorField, propDag);

    evecM    = Zero();
    evecMDag = Zero();
    tempRb   = Zero();
    ferm     = Zero();
    prop     = Zero();
    propDag  = Zero();
    evecNeg  = Zero();
    tempRbScalar   = Zero();

    envTmp(std::vector<Gamma::Algebra>,"gammaList",1,0);
    envGetTmp(std::vector<Gamma::Algebra>,gammaList);
    gammaList.clear();

    if (hasGammas_)  {
        gammaList = strToVec<Gamma::Algebra>(par().gammas);

        std::map<Gamma::Algebra,std::vector<PropagatorField>> dummy;
        envCreate(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField>>), getName(), 1, dummy);
        auto &sol = envGet(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField>>), getName());
        for (auto & gamma:gammaList) {
            sol.insert({gamma,std::vector<PropagatorField>(source.size(),envGetGrid(PropagatorField))});
            for (auto & s:sol.at(gamma)) {
                s = Zero();
            }
        }
    } else {
        envCreate(std::vector<PropagatorField>, getName(), 1, 
                  source.size(), envGetGrid(PropagatorField));
        auto &sol = envGet(std::vector<PropagatorField>,getName());
        for (auto & s:sol) {
            s = Zero();
        }
    }
}

/******************************************************************************
 *              TLMAProp execution                                     *
 ******************************************************************************/
template <typename FImpl>
void TLMAProp<FImpl>::execute(void)
{
    envGetTmp(LatticeComplex,flipOddcb);
    envGetTmp(LatticeComplex,tempRbScalar);
    envGetTmp(FermionField,evecM);
    envGetTmp(FermionField,evecMDag);
    envGetTmp(FermionField,tempRb);
    envGetTmp(FermionField,evecNeg);
    envGetTmp(FermionField,ferm);
    envGetTmp(PropagatorField,prop);
    envGetTmp(PropagatorField,propDag);

    auto &action = envGet(FMat, par().action);

    typename std::vector<Real>::iterator it_eval;

    auto &source  = envGet(std::vector<PropagatorField>, par().source);
    auto &evals   = envGet(std::vector<ComplexD>, par().lowModes+"_evalM");
    auto &evecs   = envGet(std::vector<FermionField>, par().lowModes);

    int cb = evecs[0].Checkerboard();
    int cbNeg = (cb==Even) ? Odd : Even;
    

    // Normalize vectors so that checkerboard has magnitude 1/sqrt(2)
    RealD norm = 1/::sqrt(2*norm2(evecs[0]));

    // Build scalar field that is negative on Odd cb. Used for MDag evecs
    flipOddcb = 1.;
    tempRbScalar.Checkerboard() = Odd;
    tempRbScalar = -1.;
    setCheckerboard(flipOddcb,tempRbScalar);

    evecNeg.Checkerboard() = cbNeg;
    tempRb.Checkerboard() = cbNeg;

    for (int i=0;i<evecs.size();i++) {

        ComplexD eval_D = ComplexD(0.0,evals[i].imag());

        // Build complementary checkerboard for evec
        action.Meooe(evecs[i], tempRb);
        evecNeg = (1.0/eval_D) * tempRb;

        setCheckerboard(evecM,evecNeg);
        setCheckerboard(evecM,evecs[i]);
        
        // Normalize eigenvectors to 1
        evecM *= norm;

        // Build MDag evec by flipping sign of Odd checkerboard
        evecMDag = flipOddcb*evecM;

        auto evalInv    = 1.0/evals[i];
        auto evalDagInv = conjugate(evalInv);

        if (hasGammas_) {
            envGetTmp(std::vector<Gamma::Algebra>,gammaList);
            envGetTmp(LatticeComplex,stagPhase);

            auto &sol   = envGet(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField>>), getName());
            auto &func = envGet(GammaFn, par().gammaFunc);

            for (auto &gamma:gammaList) {
                stagPhase = func(gamma);
                for (int j=0;j<source.size();j++) {
                    const PropagatorField &src = source[j];

                    for (int k=0;k<FImpl::Dimension;k++) {
                        setFerm(ferm,src,k);
                        ferm *= stagPhase;
                        auto ip = innerProduct(evecM,ferm)*evalInv;
                        auto ipDag = innerProduct(evecMDag,ferm)*evalDagInv;
                        ferm = ip*evecM;
                        setProp(prop,ferm,k);
                        ferm = ipDag*evecMDag;
                        setProp(propDag,ferm,k);
                    }
                    sol.at(gamma)[j] += prop;
                    sol.at(gamma)[j] += propDag;
                }
            }
        } else {
            auto &sol   = envGet(std::vector<PropagatorField>, getName());
            for (int j=0;j<source.size();j++) {
                const PropagatorField &src = source[j];

                for (int k=0;k<FImpl::Dimension;k++) {
                    setFerm(ferm,src,k);
                    auto ip = innerProduct(evecM,ferm)*evalInv;
                    auto ipDag = innerProduct(evecMDag,ferm)*evalDagInv;
                    ferm = ip*evecM;
                    setProp(prop,ferm,k);
                    ferm = ipDag*evecMDag;
                    setProp(propDag,ferm,k);
                }
                sol[j] += prop;
                sol[j] += propDag;
            }
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_LMAProp_hpp_
