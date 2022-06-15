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
 *                       Calculate Low Mode Average Prop                      *
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

    inline void projectHelper(FermionField& sol, const FermionField& src);
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

    in.push_back(par().lowModes);
    in.push_back(par().lowModes+"_evalM");
    in.push_back(par().gammaFunc);
    if (par().gammas.empty()) {
        HADRONS_ERROR(Logic,"Must provide a list of gammas to " + getName());
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

    envTmpLat(FermionField, "fermSrc");
    envTmpLat(FermionField, "fermSol");
    envTmpLat(PropagatorField, "prop");
    envTmpLat(LatticeComplex,"stagPhase");
    envTmp(FermionField, "rbFerm", 1, envGetRbGrid(FermionField));
    envTmp(FermionField, "rbFermNeg", 1, envGetRbGrid(FermionField));
    envTmp(FermionField, "MrbFermNeg", 1, envGetRbGrid(FermionField));
    envTmp(FermionField, "rbTemp", 1, envGetRbGrid(FermionField));
    envTmp(FermionField, "rbTempNeg", 1, envGetRbGrid(FermionField));

    envGetTmp(FermionField, fermSrc);
    envGetTmp(FermionField, fermSol);
    envGetTmp(PropagatorField, prop);
    envGetTmp(FermionField, rbTemp);
    envGetTmp(FermionField, rbTempNeg);
    envGetTmp(FermionField, rbFerm);
    envGetTmp(FermionField, rbFermNeg);
    envGetTmp(FermionField, MrbFermNeg);

    fermSrc    = Zero();
    fermSol    = Zero();
    prop       = Zero();
    rbFerm     = Zero();
    rbFermNeg  = Zero();
    MrbFermNeg = Zero();

    envTmp(std::vector<Gamma::Algebra>,"gammaList",1,0);
    envGetTmp(std::vector<Gamma::Algebra>,gammaList);
    gammaList.clear();

    gammaList = strToVec<Gamma::Algebra>(par().gammas);

    if ( envHasType(std::vector<PropagatorField>,par().source)) {
        auto &source = envGet(std::vector<PropagatorField>, par().source);

        std::map<Gamma::Algebra,std::vector<PropagatorField>> dummy;
        envCreate(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField>>), getName(), 1, dummy);

        auto &sol = envGet(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField>>), getName());
        for (auto & gamma:gammaList) {
            sol.insert({gamma,std::vector<PropagatorField>(source.size(),envGetGrid(PropagatorField))});
            for (auto & s:sol.at(gamma)) {
                s = Zero();
            }
        }
    } else if ( envHasType(std::vector<FermionField>,par().source)) {
        auto &source = envGet(std::vector<FermionField>, par().source);

        std::map<Gamma::Algebra,std::vector<FermionField>> dummy;
        envCreate(ARG(std::map<Gamma::Algebra,std::vector<FermionField>>), getName(), 1, dummy);

        auto &sol = envGet(ARG(std::map<Gamma::Algebra,std::vector<FermionField>>), getName());
        for (auto & gamma:gammaList) {
            sol.insert({gamma,std::vector<FermionField>(source.size(),envGetGrid(FermionField))});
            for (auto & s:sol.at(gamma)) {
                s = Zero();
            }
        }
    } else {
       HADRONS_ERROR(Argument, "source parameter: '" + par().source + "' must be a vector of FermionFields or PropagatorFields.");
    }
}

template <typename FImpl>
inline void TLMAProp<FImpl>::projectHelper(FermionField& sol, const FermionField& src) {

    envGetTmp(FermionField,rbTemp);
    envGetTmp(FermionField,rbTempNeg);
    envGetTmp(FermionField,rbFerm);
    envGetTmp(FermionField,rbFermNeg);
    envGetTmp(FermionField,MrbFermNeg);

    auto &action = envGet(FMat, par().action);

    auto &evals   = envGet(std::vector<ComplexD>, par().lowModes+"_evalM");
    auto &evecs   = envGet(std::vector<FermionField>, par().lowModes);

    int cb = evecs[0].Checkerboard();
    int cbNeg = (cb==Even) ? Odd : Even;

    // Normalize vectors so that checkerboard has magnitude 1/sqrt(2)
    // Extra factor of 2 accounts for contributions from M and Mdag evecs
    RealD norm = 1./::sqrt(norm2(evecs[0]));

    rbTemp = Zero();
    rbTemp.Checkerboard() = cb;
    rbTempNeg = Zero();
    rbTempNeg.Checkerboard() = cb;

    rbFerm.Checkerboard() = cb;
    rbFermNeg.Checkerboard() = cbNeg;
    MrbFermNeg.Checkerboard() = cb;

    pickCheckerboard(cb,rbFerm,src);
    pickCheckerboard(cbNeg,rbFermNeg,src);

    action.MeooeDag(rbFermNeg, MrbFermNeg); // Move cbNeg component of source to cb

    // Add up source vector projection onto provided evec checkerboard
    // [ lam*(|e> + |o>)(<e| + <o|)  +  conj(lam)*(|e> - |o>)(<e| - <o|) ] |psi>
    for (int k=evecs.size()-1;k >= 0;k--) {
        const FermionField& e = evecs[k];

        const RealD mass     = evals[k].real();
        const RealD lam_D    = evals[k].imag();
        const RealD invlam_D = 1./lam_D; 
        const RealD invmag   = 1./(pow(mass,2)+pow(lam_D,2));
        const ComplexD ip    = TensorRemove(innerProduct(e,rbFerm))*invmag;
        const ComplexD ipNeg = TensorRemove(innerProduct(e,MrbFermNeg))*invmag;

        axpy(rbTemp,    mass*ip+ipNeg,   e,rbTemp);
        axpy(rbTempNeg, mass*ipNeg*invlam_D*invlam_D-ip, e,rbTempNeg);
    }

    action.Meooe(rbTempNeg, rbFermNeg); // Move projection back to cbNeg checkerboard

    setCheckerboard(sol,rbTemp);
    setCheckerboard(sol,rbFermNeg);

    sol *= norm;
}

/******************************************************************************
 *              TLMAProp execution                                     *
 ******************************************************************************/
template <typename FImpl>
void TLMAProp<FImpl>::execute(void)
{
    envGetTmp(FermionField,fermSrc);
    envGetTmp(FermionField,fermSol);
    envGetTmp(PropagatorField,prop);
    envGetTmp(LatticeComplex,stagPhase);
    envGetTmp(std::vector<Gamma::Algebra>,gammaList);

    auto &func  = envGet(GammaFn, par().gammaFunc);

    if ( envHasType(std::vector<PropagatorField>,par().source)) {
        auto &sol   = envGet(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField>>), getName());
        auto &source  = envGet(std::vector<PropagatorField>, par().source);

        for (auto &gamma:gammaList) {
            stagPhase = func(gamma);

            for (int i=0;i<source.size();i++) {
                const PropagatorField& src = source[i];

                prop = stagPhase*src;

                for (int j=0;j<FImpl::Dimension;j++) {

                    setFerm(fermSrc,prop,j);
                    projectHelper(fermSol,fermSrc);
                    setProp(sol.at(gamma)[i],fermSol,j);
                }
            }
        }
    } else {
        auto &sol   = envGet(ARG(std::map<Gamma::Algebra,std::vector<FermionField>>), getName());
        auto &source  = envGet(std::vector<FermionField>, par().source);

        for (auto &gamma:gammaList) {
            stagPhase = func(gamma);

            for (int i=0;i<source.size();i++) {
                const FermionField& src = source[i];

                fermSrc = stagPhase*src;
                projectHelper(sol.at(gamma)[i],fermSrc);
            }
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MFermion_LMAProp_hpp_
