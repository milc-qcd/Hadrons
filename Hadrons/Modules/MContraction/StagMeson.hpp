/*
 * StagMeson.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
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

#ifndef Hadrons_MContraction_StagMeson_hpp_
#define Hadrons_MContraction_StagMeson_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Disconnected loop contractions
 ------------------------------
 
 * options:
 - q_loop: input propagator (string)
 - gammas: gammas: gamma matrices to insert
           (space-separated strings e.g. "GammaT GammaX GammaY") 

           Special values: "all" - perform all possible contractions.
*/

/******************************************************************************
 *                                StagMeson                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;
class StagMesonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagMesonPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, gammas,
                                    std::string, gammaFunc,
                                    std::string, sink,
                                    std::string, sourceShift,
                                    std::string, output);
};

template <typename FImpl>
class TStagMeson: public Module<StagMesonPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    SINK_TYPE_ALIASES(Scalar);

    typedef std::function<LatticeComplex (Gamma::Algebra gamma)> GammaFn;

    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma_snk,
                                        Gamma::Algebra, gamma_src,
                                        std::vector<Complex>, corr);
    };
public:
    // constructor
    TStagMeson(const std::string name);
    // destructor
    virtual ~TStagMeson(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    void parseGammaString();
protected:
    template<typename TField>
    EnableIf<is_lattice<TField>,void> contract(Result &ret, const TField &fSink, const TField &fSrc, Real scale = 1.0, Integer shift = 0);
    template<typename TField>
    EnableIf<is_lattice<TField>,void> contract(Result &ret, const std::vector<TField> &fSink, const std::vector<TField> &fSrc);
    template<typename TField>
    EnableIf<is_lattice<TField>,void> contract(std::vector<Result> &ret, const TField &fSink, const TField &fSrc);
    template<typename TField>
    EnableIf<is_lattice<TField>,void> contract(std::vector<Result> &ret, const std::vector<TField> &fSink, const std::vector<TField> &fSrc);
    template<typename TField>
    EnableIf<is_lattice<TField>,void> contract(std::vector<Result> &ret, const std::map<Gamma::Algebra,std::vector<TField>> &fSink, const std::map<Gamma::Algebra,std::vector<TField>> &fSrc);

    inline void buildProp(PropagatorField &ret, const FermionField &fSink, const FermionField &fSrc, Gamma::Algebra gamma) {

        envGetTmp(FermionField,left);
        envGetTmp(LatticeComplex,stagPhase);

        auto &func = envGet(GammaFn, par().gammaFunc);

        stagPhase = func(gamma);
        left = fSink*stagPhase;

        ret = outerProduct(left,fSrc);
    }
    inline void buildProp(PropagatorField &ret, const PropagatorField &fSink, const PropagatorField &fSrc, Gamma::Algebra gamma) {

        envGetTmp(PropagatorField,left);
        envGetTmp(LatticeComplex,stagPhase);

        auto &func = envGet(GammaFn, par().gammaFunc);

        stagPhase = func(gamma);
        left = fSink*stagPhase;

        ret = left*adj(fSrc);
    }

    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
};

MODULE_REGISTER_TMP(StagMeson, TStagMeson<STAGIMPL>, MContraction);

/******************************************************************************
 *                       TStagMeson implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TStagMeson<FImpl>::TStagMeson(const std::string name)
: Module<StagMesonPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TStagMeson<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q1, par().q2, par().sink, par().gammaFunc};
    if(!par().sourceShift.empty())
        in.push_back(par().sourceShift);
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TStagMeson<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TStagMeson<FImpl>::setup(void)
{
    envTmp(std::vector<GammaPair>,"gammaList",1,0);

    envTmpLat(LatticeComplex,"stagPhase");

    envTmpLat(PropagatorField, "op");
    if ((envHasType(PropagatorField, par().q1) && envHasType(PropagatorField, par().q2))
        || (envHasType(std::vector<PropagatorField>, par().q1) && envHasType(std::vector<PropagatorField>, par().q2))
        || (envHasType(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField>>), par().q1) && envHasType(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField>>), par().q2))) {

        envTmpLat(PropagatorField, "left");

    } else if ((envHasType(FermionField, par().q1) && envHasType(FermionField, par().q2))
        || (envHasType(std::vector<FermionField>, par().q1) && envHasType(std::vector<FermionField>, par().q2))
        || (envHasType(ARG(std::map<Gamma::Algebra,std::vector<FermionField>>), par().q1) && envHasType(ARG(std::map<Gamma::Algebra,std::vector<FermionField>>), par().q2))) {

        envTmpLat(FermionField, "left");

    } else {
        HADRONS_ERROR(Logic,"q1 and q2 must have the same type.");
    }

}

template <typename FImpl>
void TStagMeson<FImpl>::parseGammaString()
{
    envGetTmp(std::vector<GammaPair>,gammaList);
    gammaList.clear();
    // Parse individual contractions from input string.
    gammaList = strToVec<GammaPair>(par().gammas);
}

// execution ///////////////////////////////////////////////////////////////////
template<typename FImpl>
template<typename TField>
EnableIf<is_lattice<TField>,void> TStagMeson<FImpl>::contract(Result &ret, const TField &fSink, const TField &fSrc, Real scale, Integer shift) {

    int offset, nt = env().getDim(Tp);
    std::vector<TComplex>     buf;

    SinkFnScalar &sink = envGet(SinkFnScalar, par().sink);

    envGetTmp(PropagatorField, op);

    buildProp(op, fSink,fSrc,ret.gamma_snk);

    buf = sink(trace(op));
    for (unsigned int t = 0; t < nt; ++t)
    {
        offset = mod(t+shift,nt);
        if (scale > 0.0)
            ret.corr[t] += (TensorRemove(buf[offset])/scale);
        else 
            ret.corr[t] += TensorRemove(buf[offset]);
    }
}

template<typename FImpl>
template<typename TField>
EnableIf<is_lattice<TField>,void> TStagMeson<FImpl>::contract(Result &ret, const std::vector<TField> &fSink, const std::vector<TField> &fSrc) {

    std::vector<Integer> shifts;
    if (!par().sourceShift.empty()) {
        shifts = envGet(std::vector<Integer>, par().sourceShift);
    }

    for (int i = 0; i < fSink.size(); i++) {
        LOG(Message) << "Contracting element i = " << i << std::endl;
        if (!par().sourceShift.empty()) {
            LOG(Message) << "Shifting correlator to (t0 = " << shifts[i] << ")" << std::endl;
            contract(ret,fSink[i],fSrc[i],fSink.size(), shifts[i]);
        } else {
            contract(ret,fSink[i],fSrc[i],fSink.size());
        }
    }
}
template<typename FImpl>
template<typename TField>
EnableIf<is_lattice<TField>,void> TStagMeson<FImpl>::contract(std::vector<Result> &ret, const TField &fSink, const TField &fSrc) {

    for (unsigned int i = 0; i < ret.size(); i++) {
        LOG(Message) << "Contracting gammas: " << ret[i].gamma_snk << " (sink), " << ret[i].gamma_src << " (source) " << std::endl;
        contract(ret[i],fSink,fSrc);
    }
}

template<typename FImpl>
template<typename TField>
EnableIf<is_lattice<TField>,void> TStagMeson<FImpl>::contract(std::vector<Result> &ret, const std::vector<TField> &fSink, const std::vector<TField> &fSrc) {

    for (unsigned int i = 0; i < ret.size(); i++) {
        LOG(Message) << "Contracting gammas: " << ret[i].gamma_snk << " (sink), " << ret[i].gamma_src << " (source) " << std::endl;
        contract(ret[i],fSink,fSrc);
    }
}

template<typename FImpl>
template<typename TField>
EnableIf<is_lattice<TField>,void> TStagMeson<FImpl>::contract(std::vector<Result> &ret, const std::map<Gamma::Algebra,std::vector<TField>> &fSink, 
                                                                const std::map<Gamma::Algebra,std::vector<TField>> &fSrc) {

    for (unsigned int i = 0; i < ret.size(); i++) {
        LOG(Message) << "Contracting gammas: " << ret[i].gamma_snk << " (sink), " << ret[i].gamma_src << " (source) " << std::endl;
        // Always use Gamma5 solve for the sink quark
        contract(ret[i],fSink.at(Gamma::Algebra::Gamma5),fSrc.at(ret[i].gamma_src));
    }
}

template <typename FImpl>
void TStagMeson<FImpl>::execute(void)
{

    LOG(Message) << "Computing meson contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "' and '" << par().q2 << "'"
                 << std::endl;

    std::vector<Result> result;
    unsigned int nt = env().getDim(Tp);

    parseGammaString();

    envGetTmp(std::vector<GammaPair>,gammaList);
    result.resize(gammaList.size());

    for (unsigned int i = 0; i < result.size(); ++i)
    {
        result[i].gamma_snk = gammaList[i].first;
        result[i].gamma_src = gammaList[i].second;
        result[i].corr.resize(nt, 0.0);
    }

    if (envHasType(PropagatorField, par().q1)) {
        auto &q1  = envGet(PropagatorField, par().q1);
        auto &q2  = envGet(PropagatorField, par().q2);

        contract(result,q1,q2);

    } else if (envHasType(std::vector<PropagatorField>, par().q1)) {
        auto &q1  = envGet(std::vector<PropagatorField>, par().q1);
        auto &q2  = envGet(std::vector<PropagatorField>, par().q2);

        contract(result,q1,q2);

    } else if (envHasType(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField>>), par().q1)) {
        auto &q1  = envGet(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField>>), par().q1);
        auto &q2  = envGet(ARG(std::map<Gamma::Algebra,std::vector<PropagatorField>>), par().q2);

        contract(result,q1,q2);

    } else if (envHasType(FermionField, par().q1)) {
        auto &q1  = envGet(FermionField, par().q1);
        auto &q2  = envGet(FermionField, par().q2);

        contract(result,q1,q2);

    } else if (envHasType(std::vector<FermionField>, par().q1)) {
        auto &q1  = envGet(std::vector<FermionField>, par().q1);
        auto &q2  = envGet(std::vector<FermionField>, par().q2);

        contract(result,q1,q2);
    } else {
        auto &q1  = envGet(ARG(std::map<Gamma::Algebra,std::vector<FermionField>>), par().q1);
        auto &q2  = envGet(ARG(std::map<Gamma::Algebra,std::vector<FermionField>>), par().q2);

        contract(result,q1,q2);
    }

    saveResult(par().output, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_StagMeson_hpp_
