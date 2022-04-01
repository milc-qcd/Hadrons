/*
 * MesonMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MContraction_MesonMILC_hpp_
#define Hadrons_MContraction_MesonMILC_hpp_

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
 *                                MesonMILC                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;
typedef std::map<Gamma::Algebra, LatticeComplex> PhaseMap;

class MesonMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonMILCPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, gammas,
                                    std::string, sink,
                                    std::string, sourceShift,
                                    std::string, output);
};

template <typename FImpl>
class TMesonMILC: public Module<MesonMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    SINK_TYPE_ALIASES(Scalar);
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
    TMesonMILC(const std::string name);
    // destructor
    virtual ~TMesonMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    template<typename T = FImpl>
    IfNotStag<T,void> parseGammaString();
    template<typename T = FImpl>
    IfStag<T,void> parseGammaString();
protected:
    template<typename TField>
    EnableIf<is_lattice<TField>,void> contract(std::vector<Result> &ret, const TField &fSink, const TField &fSrc, Real scale = 1.0, Integer shift = 0);
    template<typename TField>
    EnableIf<is_lattice<TField>,void> contract(std::vector<Result> &ret, const std::vector<TField> &fSink, const std::vector<TField> &fSrc);

    template<typename TField, typename T = FImpl>
    inline IfStag<T,void> applyGamma(TField &ret, const TField &f, Gamma::Algebra gamma) {
        envGetTmp(PhaseMap,stag_phase);
        ret = stag_phase.at(gamma)*f;
    }
    template<typename TField, typename T = FImpl>
    inline IfNotStag<T,void> applyGamma(TField &ret, const TField &f,Gamma::Algebra gamma) {
        Gamma g(gamma);
        ret = g*f;
    }

    inline void buildProp(PropagatorField &ret, const FermionField &fSink, const FermionField &fSrc, int gamma) {
        envGetTmp(std::vector<GammaPair>,gammaList);
        envGetTmp(FermionField,left);
        envGetTmp(FermionField,right);
        applyGamma(left, fSink, ((Gamma::Algebra)(gammaList[gamma].first)));
        applyGamma(right, fSrc, gammaList[gamma].second);

        ret = outerProduct(left,right);
    }
    inline void buildProp(PropagatorField &ret, const PropagatorField &fSink, const PropagatorField &fSrc, int gamma) {
        envGetTmp(std::vector<GammaPair>,gammaList);
        envGetTmp(PropagatorField,left);
        envGetTmp(PropagatorField,right);
        applyGamma(left, fSink, gammaList[gamma].first);
        applyGamma(right, fSrc, gammaList[gamma].second);

        ret = left*adj(right);
    }

    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
};

MODULE_REGISTER_TMP(MesonMILC, TMesonMILC<STAGIMPL>, MContraction);

/******************************************************************************
 *                       TMesonMILC implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TMesonMILC<FImpl>::TMesonMILC(const std::string name)
: Module<MesonMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TMesonMILC<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q1, par().q2, par().sink};
    if(!par().sourceShift.empty())
        in.push_back(par().sourceShift);
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TMesonMILC<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TMesonMILC<FImpl>::setup(void)
{
    envTmp(std::vector<GammaPair>,"gammaList",1,0);

    PhaseMap tmp = {
        {Gamma::Algebra::GammaX,envGetGrid(LatticeComplex)},
        {Gamma::Algebra::GammaY,envGetGrid(LatticeComplex)},
        {Gamma::Algebra::GammaZ,envGetGrid(LatticeComplex)},
        {Gamma::Algebra::Gamma5,envGetGrid(LatticeComplex)}};

    envTmp(PhaseMap,"stag_phase",1,tmp);

    envTmpLat(PropagatorField, "op");
    if ((envHasType(PropagatorField, par().q1) && envHasType(PropagatorField, par().q2))
        || (envHasType(std::vector<PropagatorField>, par().q1) && envHasType(std::vector<PropagatorField>, par().q2))) {

        envTmpLat(PropagatorField, "left");
        envTmpLat(PropagatorField, "right");

    } else if ((envHasType(FermionField, par().q1) && envHasType(FermionField, par().q2))
        || (envHasType(std::vector<FermionField>, par().q1) && envHasType(std::vector<FermionField>, par().q2))) {

        envTmpLat(FermionField, "left");
        envTmpLat(FermionField, "right");

    } else {
        HADRONS_ERROR(Logic,"q1 and q2 must have the same type.");
    }

}

template <typename FImpl>
template <typename T>
IfStag<T,void> TMesonMILC<FImpl>::parseGammaString()
{
    std::vector<Gamma::Algebra> keys = {
        Gamma::Algebra::GammaX,
        Gamma::Algebra::GammaY,
        Gamma::Algebra::GammaZ,
        Gamma::Algebra::Gamma5
    };


    envGetTmp(PhaseMap,stag_phase);
    for (const auto &key:keys) {
        stag_phase.insert({key,envGetGrid(LatticeComplex)});
        stag_phase.at(key) = 1.0;
    }

    Lattice<iScalar<vInteger> > x(env().getGrid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(env().getGrid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(env().getGrid()); LatticeCoordinate(z,2);
    
    stag_phase.at(Gamma::Algebra::GammaX) = where( mod(x,2)==(Integer)0, stag_phase.at(Gamma::Algebra::GammaX), -stag_phase.at(Gamma::Algebra::GammaX));
    stag_phase.at(Gamma::Algebra::GammaY) = where( mod(y,2)==(Integer)0, stag_phase.at(Gamma::Algebra::GammaY), -stag_phase.at(Gamma::Algebra::GammaY));
    stag_phase.at(Gamma::Algebra::GammaZ) = where( mod(z,2)==(Integer)0, stag_phase.at(Gamma::Algebra::GammaZ), -stag_phase.at(Gamma::Algebra::GammaZ));

    envGetTmp(std::vector<GammaPair>,gammaList);
    gammaList.clear();
    // Determine gamma matrices to insert at source/sink.
    if (par().gammas.compare("all") == 0)
    {
        // Do all contractions.
        for (const auto& key1:keys)
        {
            for (const auto& key2:keys)
                gammaList.push_back(std::make_pair(key1, key2));
        }
    }
    else
    {
        // Parse individual contractions from input string.
        gammaList = strToVec<GammaPair>(par().gammas);
        for (const auto& pair:gammaList) {
            switch(pair.first) {
                case Gamma::Algebra::GammaX:
                case Gamma::Algebra::GammaY:
                case Gamma::Algebra::GammaZ:
                case Gamma::Algebra::Gamma5:
                    break;
                default:
                    HADRONS_ERROR(Implementation,"your gamma is not supported for stag fields");
            }
            switch(pair.second) {
                case Gamma::Algebra::GammaX:
                case Gamma::Algebra::GammaY:
                case Gamma::Algebra::GammaZ:
                case Gamma::Algebra::Gamma5:
                    break;
                default:
                    HADRONS_ERROR(Implementation,"your gamma is not supported for stag fields");
            }
        }
    } 
}

template <typename FImpl>
template <typename T>
IfNotStag<T,void> TMesonMILC<FImpl>::parseGammaString()
{
    envGetTmp(std::vector<GammaPair>,gammaList);
    gammaList.clear();
    // Determine gamma matrices to insert at source/sink.
    if (par().gammas.compare("all") == 0)
    {
        // Do all contractions.
        for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
        {
            for (unsigned int j = 1; j < Gamma::nGamma; j += 2)
            {
                gammaList.push_back(std::make_pair((Gamma::Algebra)i, 
                                                   (Gamma::Algebra)j));
            }
        }
    }
    else
    {
        // Parse individual contractions from input string.
        gammaList = strToVec<GammaPair>(par().gammas);
    } 
}

// execution ///////////////////////////////////////////////////////////////////
template<typename FImpl>
template<typename TField>
EnableIf<is_lattice<TField>,void> TMesonMILC<FImpl>::contract(std::vector<Result> &ret, const TField &fSink, const TField &fSrc, Real scale, Integer shift) {

    int offset, nt = env().getDim(Tp);
    std::vector<TComplex>     buf;

    LOG(Message) << "(using sink '" << par().sink << "')" << std::endl;
    SinkFnScalar &sink = envGet(SinkFnScalar, par().sink);

    envGetTmp(PropagatorField, op);
    envGetTmp(std::vector<GammaPair>,gammaList);

    for (unsigned int g = 0; g < gammaList.size(); ++g)
    {
        LOG(Message) << "Using gammas: " << gammaList[g].first << " (sink), " << gammaList[g].second << " (source) " << std::endl;

        buildProp(op, fSink,fSrc,g);

        buf = sink(trace(op));
        for (unsigned int t = 0; t < nt; ++t)
        {
            offset = mod(t+shift,nt);
            if (scale > 0.0)
                ret[g].corr[t] += (TensorRemove(buf[offset])/scale);
            else 
                ret[g].corr[t] += TensorRemove(buf[offset]);
        }
    }
}
template<typename FImpl>
template<typename TField>
EnableIf<is_lattice<TField>,void> TMesonMILC<FImpl>::contract(std::vector<Result> &ret, const std::vector<TField> &fSink, const std::vector<TField> &fSrc) {

    std::vector<Integer> shifts;
    if (!par().sourceShift.empty()) {
        LOG(Message) << "Using " << par().sourceShift << " to shift time axis." << std::endl;
        shifts = envGet(std::vector<Integer>, par().sourceShift);
    }

    for (int i = 0; i < fSink.size(); i++) {
        if (!par().sourceShift.empty()) {
            LOG(Message) << "Shifting correlator by " << shifts[i] << std::endl;
            contract(ret,fSink[i],fSrc[i],fSink.size(), shifts[i]);
        } else {
            contract(ret,fSink[i],fSrc[i],fSink.size());
        }
    }


}

template <typename FImpl>
void TMesonMILC<FImpl>::execute(void)
{

    LOG(Message) << "Computing meson contractions '" << getName() << "' using"
                 << " quarks '" << par().q1 << "' and '" << par().q2 << "'"
                 << std::endl;

    std::vector<Result> result;
    unsigned int nt = env().getDim(Tp);

    parseGammaString();

    envGetTmp(std::vector<GammaPair>,gammaList);
    result.resize(gammaList.size());

    for (unsigned int g = 0; g < result.size(); ++g)
    {
        result[g].gamma_snk = gammaList[g].first;
        result[g].gamma_src = gammaList[g].second;
        result[g].corr.resize(nt, 0.0);
    }

    if (envHasType(PropagatorField, par().q1)) {
        auto &q1  = envGet(PropagatorField, par().q1);
        auto &q2  = envGet(PropagatorField, par().q2);

        contract(result,q1,q2);

    } else if (envHasType(std::vector<PropagatorField>, par().q1)) {
        auto &q1  = envGet(std::vector<PropagatorField>, par().q1);
        auto &q2  = envGet(std::vector<PropagatorField>, par().q2);

        contract(result,q1,q2);

    } else if (envHasType(FermionField, par().q1)) {
        auto &q1  = envGet(FermionField, par().q1);
        auto &q2  = envGet(FermionField, par().q2);

        contract(result,q1,q2);

    } else {
        auto &q1  = envGet(std::vector<FermionField>, par().q1);
        auto &q2  = envGet(std::vector<FermionField>, par().q2);

        contract(result,q1,q2);
    }

    saveResult(par().output, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_MesonMILC_hpp_
