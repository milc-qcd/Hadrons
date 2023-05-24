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

/******************************************************************************
 *                                StagMeson                                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class StagMesonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagMesonPar,
                                    std::string, source,
                                    std::string, sink,
                                    std::string, sourceGammas,
                                    SpinTasteParams, sinkSpinTaste,
                                    std::string, sinkFunc,
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

    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        StagGamma::StagAlgebra, gamma_sink_spin,
                                        StagGamma::StagAlgebra, gamma_sink_taste,
                                        std::string, src_gamma,
                                        std::vector<Complex>, corr,
                                        std::vector<std::vector<Complex> >, srcCorrs,
                                        std::vector<Integer>, timeShifts,
                                        Real,                 scaling);
    };
public:
    // constructor
    TStagMeson(const std::string name);
    // destructor
    virtual ~TStagMeson(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    template<typename TField>
    EnableIf<is_lattice<TField>,void> contract(Result &result, const TField &source, const TField &sink, StagGamma& gamma);
    template<typename TField>
    EnableIf<is_lattice<TField>,void> contract(Result &result, const std::vector<TField> &source, const std::vector<TField> &sink, StagGamma& gamma);
    template<typename TField>
    void executeHelper(std::vector<Result> &results, const TField &sink);

    inline void buildProp(PropagatorField &result, const FermionField &source, const FermionField &sink) {
        result = outerProduct(sink,source);
    }
    inline void buildProp(PropagatorField &result, const PropagatorField &source, const PropagatorField &sink) {
        result = sink*adj(source);
    }

    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::string sinkSuffix_ = "";
    std::vector<std::string> sourceGammas_;
    std::vector<StagGamma::SpinTastePair> sinkGammas_;
    Integer Nt_;
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
    std::vector<std::string> in = {par().sinkFunc};

    if(!par().sourceShift.empty()) {
        in.push_back(par().sourceShift);
    }
    
    if (!par().sinkSpinTaste.gauge.empty()) {
        in.push_back(par().sinkSpinTaste.gauge);
    }

    if (!par().sourceGammas.empty()) {
        for (auto gamma : strToVec<StagGamma::SpinTastePair>(par().sourceGammas)) {
            in.push_back(par().source+StagGamma::GetName(gamma));
        }
    } else {
        in.push_back(par().source);
    }

    std::string identityName = StagGamma::GetName(StagGamma::StagAlgebra::G1,StagGamma::StagAlgebra::G1);

    if (env().hasObject(par().sink + identityName)) {
        sinkSuffix_ = identityName;
        in.push_back(par().sink+identityName);
    } else {
        in.push_back(par().sink);
    }

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

    envTmpLat(PropagatorField, "prop");

    Nt_ = env().getDim(Tp);

    sourceGammas_.clear();
    for (auto gamma : strToVec<StagGamma::SpinTastePair>(par().sourceGammas)) {
        sourceGammas_.push_back(StagGamma::GetName(gamma));
    }

    if (!par().sinkSpinTaste.gammas.empty()) {
        sinkGammas_ = strToVec<StagGamma::SpinTastePair>(par().sinkSpinTaste.gammas);

        if (par().sinkSpinTaste.applyG5) {
            StagGamma st;
            StagGamma g5(StagGamma::StagAlgebra::G5,StagGamma::StagAlgebra::G5);
            for (auto &g : sinkGammas_) {
                st.setSpinTaste(g);
                st = st*g5;
                g.first = st._spin;
                g.second = st._taste;
            }
        }
    } else {
        sinkGammas_.push_back(StagGamma::SpinTastePair(StagGamma::StagAlgebra::G1,StagGamma::StagAlgebra::G1));
    }

    if (sourceGammas_.size() > 0 && sourceGammas_.size() != sinkGammas_.size()) {
        HADRONS_ERROR(Argument,"Parameter 'sourceGammas' must be empty or have the same number of operators as 'sinkSpinTaste.gammas'.");
    }
    if (envHasType(PropagatorField, par().sink+sinkSuffix_) || envHasType(std::vector<PropagatorField>, par().sink+sinkSuffix_)) {
        envTmpLat(PropagatorField, "field");
    } else if (envHasType(FermionField, par().sink+sinkSuffix_) || envHasType(std::vector<FermionField>, par().sink+sinkSuffix_)) {
        envTmpLat(FermionField, "field");
    } else {
        HADRONS_ERROR(Argument,"Sink parameter '" + par().sink+"' must be a PropagatorField, FermionField, or a std::vector of these fields.");

    }

}

// execution ///////////////////////////////////////////////////////////////////
template<typename FImpl>
template<typename TField>
EnableIf<is_lattice<TField>,void> TStagMeson<FImpl>::contract(Result &result, const TField &source, const TField &sink, StagGamma& gamma) {

    int offset;
    std::vector<TComplex> buf;

    SinkFnScalar &sinkFunc = envGet(SinkFnScalar, par().sinkFunc);

    envGetTmp(PropagatorField, prop);
    envGetTmp(TField, field);

    gamma(field,sink);

    buildProp(prop, source, field);


    buf = sinkFunc(trace(prop));

    offset = 0;
    if (result.timeShifts.size() > 0) {
        LOG(Message) << "Shifting correlator to (t0 = " << result.timeShifts[0] << ")" << std::endl;
        offset = result.timeShifts[0];
    }

    for (unsigned int t = 0; t < Nt_; ++t) {
        auto ct = TensorRemove(buf[offset]);
        result.srcCorrs[0][t] = ct; // Save corr for individual source

        offset = mod(offset + 1,Nt_);
    }
}

template<typename FImpl>
template<typename TField>
EnableIf<is_lattice<TField>,void> TStagMeson<FImpl>::contract(Result &result, const std::vector<TField> &source, const std::vector<TField> &sink, StagGamma& gamma) {

    int offset;
    std::vector<TComplex> buf;

    SinkFnScalar &sinkFunc = envGet(SinkFnScalar, par().sinkFunc);

    result.srcCorrs.resize(sink.size(), std::vector<Complex>(Nt_,0.0));
    result.scaling = sink.size();
 
    envGetTmp(PropagatorField, prop);
    envGetTmp(TField, field);


    for (int i = 0; i < result.srcCorrs.size(); i++) {

        LOG(Message) << "Contracting element i = " << i << "." << std::endl;

        gamma(field,sink[i]);
        buildProp(prop, source[i], field);


        buf = sinkFunc(trace(prop));

        offset = 0;
        if (result.timeShifts.size() > 0) {
            LOG(Message) << "Shifting correlator " << i << " to (t0 = " << result.timeShifts[i] << ")" << std::endl;
            offset = result.timeShifts[i];
        }

        for (unsigned int t = 0; t < Nt_; ++t) {
            auto ct = TensorRemove(buf[offset]);
            result.srcCorrs[i][t] = ct; // Save corr for individual source

            offset = mod(offset + 1,Nt_);
        }
    }
}

template<typename FImpl>
template<typename TField>
void TStagMeson<FImpl>::executeHelper(std::vector<Result> &results, const TField &sink) {

    std::string srcName;
    StagGamma spinTaste;

    if (!par().sinkSpinTaste.gauge.empty()) {
        auto& U = envGet(GaugeField,par().sinkSpinTaste.gauge);
        spinTaste.setGaugeField(U);
    }

    for (int i = 0; i < results.size(); ++i)
    {
        spinTaste.setSpinTaste(results[i].gamma_sink_spin,results[i].gamma_sink_taste);

        LOG(Message) << "Contracting with gamma: " << spinTaste.getName() << std::endl;

        srcName = par().source;

        if (!par().sourceGammas.empty()) {
            srcName += sourceGammas_[i];
            results[i].src_gamma = sourceGammas_[i];
            LOG(Message) << "Using source gamma: '" << results[i].src_gamma << "'." << std::endl;
        }

        auto & source  = envGet(TField,srcName);

        contract(results[i],source,sink,spinTaste);

        for (int j = 0; j < results[i].srcCorrs.size(); j++) {
            for (int t = 0; t < Nt_; t++) {
                results[i].corr[t] += (results[i].srcCorrs[j])[t]/results[i].scaling;
            }
        }
    }
}

template <typename FImpl>
void TStagMeson<FImpl>::execute(void)
{

    LOG(Message) << "Computing meson contractions '" << getName() << "' using"
                 << " quarks '" << par().source << "' and '" << par().sink << "'"
                 << std::endl;

    std::vector<Result> results;

   results.resize(sinkGammas_.size());

    for (unsigned int i = 0; i < results.size(); ++i)
    {
        results[i].gamma_sink_spin = sinkGammas_[i].first;
        results[i].gamma_sink_taste = sinkGammas_[i].second;

        results[i].srcCorrs.resize(1, std::vector<Complex>(Nt_,0.0));
        results[i].corr.resize(Nt_, 0.0);
        results[i].scaling = 1.0;
        if (!par().sourceShift.empty()) {
            results[i].timeShifts = envGet(std::vector<Integer>, par().sourceShift);
        }
    }

    if (envHasType(PropagatorField, par().sink+sinkSuffix_)) {
        auto &sink  = envGet(PropagatorField, par().sink+sinkSuffix_);
        executeHelper(results,sink);

    } else if (envHasType(std::vector<PropagatorField>, par().sink+sinkSuffix_)) {
        auto &sink  = envGet(std::vector<PropagatorField>, par().sink+sinkSuffix_);
        executeHelper(results,sink);

    } else if (envHasType(FermionField, par().sink+sinkSuffix_)) {
        auto &sink  = envGet(FermionField, par().sink+sinkSuffix_);
        executeHelper(results,sink);

    } else if (envHasType(std::vector<FermionField>, par().sink+sinkSuffix_)) {
        auto &sink  = envGet(std::vector<FermionField>, par().sink+sinkSuffix_);
        executeHelper(results,sink);
    }

    saveResult(par().output, "meson", results);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_StagMeson_hpp_
