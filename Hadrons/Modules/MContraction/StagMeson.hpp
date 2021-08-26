/*
 * StagMeson.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
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

class StagMesonPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(StagMesonPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, gammas,
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
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma,
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
    template<typename TFImpl = FImpl, IfStag<TFImpl> = 0>
    void parseGammaString();
    template<typename TFImpl = FImpl, IfNotStag<TFImpl> = 0>
    void parseGammaString();
protected:
    template<typename TField,  EnableIf<is_lattice<TField>, int> = 0 >
    void contract(std::vector<Result> &ret, const TField &fSink, const TField &fSrc, Real scale = 1.0, Integer shift = 0);
    template<typename TField,  EnableIf<is_lattice<TField>, int> = 0 >
    void contract(std::vector<Result> &ret, const std::vector<TField> &fSink, const std::vector<TField> &fSrc);
    inline void buildProp(PropagatorField &ret, const FermionField &fSink, const FermionField &fSrc) {
        ret = outerProduct(fSink,fSrc);
    }
    inline void buildProp(PropagatorField &ret, const PropagatorField &fSink, const PropagatorField &fSrc) {
        ret = fSink*adj(fSrc);
    }
    template<typename TFImpl = FImpl, IfStag<TFImpl> = 0>
    inline void applyGamma(PropagatorField &ret, int gamma) {
        ret = stag_phase_sink_[gamma]*ret;
    }
    template<typename TFImpl = FImpl, IfNotStag<TFImpl> = 0>
    inline void applyGamma(PropagatorField &ret, int gamma) {
        Gamma g(gammaList_[gamma]);
        ret = g*ret;
    }
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<LatticeComplex> stag_phase_sink_;
    std::vector<Gamma::Algebra>        gammaList_;
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
    std::vector<std::string> in = {par().q1, par().q2, par().sink};
    if(!par().sourceShift.empty())
        in.push_back(par().sourceShift+"_shift");
    
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
    envTmpLat(PropagatorField, "op");
}

template <typename FImpl>
template <typename TFImpl, IfStag<TFImpl> >
void TStagMeson<FImpl>::parseGammaString()
{
    gammaList_.clear();
    // Determine gamma matrices to insert at source/sink.
    if (par().gammas.compare("all") == 0)
    {
        // Do all contractions.
        for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
        {
            gammaList_.push_back((Gamma::Algebra)i);
        }
    }
    else
    {
        // Parse individual contractions from input string.
        gammaList_ = strToVec<Gamma::Algebra>(par().gammas);
    } 

    int Ngam=gammaList_.size();

    stag_phase_sink_.resize(Ngam, env().getGrid());

    Lattice<iScalar<vInteger> > x(env().getGrid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(env().getGrid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(env().getGrid()); LatticeCoordinate(z,2);
    
    // local taste non-singlet ops, including ``Hermiticity" phase,
    // see Tab. 11.2 in Degrand and Detar
    for(int i=0; i < gammaList_.size(); i++){

        stag_phase_sink_[i] = 1.0;
        
        LOG(Message) << "Using gamma: " << gammaList_[i] << std::endl;
        switch(gammaList_[i]) {
                
            case Gamma::Algebra::GammaX  :
                stag_phase_sink_[i] = where( mod(x,2)==(Integer)0, stag_phase_sink_[i], -stag_phase_sink_[i]);
                break;
                
            case Gamma::Algebra::GammaY  :
                stag_phase_sink_[i] = where( mod(y,2)==(Integer)0, stag_phase_sink_[i], -stag_phase_sink_[i]);
                break;
                
            case Gamma::Algebra::GammaZ  :
                stag_phase_sink_[i] = where( mod(z,2)==(Integer)0, stag_phase_sink_[i], -stag_phase_sink_[i]);
                break;

            case Gamma::Algebra::Gamma5  :
                break;

            default :
                std::cout << "your gamma is not supported for stag fields" << std::endl;
                assert(0);
        }
    }
}

template <typename FImpl>
template <typename TFImpl, IfNotStag<TFImpl> >
void TStagMeson<FImpl>::parseGammaString()
{
    gammaList_.clear();
    // Determine gamma matrices to insert at source/sink.
    if (par().gammas.compare("all") == 0)
    {
        // Do all contractions.
        for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
        {
            gammaList_.push_back((Gamma::Algebra)i);
        }
    }
    else
    {
        // Parse individual contractions from input string.
        gammaList_ = strToVec<Gamma::Algebra>(par().gammas);
    } 
}

// execution ///////////////////////////////////////////////////////////////////
template<typename FImpl>
template<typename TField, EnableIf<is_lattice<TField>, int> >
void TStagMeson<FImpl>::contract(std::vector<Result> &ret, const TField &fSink, const TField &fSrc, Real scale, Integer shift) {

    int offset, nt = env().getDim(Tp);
    std::vector<TComplex>     buf;

    SinkFnScalar &sink = envGet(SinkFnScalar, par().sink);

    envGetTmp(PropagatorField, op);

    buildProp(op, fSink,fSrc);

    for (unsigned int g = 0; g < gammaList_.size(); ++g)
    {
        applyGamma(op,g);

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
template<typename TField, EnableIf<is_lattice<TField>, int> >
void TStagMeson<FImpl>::contract(std::vector<Result> &ret, const std::vector<TField> &fSink, const std::vector<TField> &fSrc) {

    std::vector<Integer> shifts;
    if (!par().sourceShift.empty()) {
        shifts = envGet(std::vector<Integer>, par().sourceShift+"_shift");
    }

    for (int i = 0; i < fSink.size(); i++) {
        if (!par().sourceShift.empty())
            contract(ret,fSink[i],fSrc[i],fSink.size(), shifts[i]);
        else
            contract(ret,fSink[i],fSrc[i],fSink.size());
    }


}

template <typename FImpl>
void TStagMeson<FImpl>::execute(void)
{

    std::vector<Result> result;
    unsigned int nt = env().getDim(Tp);

    parseGammaString();

    result.resize(gammaList_.size());

    for (unsigned int g = 0; g < result.size(); ++g)
    {
        result[g].gamma = gammaList_[g];
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

#endif // Hadrons_MContraction_StagMeson_hpp_
