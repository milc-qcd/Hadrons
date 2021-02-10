/*
 * Meson.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Vera Guelpers <vmg1n14@soton.ac.uk>
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

#ifndef Hadrons_MContraction_MesonF_hpp_
#define Hadrons_MContraction_MesonF_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Meson contractions
 -----------------------------
 
 * options:
 - q1: input propagator 1 (string)
 - q2: input propagator 2 (string)
 - gammas: gamma products to insert at sink & source, pairs of gamma matrices 
           (space-separated strings) in round brackets (i.e. (g_sink g_src)),
           in a sequence (e.g. "(Gamma5 Gamma5)(Gamma5 GammaT)").

           Special values: "all" - perform all possible contractions.
 - sink: module to compute the sink to use in contraction (string).
*/

/******************************************************************************
 *                                TMeson                                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

typedef std::pair<Gamma::Algebra, Gamma::Algebra> GammaPair;

class MesonFPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFPar,
                                    std::string, q1,
                                    std::string, q2,
                                    std::string, gammas,
                                    std::string, sink,
                                    std::string, output);
};

template <typename FImpl1, typename FImpl2>
class TStagMesonF: public Module<MesonFPar>
{
public:
    FERM_TYPE_ALIASES(FImpl1, 1);
    FERM_TYPE_ALIASES(FImpl2, 2);
    BASIC_TYPE_ALIASES(ScalarImplCR, Scalar);
    SINK_TYPE_ALIASES(Scalar);

    class Result: Serializable {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        Gamma::Algebra, gamma_snk,
                                        Gamma::Algebra, gamma_src,
                                        std::vector<std::vector<Complex>>, corr);
    };
public:
    // constructor
    TStagMesonF(const std::string name);
    // destructor
    virtual ~TStagMesonF(void) {};
    // dependencies/products
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void parseGammaString();
protected:
    void contractionHelper(std::vector<Complex>&,
        const FermionField1&, const FermionField2&, const LatticeComplex&);
    // execution
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<Gamma::Algebra>        gammaList;
};

MODULE_REGISTER_TMP(StagMesonF, ARG(TStagMesonF<STAGIMPL, STAGIMPL>), MContraction);


/******************************************************************************
 *                           TStagMesonF implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
TStagMesonF<FImpl1, FImpl2>::TStagMesonF(const std::string name)
: Module<MesonFPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonF<FImpl1, FImpl2>::getInput(void)
{
    std::vector<std::string> input = {par().q1, par().q2, par().sink};
    
    return input;
}

template <typename FImpl1, typename FImpl2>
std::vector<std::string> TStagMesonF<FImpl1, FImpl2>::getOutput(void)
{
    std::vector<std::string> output = {};
    
    return output;
}

template <typename FImpl1, typename FImpl2>
void TStagMesonF<FImpl1, FImpl2>::parseGammaString()
{
    gammaList.clear();
    
    // Determine gamma matrices to insert at source/sink.
    // only gamma-src = gamma-snk supported for now
    if (par().gammas.compare("all") == 0)
    {
        // Do diagonal contractions, g5,gi(i=1,2,3).
        gammaList.push_back(Gamma::Algebra::Gamma5);
        gammaList.push_back(Gamma::Algebra::GammaX);
        gammaList.push_back(Gamma::Algebra::GammaY);
        gammaList.push_back(Gamma::Algebra::GammaZ);
    }
    else
    {
        // Parse individual contractions from input string.
        gammaList = strToVec<Gamma::Algebra>(par().gammas);
    }
}

// setup ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonF<FImpl1, FImpl2>::setup(void)
{
    envTmpLat(LatticeComplex, "c");

    parseGammaString();

    int Ngammas = gammaList.size();

    if ((envHasType(FermionField1, par().q1) and envHasType(FermionField2, par().q2)) or
        (envHasType(std::vector<FermionField1>, par().q1) and envHasType(std::vector<FermionField2>, par().q2))) {

        envTmp(std::vector<Result>,  "result", 1, Ngammas, Result());

    } else {
        HADRONS_ERROR(Argument, "incompatible field types '" 
            + env().getObjectType(par().q1) + "', and '" + env().getObjectType(par().q2)+ "'.")  
    }

    // grid can't handle real * prop, so use complex
    envTmp(std::vector<LatticeComplex>,  "stag_phase_sink", 1, Ngammas,
           LatticeComplex(env().getGrid()));
    
    envGetTmp(std::vector<LatticeComplex>,stag_phase_sink);

    Lattice<iScalar<vInteger> > x(env().getGrid()); LatticeCoordinate(x,0);
    Lattice<iScalar<vInteger> > y(env().getGrid()); LatticeCoordinate(y,1);
    Lattice<iScalar<vInteger> > z(env().getGrid()); LatticeCoordinate(z,2);
    
    // local taste non-singlet ops, including ``Hermiticity" phase,
    // see Tab. 11.2 in Degrand and Detar
    for(int i=0; i < Ngammas; i++){

        stag_phase_sink[i] = 1.0;

        LOG(Message) << "Using gamma: " << gammaList[i] << std::endl;

        switch(gammaList[i]) {
                
            case Gamma::Algebra::GammaX  :
                stag_phase_sink[i] = where( mod(x,2)==(Integer)0, stag_phase_sink[i], -stag_phase_sink[i]);
                break;
                
            case Gamma::Algebra::GammaY  :
                stag_phase_sink[i] = where( mod(y,2)==(Integer)0, stag_phase_sink[i], -stag_phase_sink[i]);
                break;
                
            case Gamma::Algebra::GammaZ  :
                stag_phase_sink[i] = where( mod(z,2)==(Integer)0, stag_phase_sink[i], -stag_phase_sink[i]);
                break;

            case Gamma::Algebra::Gamma5  :
                break;

            default :
                HADRONS_ERROR(Implementation, "requested gamma operation not implemented")
        }
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl1, typename FImpl2>
void TStagMesonF<FImpl1, FImpl2>::contractionHelper(std::vector<Complex> &corr, 
    const FermionField1 &f1, const FermionField2 &f2, const LatticeComplex &phaseSink) {

    LOG(Message) << "(using sink '" << par().sink << "')" << std::endl;

    std::string ns = vm().getModuleNamespace(env().getObjectModule(par().sink));

    
    if (ns == "MSink") {

        envGetTmp(LatticeComplex, c);

        SinkFnScalar& sink = envGet(SinkFnScalar, par().sink);

        c = localInnerProduct<typename FermionField1::vector_object>(f2,f1);

        c *= phaseSink;

        std::vector<TComplex>  buf = sink(c);

        for (unsigned int i=0; i<corr.size();i++) {
            corr[i] = TensorRemove(buf[i]);
        }

    } else {

        FermionField1 temp = f1*(phaseSink);

        sliceInnerProductVector<typename FermionField1::vector_object>(corr,temp,f2,Tp);
    }
}

template <typename FImpl1, typename FImpl2>
void TStagMesonF<FImpl1, FImpl2>::execute(void)
{
    LOG(Message) << "Computing meson contractions '" << getName() << "' using"
    << " quarks '" << par().q1 << "' and '" << par().q2 << "'"
    << std::endl;
    
    int nt = env().getDim(Tp);

    // staggered gammas
    envGetTmp(std::vector<LatticeComplex>, stag_phase_sink);
    envGetTmp(std::vector<Result>,result);

    for (unsigned int i=0;i<result.size();i++) {

        result[i].gamma_snk = gammaList[i];
        result[i].gamma_src = Gamma::Algebra::Identity;

        if (envHasType(FermionField1, par().q1)) {

            auto& q1 = envGet(FermionField1, par().q1);
            auto& q2 = envGet(FermionField2, par().q2);

            result[i].corr.resize(1,std::vector<Complex>(nt));

            contractionHelper(result[i].corr[0],q1,q2,stag_phase_sink[i]);

        } else {

            auto& q1 = envGet(std::vector<FermionField1>, par().q1);
            auto& q2 = envGet(std::vector<FermionField2>, par().q2);

            if (q1.size() != q2.size()) {
                HADRONS_ERROR(Implementation, "Propagators " + par().q1 + " and " 
                    + par().q2 + " must have the same number of vectors")
            }

            result[i].corr.resize(q1.size(),std::vector<Complex>(nt));

            for (unsigned int j=0;j<q1.size();j++){
                contractionHelper(result[i].corr[j],q1[j],q2[j],stag_phase_sink[i]);
            }
        }
    }
    saveResult(par().output, "meson", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Meson_hpp_
