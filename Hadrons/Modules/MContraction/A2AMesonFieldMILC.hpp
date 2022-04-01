/*
 * A2AMesonFieldMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: ferben <ferben@debian.felix.com>
 * Author: paboyle <paboyle@ph.ed.ac.uk>
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
#ifndef Hadrons_MContraction_A2AMesonFieldMILC_hpp_
#define Hadrons_MContraction_A2AMesonFieldMILC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/A2AMatrix.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     All-to-all meson field creation                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AMesonFieldMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMesonFieldMILCPar,
                                    int, cacheBlock,
                                    int, block,
                                    std::string, lowModes,
                                    std::string, left,
                                    std::string, action,
                                    std::string, right,
                                    std::string, output,
                                    std::string, gammas,
                                    std::vector<std::string>, mom);
};

class A2AMesonFieldMILCMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMesonFieldMILCMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra, gamma);
};

template <typename T, typename FImpl>
class MesonFieldMILCKernel: public A2AKernel<T, typename FImpl::FermionField>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    MesonFieldMILCKernel(const std::vector<Gamma::Algebra> &gamma,
                     const std::vector<LatticeComplex> &mom,
                     GridBase *grid)
    : gamma_(gamma), mom_(mom), grid_(grid)
    {
        vol_ = 1.;
        for (auto &d: grid_->GlobalDimensions())
        {
            vol_ *= d;
        }
    }

    virtual ~MesonFieldMILCKernel(void) = default;
    virtual void operator()(A2AMatrixSet<T> &m, const FermionField *left, 
                            const FermionField *right,
                            const unsigned int orthogDim, double &t, double &tg)
    {
        MesonFunction<FImpl>(m, left, right, gamma_, mom_, orthogDim, &t, &tg);
    }

    virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        return vol_*(2*8.0+6.0+8.0*mom_.size())*blockSizei*blockSizej*gamma_.size();
    }

    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej)
    {
        return vol_*(12.0*sizeof(T))*blockSizei*blockSizej
               +  vol_*(2.0*sizeof(T)*mom_.size())*blockSizei*blockSizej*gamma_.size();
    }
private:
 template<typename TFImpl, typename ... Args>
 IfNotStag<TFImpl,void> MesonFunction(Args && ... args){
     A2Autils<FImpl>::MesonField(args...);
 }

 template<typename TFImpl, typename ... Args>
 IfStag<TFImpl,void> MesonFunction(Args && ... args){
     A2Autils<FImpl>::StagMesonFieldMILC(args...);
 }

private:
    const std::vector<Gamma::Algebra> &gamma_;
    const std::vector<LatticeComplex> &mom_;
    GridBase                          *grid_;
    double                            vol_;
};

template <typename FImpl>
class TA2AMesonFieldMILC : public Module<A2AMesonFieldMILCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef A2AMatrixBlockComputation<Complex, 
                                      FermionField, 
                                      A2AMesonFieldMILCMetadata, 
                                      HADRONS_A2AM_IO_TYPE> Computation;
    typedef MesonFieldMILCKernel<Complex, FImpl> Kernel;
public:
    // constructor
    TA2AMesonFieldMILC(const std::string name);
    // destructor
    virtual ~TA2AMesonFieldMILC(void){};

    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool                               hasPhase_{false};
    std::string                        momphName_;
    std::vector<Gamma::Algebra>        gamma_;
    std::vector<std::vector<Real>>     mom_;
};

MODULE_REGISTER(A2AMesonFieldMILC, ARG(TA2AMesonFieldMILC<STAGIMPL>), MContraction);

/******************************************************************************
*                  TA2AMesonFieldMILC implementation                             *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AMesonFieldMILC<FImpl>::TA2AMesonFieldMILC(const std::string name)
: Module<A2AMesonFieldMILCPar>(name)
, momphName_(name + "_momph")
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AMesonFieldMILC<FImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    if (!par().left.empty()) 
        in.push_back(par().left);
    if (!par().right.empty()) 
        in.push_back(par().right);

    if (!par().lowModes.empty()) {
        if (!par().action.empty())
           in.push_back(par().action);
       in.push_back(par().lowModes+"_evec");
       in.push_back(par().lowModes+"_evalM");
    }

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AMesonFieldMILC<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMesonFieldMILC<FImpl>::setup(void)
{
    if (!par().action.empty()) {
    }

    gamma_.clear();
    mom_.clear();
    if (par().gammas == "all")
    {
        gamma_ = {
            Gamma::Algebra::Gamma5,
            Gamma::Algebra::Identity,    
            Gamma::Algebra::GammaX,
            Gamma::Algebra::GammaY,
            Gamma::Algebra::GammaZ,
            Gamma::Algebra::GammaT,
            Gamma::Algebra::GammaXGamma5,
            Gamma::Algebra::GammaYGamma5,
            Gamma::Algebra::GammaZGamma5,
            Gamma::Algebra::GammaTGamma5,
            Gamma::Algebra::SigmaXY,
            Gamma::Algebra::SigmaXZ,
            Gamma::Algebra::SigmaXT,
            Gamma::Algebra::SigmaYZ,
            Gamma::Algebra::SigmaYT,
            Gamma::Algebra::SigmaZT
        };
    }
    else
    {
        gamma_ = strToVec<Gamma::Algebra>(par().gammas);
    }
    for (auto &pstr: par().mom)
    {
        auto p = strToVec<Real>(pstr);

        if (p.size() != env().getNd() - 1)
        {
            HADRONS_ERROR(Size, "Momentum has " + std::to_string(p.size())
                                + " components instead of " 
                                + std::to_string(env().getNd() - 1));
        }
        mom_.push_back(p);
    }
    envCache(std::vector<ComplexField>, momphName_, 1, 
             par().mom.size(), envGetGrid(ComplexField));
    envTmpLat(ComplexField, "coor");
    envTmp(Computation, "computation", 1, envGetGrid(FermionField), 
           env().getNd() - 1, mom_.size(), gamma_.size(), par().block, 
           par().cacheBlock, this);
    envTmp(std::vector<FermionField>, "dummy", 1, 0, envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AMesonFieldMILC<FImpl>::execute(void)
{
    bool hasHighModes = (!par().left.empty() && !par().right.empty());
    bool hasLowModes = (!par().lowModes.empty());
    bool isCheckerBoarded = (!par().action.empty());

    std::vector<FermionField> *left, *right;
    if (hasHighModes) {
        left  = &(envGet(std::vector<FermionField>, par().left));
        right = &(envGet(std::vector<FermionField>, par().right));
    } else {
        envGetTmp(std::vector<FermionField>, dummy);
        left = &dummy;
        right = &dummy;
    }

    int nt         = env().getDim().back();
    int N_i        = left->size();
    int N_j        = right->size();

    if (hasLowModes)
    {
        auto &lowModeVec = envGet(std::vector<FermionField>, par().lowModes+"_evec");
        N_i += (isCheckerBoarded?2:1)*lowModeVec.size();
        N_j += (isCheckerBoarded?2:1)*lowModeVec.size();
    }
    int ngamma     = gamma_.size();
    int nmom       = mom_.size();
    int block      = par().block;
    int cacheBlock = par().cacheBlock;

    if (N_i < block || N_j < block)
    {
        HADRONS_ERROR(Range, "blockSize must not exceed size of input vector.");
    }

    LOG(Message) << "Computing all-to-all meson fields" << std::endl;
    if (hasLowModes)
        LOG(Message) << "Low Modes: '" << par().lowModes << std::endl;
    if (hasHighModes)
        LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right << "'" << std::endl;
    LOG(Message) << "Momenta:" << std::endl;

    for (auto &p: mom_)
    {
        LOG(Message) << "  " << p << std::endl;
    }
    LOG(Message) << "Spin bilinears:" << std::endl;

    for (auto &g: gamma_)
    {
        LOG(Message) << "  " << g << std::endl;
    }

    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j 
                 << " (filesize " << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE)) 
                 << "/momentum/bilinear)" << std::endl;

    auto &ph = envGet(std::vector<ComplexField>, momphName_);

    if (!hasPhase_)
    {
        startTimer("Momentum phases");
        for (unsigned int j = 0; j < nmom; ++j)
        {
            Complex           i(0.0,1.0);
            std::vector<Real> p;

            envGetTmp(ComplexField, coor);
            ph[j] = Zero();
            for(unsigned int mu = 0; mu < mom_[j].size(); mu++)
            {
                LatticeCoordinate(coor, mu);
                ph[j] = ph[j] + (mom_[j][mu]/env().getDim(mu))*coor;
            }
            ph[j] = exp((Real)(2*M_PI)*i*ph[j]);
        }
        hasPhase_ = true;
        stopTimer("Momentum phases");
    }

    auto ionameFn = [this](const unsigned int m, const unsigned int g)
    {
        std::stringstream ss;

        ss << gamma_[g] << "_";
        for (unsigned int mu = 0; mu < mom_[m].size(); ++mu)
        {
            ss << mom_[m][mu] << ((mu == mom_[m].size() - 1) ? "" : "_");
        }

        return ss.str();
    };

    auto filenameFn = [this, &ionameFn](const unsigned int m, const unsigned int g)
    {
        return par().output + "." + std::to_string(vm().getTrajectory()) 
               + "/" + ionameFn(m, g) + ".h5";
    };

    auto metadataFn = [this](const unsigned int m, const unsigned int g)
    {
        A2AMesonFieldMILCMetadata md;

        for (auto pmu: mom_[m])
        {
            md.momentum.push_back(pmu);
        }
        md.gamma = gamma_[g];
        
        return md;
    };

    envGetTmp(Computation, computation);

    Kernel      kernel(gamma_, ph, envGetGrid(FermionField));

    if(hasLowModes) {
        if (isCheckerBoarded) {
            auto &action      = envGet(FMat, par().action);
            auto &lowModeVec = envGet(std::vector<FermionField>, par().lowModes+"_evec");
            auto &lowModeVal = envGet(std::vector<ComplexD>, par().lowModes+"_evalM");

            std::function<void(int)> swapEvecCheckerFn = [this,&action, &lowModeVec, &lowModeVal](int index)
            {
                ComplexD eval_D = ComplexD(0.0,lowModeVal[index].imag());
                int cb = lowModeVec[index].Checkerboard();
                int cbNeg = (cb==Even) ? Odd : Even;

                FermionField temp(lowModeVec[index].Grid());
                temp.Checkerboard() = cbNeg;
                action.Meooe(lowModeVec[index], temp);
                lowModeVec[index].Checkerboard() = cbNeg;
                lowModeVec[index] = (1.0/eval_D) * temp;
            };

            computation.execute(*left, *right, kernel, ionameFn, filenameFn, metadataFn, &lowModeVec, lowModeVal, &swapEvecCheckerFn);
        } else{
            auto &lowModeVec = envGet(std::vector<FermionField>, par().lowModes+"_evec");
            auto &lowModeVal = envGet(std::vector<ComplexD>, par().lowModes+"_evalM");

            computation.execute(*left, *right, kernel, ionameFn, filenameFn, metadataFn, &lowModeVec, lowModeVal);

        }
    } else {
        computation.execute(*left, *right, kernel, ionameFn, filenameFn, metadataFn);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AMesonFieldMILC_hpp_
