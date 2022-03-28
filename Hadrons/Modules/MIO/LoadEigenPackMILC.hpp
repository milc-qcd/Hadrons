/*
 * LoadEigenPackMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MIO_LoadEigenPackMILC_hpp_
#define Hadrons_MIO_LoadEigenPackMILC_hpp_

#include <typeinfo>
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Load eigen vectors/values package                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadEigenPackMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadEigenPackMILCPar,
                                    std::string,  filestem,
                                    bool,         multiFile,
                                    unsigned int, size,
                                    unsigned int, Ls,
                                    std::string,  gaugeXform,
                                    bool,         evenEigen,
                                    double,       mass);
};

template <typename Pack, typename GImpl>
class TLoadEigenPackMILC: public Module<LoadEigenPackMILCPar>
{
public:
    typedef typename Pack::Field   Field;
    typedef typename Pack::FieldIo FieldIo;
    typedef typename Pack::EvalType EvalType;
    typedef BaseEigenPack<Field,EvalType>   BasePack;

public:
    GAUGE_TYPE_ALIASES(GImpl, );
    typedef typename GImpl::GaugeLinkField GaugeMat;
public:
    // constructor
    TLoadEigenPackMILC(const std::string name);
    // destructor
    virtual ~TLoadEigenPackMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadFermionEigenPackMILC, ARG(TLoadEigenPackMILC<FermionEigenPack<STAGIMPL,STAGIMPL,Complex>, GIMPL>), MIO);
MODULE_REGISTER_TMP(LoadFermionEigenPackMILCHermitian, ARG(TLoadEigenPackMILC<FermionEigenPack<STAGIMPL>, GIMPL>), MIO);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(LoadFermionEigenPackMILCF, ARG(TLoadEigenPackMILC<FermionEigenPack<STAGIMPLF, STAGIMPLF, Complex>, GIMPLF>), MIO);
MODULE_REGISTER_TMP(LoadFermionEigenPackMILCIo32, ARG(TLoadEigenPackMILC<FermionEigenPack<STAGIMPL, STAGIMPLF, Complex>, GIMPL>), MIO);
MODULE_REGISTER_TMP(LoadFermionEigenPackMILCHermitianF, ARG(TLoadEigenPackMILC<FermionEigenPack<STAGIMPLF>, GIMPLF>), MIO);
MODULE_REGISTER_TMP(LoadFermionEigenPackMILCHermitianIo32, ARG(TLoadEigenPackMILC<FermionEigenPack<STAGIMPL, STAGIMPLF>, GIMPL>), MIO);
#endif

/******************************************************************************
 *                    TLoadEigenPackMILC implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
TLoadEigenPackMILC<Pack, GImpl>::TLoadEigenPackMILC(const std::string name)
: Module<LoadEigenPackMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
std::vector<std::string> TLoadEigenPackMILC<Pack, GImpl>::getInput(void)
{
    std::vector<std::string> in;

    if (!par().gaugeXform.empty())
    {
        in = {par().gaugeXform};
    }
    
    return in;
}

template <typename Pack, typename GImpl>
std::vector<std::string> TLoadEigenPackMILC<Pack, GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName() + "_evec", getName() + "_eval"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
void TLoadEigenPackMILC<Pack, GImpl>::setup(void)
{
    GridBase *gridIo = nullptr;

    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        gridIo = envGetRbGrid(FieldIo, par().Ls);
    }

    envCreate(std::vector<Field>,getName() + "_evec", par().Ls, par().size, envGetRbGrid(Field, par().Ls));
    envCreate(std::vector<EvalType>,getName() + "_eval", par().Ls, par().size);

    auto &evecOut = envGet(std::vector<Field>,getName() + "_evec");
    auto &evalOut = envGet(std::vector<EvalType>,getName() + "_eval");

    envCreateDerived(BasePack, Pack, getName(), par().Ls, evecOut, evalOut, gridIo);

    if (!par().gaugeXform.empty())
    {
        if (par().Ls > 1)
        {
            LOG(Message) << "Setup 5d GaugeMat for Ls = " << par().Ls << std::endl;
            envTmp(GaugeMat,    "tmpXform", par().Ls, envGetGrid5(Field, par().Ls));
            envTmp(GaugeMat, "tmpXformOdd", par().Ls, envGetRbGrid5(Field, par().Ls));
        }
        else
        {
            LOG(Message) << "Setup 4d GaugeMat for Ls = " << par().Ls << std::endl;
            envTmp(GaugeMat,    "tmpXform", par().Ls, envGetGrid(Field));
            envTmp(GaugeMat, "tmpXformOdd", par().Ls, envGetRbGrid(Field));
        }
        
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
void TLoadEigenPackMILC<Pack, GImpl>::execute(void)
{
    auto &epack = envGetDerived(BasePack, Pack, getName());

    epack.read(par().filestem, par().multiFile, vm().getTrajectory());
    epack.eval.resize(par().size);

    if (par().mass > 0.0) {
        Real m = 2*par().mass;

        if (typeid(decltype(epack.eval[0])) == typeid(Real)) {
            m = ::pow(m,2);
            epack.record.operatorXml = "<!-- WARNING! This EigenPack has been altered! metadata may be inaccurate. Shifted eigenvalues by m^2; m = " 
                + std::to_string(2*par().mass) + ". -->" + epack.record.operatorXml; 
        } else {
            LOG(Warning) << "The LoadEigenPackMILC module provides eigenvalues of the Dirac operator, i.e. mass + i*lambda_D." << std::endl;

            epack.record.operatorXml = "<!-- WARNING! This EigenPack has been altered! metadata may be inaccurate. Changed evals to m+i*lambda_D; m = " 
                + std::to_string(2*par().mass) + ". -->" + epack.record.operatorXml; 

        }
        LOG(Message) << "Shifting eigenvalues by mass (including MILC factor of 2) = " << m << std::endl;

        for (auto &lam:epack.eval) {
            lam += m;
        }        
        LOG(Message) << epack.record.operatorXml << std::endl;
    }


    ComplexD norm(1.0,0.0);
    if (typeid(decltype(epack.eval[0])) != typeid(Real)) {
        norm *= 1.0/sqrt(2.0);
        LOG(Message) << "Normalizing eigenvectors by 1/sqrt(2)" << std::endl;
    }
    for (auto &e:epack.evec) {
        e *= norm;
        e.Checkerboard() = (par().evenEigen ? Even : Odd);
    }

    if (!par().gaugeXform.empty())
    {

        LOG(Message) << "Applying gauge transformation to eigenvectors " << getName()
                     << " using " << par().gaugeXform << std::endl;
        auto &xform = envGet(GaugeMat, par().gaugeXform);
        envGetTmp(GaugeMat,    tmpXform);
        envGetTmp(GaugeMat, tmpXformOdd);

        if (par().Ls > 1) 
        {
            LOG(Message) << "Creating 5d GaugeMat from " << par().gaugeXform << std::endl;
            startTimer("5-d gauge transform creation");
            for (unsigned int j = 0; j < par().Ls; j++)
            {
                InsertSlice(xform, tmpXform, j, 0);
            }
            stopTimer("5-d gauge transform creation");
        }
        else
        {
            tmpXform = xform;
        }

        pickCheckerboard(Odd, tmpXformOdd, tmpXform);
        startTimer("Transform application");
        for (unsigned int i = 0; i < par().size; i++)
        {
            LOG(Message) << "Applying gauge transformation to eigenvector i = " << i+1 << "/" << par().size << std::endl;
            epack.evec[i].Checkerboard() = Odd;
            epack.evec[i] = tmpXformOdd * epack.evec[i];
        }
        stopTimer("Transform application");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadEigenPackMILC_hpp_
