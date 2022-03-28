/*
 * LoadCoarseEigenPackMILC.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MIO_LoadCoarseEigenPackMILC_hpp_
#define Hadrons_MIO_LoadCoarseEigenPackMILC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Load local coherence eigen vectors/values package             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadCoarseEigenPackMILCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadCoarseEigenPackMILCPar,
                                    std::string, filestem,
                                    bool,         multiFile,
                                    unsigned int, sizeFine,
                                    unsigned int, sizeCoarse,
                                    unsigned int, Ls,
                                    bool,         evenEigen,
                                    double,       mass,
                                    std::vector<int>, blockSize);
};

template <typename Pack>
class TLoadCoarseEigenPackMILC: public Module<LoadCoarseEigenPackMILCPar>
{
public:
    typedef typename Pack::Field                Field;
    typedef typename Pack::FieldIo              FieldIo;
    typedef typename Pack::CoarseField          CoarseField;
    typedef typename Pack::CoarseFieldIo        CoarseFieldIo;
    typedef CoarseEigenPack<Field, CoarseField, FieldIo, CoarseFieldIo> BasePack;
    template <typename vtype> 
    using iImplScalar = iScalar<iScalar<iScalar<vtype>>>;
    typedef iImplScalar<typename Pack::Field::vector_type> SiteComplex;
public:
    // constructor
    TLoadCoarseEigenPackMILC(const std::string name);
    // destructor
    virtual ~TLoadCoarseEigenPackMILC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadCoarseFermionEigenPackMILC, 
                    ARG(TLoadCoarseEigenPackMILC<CoarseFermionEigenPack<STAGIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS>>), MIO);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(LoadCoarseFermionEigenPackMILCIo32, 
                    ARG(TLoadCoarseEigenPackMILC<CoarseFermionEigenPack<STAGIMPL, HADRONS_DEFAULT_LANCZOS_NBASIS, STAGIMPLF>>), MIO);
#endif

/******************************************************************************
 *                 TLoadCoarseEigenPackMILC implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack>
TLoadCoarseEigenPackMILC<Pack>::TLoadCoarseEigenPackMILC(const std::string name)
: Module<LoadCoarseEigenPackMILCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack>
std::vector<std::string> TLoadCoarseEigenPackMILC<Pack>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename Pack>
std::vector<std::string> TLoadCoarseEigenPackMILC<Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName()+"_mass", getName() + "_evenEigen"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack>
void TLoadCoarseEigenPackMILC<Pack>::setup(void)
{
    GridBase     *gridIo = nullptr, *gridCoarseIo = nullptr;

    envCreate(std::vector<Real>, getName()+"_mass", 1, 1, 2.*par().mass);
    envCreate(std::vector<bool>, getName()+"_evenEigen", 1, 1, par().evenEigen == true);

    if (par().mass > 0) {
        LOG(Warning) << "The LoadEigenPackMILC module assumes MASSLESS eigenvalues of the Dirac Operator squarred." << std::endl;
    }

    if (typeHash<Field>() != typeHash<FieldIo>())
    {
        gridIo = (par().Ls>1?envGetRbGrid(FieldIo, par().Ls):envGetRbGrid(FieldIo));
    }
    if (typeHash<CoarseField>() != typeHash<CoarseFieldIo>())
    {
        gridCoarseIo = envGetCoarseGrid(CoarseFieldIo, par().blockSize, par().Ls);
    }
    envCreate(std::vector<Field>,getName() + "_evec_fine", par().Ls, par().sizeFine, envGetRbGrid(Field, par().Ls));
    envCreate(std::vector<RealD>,getName() + "_eval_fine", par().Ls, par().sizeFine);
    envCreate(std::vector<CoarseField>,getName() + "_evec_coarse", par().Ls, par().sizeCoarse, envGetCoarseGrid(CoarseField, par().blockSize, par().Ls));
    envCreate(std::vector<RealD>,getName() + "_eval_coarse", par().Ls, par().sizeCoarse);

    auto &evecOut       = envGet(std::vector<Field>,getName() + "_evec_fine");
    auto &evalOut       = envGet(std::vector<RealD>,getName() + "_eval_fine");
    auto &evecCoarseOut = envGet(std::vector<CoarseField>,getName() + "_evec_coarse");
    auto &evalCoarseOut = envGet(std::vector<RealD>,getName() + "_eval_coarse");

    envCreateDerived(BasePack, Pack, getName(), par().Ls, evecOut, evalOut, evecCoarseOut, evalCoarseOut, gridIo, gridCoarseIo);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack>
void TLoadCoarseEigenPackMILC<Pack>::execute(void)
{
    auto                 cg     = envGetCoarseGrid(CoarseField, par().blockSize, par().Ls);
    auto                 &epack = envGetDerived(BasePack, Pack, getName());
    Lattice<SiteComplex> dummy(cg);

    epack.read(par().filestem, par().multiFile, vm().getTrajectory());

    if (par().mass > 0.0) {
        Real m2 = pow(2*par().mass,2);

        LOG(Message) << "Shifting eigenvalues by mass^2 (including MILC factor of 2) = " << m2 << std::endl;

        for (auto &lam:epack.eval) {
            lam += m2;
        }        

        epack.record.operatorXml = "<!-- WARNING! This EigenPack has been altered! metadata may be inaccurate. Added m^2 to evals; m = " 
            + std::to_string(2*par().mass) + ". -->" + epack.record.operatorXml; 
        LOG(Message) << epack.record.operatorXml << std::endl;
    }

    LOG(Message) << "Block Gramm-Schmidt pass 1"<< std::endl;
    blockOrthogonalise(dummy, epack.evec);
    LOG(Message) << "Block Gramm-Schmidt pass 2"<< std::endl;
    blockOrthogonalise(dummy, epack.evec);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadCoarseEigenPackMILC_hpp_
