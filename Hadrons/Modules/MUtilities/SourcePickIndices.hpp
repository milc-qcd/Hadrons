/*
 * SourcePickIndices.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
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
#ifndef Hadrons_MUtilities_SourcePickIndices_hpp_
#define Hadrons_MUtilities_SourcePickIndices_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                Utility module to unpack a vector of fields                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class SourcePickIndicesPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SourcePickIndicesPar,
                                    std::string,  source,
                                    std::string,  indices);
};

template <typename Field>
class TSourcePickIndices: public Module<SourcePickIndicesPar>
{
public:
    // constructor
    TSourcePickIndices(const std::string name);
    // destructor
    virtual ~TSourcePickIndices(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    DEFINE_ENV_ALIAS;
    std::vector<Integer> indices_;
    bool includeShifts_{false};
};

MODULE_REGISTER_TMP(StagComplexSourcePickIndices, TSourcePickIndices<STAGIMPL::ComplexField>, MUtilities);
MODULE_REGISTER_TMP(StagFermionSourcePickIndices, TSourcePickIndices<STAGIMPL::FermionField>, MUtilities);
MODULE_REGISTER_TMP(StagSourcePickIndices, TSourcePickIndices<STAGIMPL::PropagatorField>, MUtilities);

/******************************************************************************
 *                       TSourcePickIndices implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TSourcePickIndices<Field>::TSourcePickIndices(const std::string name)
: Module<SourcePickIndicesPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TSourcePickIndices<Field>::getInput(void)
{
    std::vector<std::string> in = {par().source};

    if (env().hasObject(par().source+"_shift")) {
        in.push_back(par().source+"_shift");
        includeShifts_ = true;
    }
    
    return in;
}

template <typename Field>
std::vector<std::string> TSourcePickIndices<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    if (env().hasObject(par().source+"_shift")) {
        out.push_back(getName()+"_shift");
    }

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TSourcePickIndices<Field>::setup(void)
{
    indices_ = strToVec<Integer>(par().indices);

    envCreate(std::vector<Field>, getName(), 1, indices_.size(), envGetGrid(Field));
    if (includeShifts_) {
        envCreate(std::vector<Integer>, getName()+"_shift", 1, indices_.size(),0);
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TSourcePickIndices<Field>::execute(void)
{
    std::vector<Field> *src;
    if (envHasType(std::vector<Field>,par().source)) {
        auto &vec = envGet(std::vector<Field>, par().source);
        src = &vec;
    } else if (envHasType(ARG(std::map<Gamma::Algebra,std::vector<Field> >),par().source)) {
        auto &vec = envGet(ARG(std::map<Gamma::Algebra,std::vector<Field> >), par().source);
        src = &(vec.at(Gamma::Algebra::Gamma5));
    }

    auto &out = envGet(std::vector<Field>,getName());

    std::vector<Integer> *shift_in, *shift_out;
    if (includeShifts_) {
        auto &vec1 = envGet(std::vector<Integer>,par().source+"_shift");
        shift_in = &vec1;
        auto &vec2 = envGet(std::vector<Integer>,getName()+"_shift");
        shift_out = &vec2;
    }

    for (int i=0;i<indices_.size();i++) {
        LOG(Message) << "Adding source '" << indices_[i] << " of " << par().source << " to " << getName() << std::endl;
        out[i] = src->at(indices_[i]);
        if (includeShifts_) {
            shift_out->at(i) = shift_in->at(indices_[i]);
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_SourcePickIndices_hpp_
