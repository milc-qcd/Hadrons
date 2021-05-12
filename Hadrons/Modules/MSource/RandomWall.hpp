/*
 * RandomWall.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: Michael Marshall <43034299+mmphys@users.noreply.github.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: fionnoh <fionnoh@gmail.com>
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

#ifndef Hadrons_MSource_RandomWall_hpp_
#define Hadrons_MSource_RandomWall_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Random Wall source
 -----------------------------
 
 * options:
 - tW:   source timeslice (integer)
 - size: number of sources (integer)
 
 */

/******************************************************************************
 *                         Random Wall                                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class RandomWallPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RandomWallPar,
                                    unsigned int, tW,
                                    unsigned int, size);
};

template <typename FImpl>
class TRandomWall: public Module<RandomWallPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TRandomWall(const std::string name);
    // destructor
    virtual ~TRandomWall(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool        hasT_{false};
    std::string tName_;
};

// MODULE_REGISTER_TMP(RandomWall, TRandomWall<FIMPL>, MSource);
MODULE_REGISTER_TMP(StagRandomWall, TRandomWall<STAGIMPL>, MSource);

/******************************************************************************
 *                 TRandomWall implementation                                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TRandomWall<FImpl>::TRandomWall(const std::string name)
: Module<RandomWallPar>(name)
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TRandomWall<FImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TRandomWall<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandomWall<FImpl>::setup(void)
{
    if (par().size && par().size > 1) {
        envCreate(std::vector<FermionField>, getName(), 1, par().size, envGetGrid(FermionField));
    } else {
        envCreate(std::vector<FermionField>, getName(), 1, 1, envGetGrid(FermionField));
    }
    envCache(Lattice<iScalar<vInteger>>, tName_,    1, envGetGrid(LatticeComplex));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TRandomWall<FImpl>::execute(void)
{    
    typedef typename FermionField::scalar_object scalar_object;
    typedef typename FermionField::scalar_type scalar_type;


    LOG(Message) << "Generating random wall source at t = " << par().tW 
                 << std::endl;
    
    auto  &t   = envGet(Lattice<iScalar<vInteger>>, tName_);
    auto  nc   = FImpl::Dimension;
    auto  &vec = envGet(std::vector<FermionField>, getName());
    
    if (!hasT_)
    {
        LatticeCoordinate(t, Tp);
        hasT_ = true;
    }

    auto &rng      = rng4d();
    GridBase *grid = rng.Grid();

    int multiplicity = RNGfillable_general(grid, vec[0].Grid()); // src has finer or same grid
    int Nsimd        = grid->Nsimd();  // guaranteed to be the same for src.Grid() too
    int osites       = grid->oSites();  // guaranteed to be <= src.Grid()->oSites() by a factor multiplicity
    int words        = sizeof(scalar_object) / sizeof(scalar_type);

    const Integer tW(par().tW);

    autoView(t_v  , t, CpuRead);

    for (auto& src: vec) {

        autoView(src_v, src, CpuWrite);

        thread_for( ss, osites, {

            ExtractBuffer<scalar_object> buf(Nsimd);
            ExtractBuffer<Integer> tbuf(Nsimd);

            for (int m = 0; m < multiplicity; m++) {  // Draw from same generator multiplicity times

                int sm = multiplicity * ss + m;  // Maps the generator site to the fine site

                extract(t_v[sm],tbuf);

                for (int si = 0; si < Nsimd; si++) {

                    scalar_type *pointer = (scalar_type *)&buf[si];

                    if (tbuf[si] == tW) {
                        int gdx = rng.generator_idx(ss, si);  // index of generator state
                        rng._gaussian[gdx].reset();
                        for (int idx = 0; idx < words; idx++) {

                            fillScalar(pointer[idx], rng._gaussian[gdx], rng._generators[gdx]);

                            // Normalize complex number
                            Complex c = pointer[idx];
                            pointer[idx] = c/sqrt(c*adj(c));
                        }
                    } else {
                        *pointer = 0.;
                    }
                }
                merge(src_v[sm], buf);
            }
        });
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_RandomWall_hpp_
