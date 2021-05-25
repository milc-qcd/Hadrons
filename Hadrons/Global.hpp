/*
 * Global.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: ferben <ferben@debian.felix.com>
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

#ifndef Hadrons_Global_hpp_
#define Hadrons_Global_hpp_

#include <atomic>
#include <set>
#include <stack>
#include <thread>
#include <regex>
#include <Grid/Grid.h>
#include <cxxabi.h>

#ifndef SITE_SIZE_TYPE
#define SITE_SIZE_TYPE size_t
#endif

#ifndef DEFAULT_ASCII_PREC
#define DEFAULT_ASCII_PREC 16
#endif

#ifndef HADRONS_XML_TOPLEV
#define HADRONS_XML_TOPLEV "grid"
#endif 

#define ARG(...) __VA_ARGS__

/* the 'using Grid::operator<<;' statement prevents a very nasty compilation
 * error with GCC 5 (clang & GCC 6 compile fine without it).
 */

#define BEGIN_HADRONS_NAMESPACE \
namespace Grid {\
namespace Hadrons {\
using Grid::operator<<;\
using Grid::operator>>;
#define END_HADRONS_NAMESPACE }}

#define BEGIN_MODULE_NAMESPACE(name)\
namespace name {\
using Grid::operator<<;\
using Grid::operator>>;

#define END_MODULE_NAMESPACE }

#define _HADRONS_IMPL(impl, sub) impl##sub
#define HADRONS_IMPL(impl, sub)   _HADRONS_IMPL(impl, sub)

#ifndef FIMPLBASE
#define FIMPLBASE WilsonImpl
#endif
#define FIMPL  HADRONS_IMPL(FIMPLBASE, R)
#define FIMPLF HADRONS_IMPL(FIMPLBASE, F)
#define FIMPLD HADRONS_IMPL(FIMPLBASE, D)

#ifndef ZFIMPLBASE
#define ZFIMPLBASE ZWilsonImpl
#endif
#define ZFIMPL  HADRONS_IMPL(ZFIMPLBASE, R)
#define ZFIMPLF HADRONS_IMPL(ZFIMPLBASE, F)
#define ZFIMPLD HADRONS_IMPL(ZFIMPLBASE, D)

#ifndef STAGIMPLBASE
#define STAGIMPLBASE StaggeredImpl // use 4d for impl for now since 5d is Lsvectorised
#endif
#define STAGIMPL  HADRONS_IMPL(STAGIMPLBASE, R)
#define STAGIMPLF HADRONS_IMPL(STAGIMPLBASE, F)
#define STAGIMPLD HADRONS_IMPL(STAGIMPLBASE, D)

#ifndef SIMPLBASE
#define SIMPLBASE ScalarImplC
#endif
#define SIMPL  HADRONS_IMPL(SIMPLBASE, R)
#define SIMPLF HADRONS_IMPL(SIMPLBASE, F)
#define SIMPLD HADRONS_IMPL(SIMPLBASE, D)

#ifndef GIMPLBASE
#define GIMPLBASE PeriodicGimpl
#endif
#define GIMPL  HADRONS_IMPL(GIMPLBASE, R)
#define GIMPLF HADRONS_IMPL(GIMPLBASE, F)
#define GIMPLD HADRONS_IMPL(GIMPLBASE, D)

#ifndef STAGIMPLBASE
#define STAGIMPLBASE StaggeredImpl // use 4d for impl for now since 5d is Lsvectorised
#endif
#define STAGIMPL  HADRONS_IMPL(STAGIMPLBASE, R)
#define STAGIMPLF HADRONS_IMPL(STAGIMPLBASE, F)
#define STAGIMPLD HADRONS_IMPL(STAGIMPLBASE, D)

BEGIN_HADRONS_NAMESPACE

// type aliases
#define BASIC_TYPE_ALIASES(Impl, suffix)\
typedef typename Impl::Field                         ScalarField##suffix;\
typedef typename Impl::PropagatorField               PropagatorField##suffix;\
typedef typename Impl::SitePropagator::scalar_object SitePropagator##suffix;\
typedef typename Impl::ComplexField                  ComplexField##suffix;\
typedef std::vector<SitePropagator##suffix>          SlicedPropagator##suffix;\
typedef std::vector<typename ComplexField##suffix::vector_object::scalar_object> SlicedComplex##suffix;

#define FERM_TYPE_ALIASES(FImpl, suffix)\
BASIC_TYPE_ALIASES(FImpl, suffix);\
typedef FermionOperator<FImpl>                     FMat##suffix;\
typedef typename FImpl::FermionField               FermionField##suffix;\
typedef typename FImpl::GaugeField                 GaugeField##suffix;\
typedef typename FImpl::DoubledGaugeField          DoubledGaugeField##suffix;\
typedef Lattice<iSpinMatrix<typename FImpl::Simd>> SpinMatrixField##suffix;\
typedef Lattice<iColourVector<typename FImpl::Simd>> ColourVectorField##suffix;

#define GAUGE_TYPE_ALIASES(GImpl, suffix)\
typedef typename GImpl::GaugeField GaugeField##suffix;

#define SOLVER_TYPE_ALIASES(FImpl, suffix)\
typedef Solver<FImpl> Solver##suffix;

#define SINK_TYPE_ALIASES(suffix)\
typedef std::function<SlicedPropagator##suffix\
                      (const PropagatorField##suffix &)> SinkFn##suffix;

// logger
class HadronsLogger: public Logger
{
public:
    HadronsLogger(int on, std::string nm): Logger("Hadrons", on, nm,
                                                  GridLogColours, "BLACK"){};
};

#define LOG(channel) std::cout << HadronsLog##channel
#define HADRONS_DEBUG_VAR(var) LOG(Debug) << #var << "= " << (var) << std::endl;

extern HadronsLogger HadronsLogError;
extern HadronsLogger HadronsLogWarning;
extern HadronsLogger HadronsLogMessage;
extern HadronsLogger HadronsLogIterative;
extern HadronsLogger HadronsLogDebug;
extern HadronsLogger HadronsLogIRL;

void initLogger(void);

// singleton pattern
#define SINGLETON(name)\
public:\
    name(const name &e) = delete;\
    void operator=(const name &e) = delete;\
    static name & getInstance(void)\
    {\
        static name e;\
        return e;\
    }\
private:\
    name(void);

#define SINGLETON_DEFCTOR(name)\
public:\
    name(const name &e) = delete;\
    void operator=(const name &e) = delete;\
    static name & getInstance(void)\
    {\
        static name e;\
        return e;\
    }\
private:\
    name(void) = default;

// type utilities
template <typename T>
const std::type_info * typeIdPt(const T &x)
{
    return &typeid(x);
}

template <typename T>
const std::type_info * typeIdPt(void)
{
    return &typeid(T);
}

size_t typeHash(const std::type_info *info);

template <typename T>
size_t typeHash(const T &x)
{
    return typeHash(typeIdPt(x));
}

template <typename T>
size_t typeHash(void)
{
    return typeHash(typeIdPt<T>());
}

template <typename T, typename U>
bool sameType(const T &x, const U &y)
{
    return (typeHash(x) == typeHash(y));
}

template <typename T, typename U>
bool sameType(void)
{
    return (typeHash<T>() == typeHash<U>());
}

std::string typeName(const std::type_info *info);

template <typename T>
std::string typeName(const T &x)
{
    return typeName(typeIdPt(x));
}

template <typename T>
std::string typeName(void)
{
    return typeName(typeIdPt<T>());
}

// default writers/readers
extern const std::string resultFileExt;

#ifdef HAVE_HDF5
typedef Hdf5Reader ResultReader;
typedef Hdf5Writer ResultWriter;
#else
typedef XmlReader ResultReader;
typedef XmlWriter ResultWriter;
#endif

// recursive mkdir
#define MAX_PATH_LENGTH 512u
int         mkdir(const std::string dirName);
std::string basename(const std::string &s);
std::string dirname(const std::string &s);
void        makeFileDir(const std::string filename, GridBase *g = nullptr);

#define HADRONS_TYPEDEF_TEMPLATE_BRANCH(type_name,TImpl,condition,true_type,false_type)\
template<typename T, typename... Args>\
using type_name = typename std::conditional<condition, true_type<Args...>, false_type<Args...> >::type;

#define HADRONS_TYPEDEF_BRANCH_STAGGERED(type_name,FImpl,staggered_type,non_staggered_type)\
HADRONS_TYPEDEF_TEMPLATE_BRANCH(type_name,FImpl,HADRONS_IS_STAGGERED_IMPLEMENTATION(T), staggered_type, non_staggered_type)

#define HADRONS_IS_STAGGERED_IMPLEMENTATION(FImpl)\
(std::is_same<FImpl,STAGIMPLD>::value || std::is_same<FImpl,STAGIMPLF>::value)

// default Schur convention
#ifndef HADRONS_DEFAULT_SCHUR
#define HADRONS_DEFAULT_SCHUR DiagTwo
#endif

#ifndef HADRONS_DEFAULT_SCHUR_STAGGERED
#define HADRONS_DEFAULT_SCHUR_STAGGERED Staggered
#endif

#define _HADRONS_SCHUR_OP_(conv) Schur##conv##Operator
#define HADRONS_SCHUR_OP(conv) _HADRONS_SCHUR_OP_(conv)
#define HADRONS_DEFAULT_SCHUR_OP HADRONS_SCHUR_OP(HADRONS_DEFAULT_SCHUR)
#define _HADRONS_SCHUR_SOLVE_(conv) SchurRedBlack##conv##Solve
#define HADRONS_SCHUR_SOLVE(conv) _HADRONS_SCHUR_SOLVE_(conv)
#define HADRONS_DEFAULT_SCHUR_SOLVE HADRONS_SCHUR_SOLVE(HADRONS_DEFAULT_SCHUR)
// #define _HADRONS_SCHUR_A2A_(conv) A2AVectorsSchur##conv
// #define HADRONS_SCHUR_A2A(conv) _HADRONS_SCHUR_A2A_(conv)
// #define HADRONS_DEFAULT_SCHUR_A2A HADRONS_SCHUR_A2A(HADRONS_DEFAULT_SCHUR)
// #define HADRONS_DEFAULT_SCHUR_A2A_STAGGERED HADRONS_SCHUR_A2A(HADRONS_DEFAULT_SCHUR_STAGGERED)
#define HADRONS_DEFAULT_SCHUR_OP_STAGGERED HADRONS_SCHUR_OP(HADRONS_DEFAULT_SCHUR_STAGGERED)
#define HADRONS_DEFAULT_SCHUR_SOLVE_STAGGERED HADRONS_SCHUR_SOLVE(HADRONS_DEFAULT_SCHUR_STAGGERED)

#define HADRONS_DEFINE_SCHUR_OP(name,FImpl)\
HADRONS_TYPEDEF_BRANCH_STAGGERED(name##_macro,FImpl,HADRONS_DEFAULT_SCHUR_OP_STAGGERED,HADRONS_DEFAULT_SCHUR_OP)\
template <typename... Args>\
using name = name##_macro<FImpl,Args...>;

#define HADRONS_DEFINE_SCHUR_SOLVE(name,FImpl)\
HADRONS_TYPEDEF_BRANCH_STAGGERED(name##_macro,FImpl,HADRONS_DEFAULT_SCHUR_SOLVE_STAGGERED,HADRONS_DEFAULT_SCHUR_SOLVE)\
template <typename... Args>\
using name = name##_macro<FImpl,Args...>;

// stringify macro
#define _HADRONS_STR(x) #x
#define HADRONS_STR(x) _HADRONS_STR(x)

// pretty print time profile
void printTimeProfile(const std::map<std::string, GridTime> &timing, GridTime total);

// token replacement utility
template <typename T>
void tokenReplace(std::string &str, const std::string token,
                  const T &x, const std::string mark = "@")
{
    std::string fullToken = mark + token + mark;
    
    auto pos = str.find(fullToken);
    if (pos != std::string::npos)
    {
        str.replace(pos, fullToken.size(), std::to_string(x));
    }
}

// generic correlator class
template <typename Metadata, typename Scalar = Complex>
struct Correlator: Serializable
{
    GRID_SERIALIZABLE_CLASS_MEMBERS(ARG(Correlator<Metadata, Scalar>),
                                    Metadata,             info,
                                    std::vector<Scalar>, corr);
};

// check if grid is initlialised
bool isGridInit(void);

// thread random wait
void randomWait(const unsigned int maxMs, GridBase *grid = nullptr);

END_HADRONS_NAMESPACE

#include <Hadrons/Exceptions.hpp>

#endif // Hadrons_Global_hpp_
