#ifndef Hadrons_MSource_SparseSpinColorDiagonal_hpp_
#define Hadrons_MSource_SparseSpinColorDiagonal_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         SparseSpinColorDiagonal                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class SparseSpinColorDiagonalPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SparseSpinColorDiagonalPar,
                                    unsigned int, nsparse);
};

template <typename FImpl>
class TSparseSpinColorDiagonal: public Module<SparseSpinColorDiagonalPar>
{
public:
  FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSparseSpinColorDiagonal(const std::string name);
    // destructor
    virtual ~TSparseSpinColorDiagonal(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SparseSpinColorDiagonal, TSparseSpinColorDiagonal<STAGIMPL>, MSource);

/******************************************************************************
 *                 TSparseSpinColorDiagonal implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSparseSpinColorDiagonal<FImpl>::TSparseSpinColorDiagonal(const std::string name)
: Module<SparseSpinColorDiagonalPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSparseSpinColorDiagonal<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSparseSpinColorDiagonal<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName()+"_shift"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSparseSpinColorDiagonal<FImpl>::setup(void)
{
    int nSparseSpinColorDiagonal = par().nsparse;

    
    envCreate(std::vector<PropagatorField>, getName(), 1, pow(par().nsparse,Nd),
              envGetGrid(PropagatorField));
    envCreate(std::vector<Integer>, getName()+"_shift", 1, pow(par().nsparse,Nd), 0);
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSparseSpinColorDiagonal<FImpl>::execute(void)
{
    std::div_t divs;
    LatticeInteger coor(envGetGrid(PropagatorField));

    int nSparse = par().nsparse;
    auto &fields = envGet(std::vector<PropagatorField>,getName());
    auto &shifts = envGet(std::vector<Integer>,getName()+"_shift");

    // Create the sparse source pattern starting at the origin
    fields[0] = 1.;
    for(int d = 0; d < Nd; ++d) 
    {
        LatticeCoordinate(coor, d);
        fields[0] = where(mod(coor,nSparse),0.*fields[0],fields[0]);
    }

    auto norm = norm2(fields[0]);
    fields[0] = sqrt(Nc/norm)*fields[0];

    for (int i = 1; i < fields.size(); i++) {
        fields[i] = fields[0];
        for (int d = 0; d < Nd; ++d) {
            divs = std::div(i, pow(nSparse, Nd-(d+1)));
            if (divs.quot != 0) {
                fields[i] = Cshift(fields[i], d, divs.quot);

                // If we're shifting in the time direction
                if (d == Tp) {
                    shifts[i] = (Integer)divs.quot;
                }
            }
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_SparseSpinColorDiagonal_hpp_
