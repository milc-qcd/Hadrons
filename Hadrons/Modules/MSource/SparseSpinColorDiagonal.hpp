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

MODULE_REGISTER_TMP(SparseColorDiagonal, TSparseSpinColorDiagonal<STAGIMPL>, MSource);

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

    
    envCreate(std::vector<PropagatorField>, getName(), 1, pow(par().nsparse,Nd-1)*env().getDim(Tp),
              envGetGrid(PropagatorField));
    envCreate(std::vector<Integer>, getName()+"_shift", 1, pow(par().nsparse,Nd-1)*env().getDim(Tp), 0);
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSparseSpinColorDiagonal<FImpl>::execute(void)
{
    LatticeInteger coor(envGetGrid(PropagatorField));
    std::div_t     divs;
    int            nSparse   = par().nsparse;
    auto           &fields   = envGet(std::vector<PropagatorField>,getName());
    auto           &shifts   = envGet(std::vector<Integer>,getName()+"_shift");

    // Create the sparse source pattern starting at the origin
    fields[0] = 1.;
    for(int d = 0; d < Nd; ++d) 
    {
        if (d != Tp) {
            LatticeCoordinate(coor, d);
            fields[0] = where(mod(coor,nSparse),0.*fields[0],fields[0]);
       }
    }

    // Normalize field
    Real norm = norm2(fields[0])/(FImpl::Dimension*env().getDim(Tp));
    norm = 1/sqrt(norm);

    fields[0] = ComplexD(norm,0.)*fields[0];

    int repeatSet = 1;
    // Shift in each direction
    for (int d = 0; d < Nd; ++d) {
        if (d == Tp) 
            continue;

        // Loop over all previous fields
        for (int i = 0; i < repeatSet; i++) {
            int nShifts = nSparse-1;
            // Iteratively Shift and save fields in direction d
            for (int n = 0; n < nShifts; n++) {
                if (n==0)
                    fields[repeatSet + i*nShifts + n] = Cshift(fields[i], d, 1);
                else
                    fields[repeatSet + i*nShifts + n] = Cshift(fields[repeatSet + i*nShifts + n-1], d, 1);
            }
        }
        repeatSet = repeatSet*nSparse;
    }

    LatticeCoordinate(coor, Tp);

    // Dilute in time direction
    // Loop backwards through fields so we're not overwriting the source fields generated above
    Integer t = env().getDim(Tp)-1;
    int j = repeatSet-1;
    for (int k = fields.size()-1; k >= 0; k--) {
        if (j < 0) {
            j = repeatSet-1;
            t--;
        }

        fields[k] = where(coor==t, fields[j], 0.*fields[j]);
        shifts[k] = t;
        j--;

    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_SparseSpinColorDiagonal_hpp_
