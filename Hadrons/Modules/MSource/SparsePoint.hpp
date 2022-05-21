#ifndef Hadrons_MSource_SparsePoint_hpp_
#define Hadrons_MSource_SparsePoint_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         SparsePoint                                 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class SparsePointPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SparsePointPar,
                                    unsigned int, nsparse);
};

template <typename FImpl>
class TSparsePoint: public Module<SparsePointPar>
{
public:
  FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TSparsePoint(const std::string name);
    // destructor
    virtual ~TSparsePoint(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    std::vector<int> nFields_;
};

MODULE_REGISTER_TMP(StagSparsePoint, TSparsePoint<STAGIMPL>, MSource);

/******************************************************************************
 *                 TSparsePoint implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSparsePoint<FImpl>::TSparsePoint(const std::string name)
: Module<SparsePointPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSparsePoint<FImpl>::getInput(void)
{
    std::vector<std::string> in;
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSparsePoint<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName(), getName()+"_shift"};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSparsePoint<FImpl>::setup(void)
{
    int nSparsePoint = par().nsparse;

    std::div_t divs;
    nFields_.resize(Nd,0);
    int total = 1;
    for (int d = 0 ; d < Nd ; d++) {
        divs = std::div(env().getDim(d),par().nsparse);
        nFields_[d] = std::max(1,divs.quot);
        total *= nFields_[d];
    }
    envCreate(std::vector<PropagatorField>, getName(), 1, total,
              envGetGrid(PropagatorField));
    envCreate(std::vector<Integer>, getName()+"_shift", 1, total, 0);
    
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSparsePoint<FImpl>::execute(void)
{
    typedef typename PropagatorField::scalar_object siteField;
    LatticeInteger coor(envGetGrid(PropagatorField));
    int            nSparse   = par().nsparse;
    auto           &fields   = envGet(std::vector<PropagatorField>,getName());
    auto           &time_shift   = envGet(std::vector<Integer>,getName()+"_shift");
    int i = 0;
    Coordinate gcoor;
    bool desiredPoint;

    auto grid = fields[0].Grid();


    siteField pointSource = siteField(1.0);

    // Iterator over every site on the lattice
    for (int g=0; g < grid->_gsites; g++){

        grid->GlobalIndexToGlobalCoor(g,gcoor);

        desiredPoint = true;
        for (auto component:gcoor) {
            if (mod(component,par().nsparse) != 0) {
                desiredPoint = false;
                break;
            }
        }

        // If all the site's coordinates are multiples of nsparse create point source
        if (desiredPoint) {
            fields[i] = 0.;
            pokeSite(pointSource,fields[i],gcoor);
            time_shift[i] = gcoor[Tp];

            i++;
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_SparsePoint_hpp_
