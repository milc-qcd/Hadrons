#ifndef Hadrons_MIO_SaveEigenPack_hpp_
#define Hadrons_MIO_SaveEigenPack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Implicitly Restarted Lanczos module                     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class SaveEigenPackPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(SaveEigenPackPar,
                                    std::string,  eigenPack,
                                    std::string,  filename,
                                    bool,         multiFile);
};

template <typename Pack>
class TSaveEigenPack: public Module<SaveEigenPackPar>
{
public:
    typedef typename Pack::Field      Field;
    typedef BaseEigenPack<Field>      BasePack;
public:
    // constructor
    TSaveEigenPack(const std::string name);
    // destructor
    virtual ~TSaveEigenPack(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(SaveFermionEigenPack, TSaveEigenPack<BaseFermionEigenPack<FIMPL>>, MIO);
MODULE_REGISTER_TMP(StagSaveFermionEigenPack, TSaveEigenPack<BaseFermionEigenPack<STAGIMPL>>, MIO);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(SaveFermionEigenPackF, TSaveEigenPack<BaseFermionEigenPack<FIMPLF>>, MIO);
MODULE_REGISTER_TMP(StagSaveFermionEigenPackF, TSaveEigenPack<BaseFermionEigenPack<STAGIMPLF>>, MIO);
#endif

/******************************************************************************
 *                 TSaveEigenPack implementation                 *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack>
TSaveEigenPack<Pack>::TSaveEigenPack(const std::string name)
: Module<SaveEigenPackPar>(name)
{}

// input ///////////////////////////////////////////////////////
template <typename Pack>
std::vector<std::string> TSaveEigenPack<Pack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};

    return in;
}

// output ///////////////////////////////////////////////////////
template <typename Pack>
std::vector<std::string> TSaveEigenPack<Pack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack>
DependencyMap TSaveEigenPack<Pack>::getObjectDependencies(void)
{
    DependencyMap dep;
    
    dep.insert({par().eigenPack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack>
void TSaveEigenPack<Pack>::setup(void)
{
  LOG(Message) << "Setting up eigenpack output to file " << par().filename << std::endl;
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack>
void TSaveEigenPack<Pack>::execute(void)
{
  auto &epack = envGet(Pack, par().eigenPack);
  GridBase *gridIo = epack.evec[0].Grid();

  if (par().filename.empty())
    {
      LOG(Message) << "Must specify an output file" << std::endl;
      assert(!par().filename.empty());
    }
  
  EigenPackIo::writePack<Field, Field>(par().filename,
				       epack.evec, epack.eval, epack.record, 
				       epack.evec.size(), par().multiFile, gridIo);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_SaveEigenPack_hpp_
