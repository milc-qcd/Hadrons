/*
 * Environment.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: ferben <ferben@debian.felix.com>
 * Author: nelsonlachini <nelsonlachini@gmail.com>
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

#ifndef Hadrons_Environment_hpp_
#define Hadrons_Environment_hpp_

#include <Hadrons/Global.hpp>

BEGIN_HADRONS_NAMESPACE

template <typename T>
struct IsLattice : std::false_type {};
template <typename Vobj>
struct IsLattice<Lattice<Vobj>> : std::true_type {};
// std::vector<Lattice<Vobj>> — e.g. time-diluted noise sources, multi-gamma
// propagator outputs. These are scatter-leaves exactly like scalar lattices:
// each element is Grid_split onto the subgrid independently, and the producing
// module is NOT rebuilt (the data, not the computation, is what's needed).
template <typename T>
struct IsLatticeVector : std::false_type {};
template <typename Vobj>
struct IsLatticeVector<std::vector<Lattice<Vobj>>> : std::true_type {};
// anything scatterable onto a subgrid (scalar lattice OR vector of lattices)
template <typename T>
struct IsScatterable
: std::integral_constant<bool, IsLattice<T>::value || IsLatticeVector<T>::value> {};

// Grid-independent metadata that can be cheaply deep-copied into the shadow
// store instead of triggering a full producer rebuild. E.g. std::vector<Integer>
// (time-dilution source shifts), scalar parameters. These contain no Grid
// lattices and no RNG state, so a copy is valid on any subgrid.
template <typename T>
struct IsCloneableMetadata : std::false_type {};
template <typename T>
struct IsCloneableMetadata<std::vector<T>> : std::is_arithmetic<T> {};

/******************************************************************************
 *                         Global environment                                 *
 ******************************************************************************/
class Object
{
public:
    Object(void) = default;
    virtual ~Object(void) = default;
    // type-erased scatter into the active per-VType subgrid; returns nullptr
    // if the held object is not a scatterable Lattice field
    virtual std::unique_ptr<Object> splitTo(void) const
    {
        return nullptr;
    }
    // lightweight (non-mutating) lattice test — true iff the held object is a
    // scatterable Lattice field. Used by the GC scheduler to tell scatter-leaf
    // objects (lattices) from rebuildable objects (solvers/actions) without
    // performing a Grid_split.
    virtual bool isLattice(void) const
    {
        return false;
    }
    // type-erased deep copy for grid-independent metadata (e.g.
    // std::vector<Integer>). Returns nullptr for types that are NOT safe to
    // clone (solvers, actions, sink functions — these are grid-bound and must
    // be rebuilt on the subgrid via setup()+execute()).
    virtual std::unique_ptr<Object> clone(void) const
    {
        return nullptr;
    }
};

template <typename T>
class Holder: public Object
{
public:
    Holder(void) = default;
    Holder(T *pt);
    virtual ~Holder(void) = default;
    T &       get(void) const;
    T *       getPt(void) const;
    void      reset(T *pt);
    std::unique_ptr<Object> splitTo(void) const override;
    bool      isLattice(void) const override
    {
        return IsScatterable<T>::value;
    }
    std::unique_ptr<Object> clone(void) const override
    {
        if constexpr (IsCloneableMetadata<T>::value)
        {
            return std::make_unique<Holder<T>>(new T(*objPt_));
        }
        else
        {
            return nullptr;
        }
    }
private:
    std::unique_ptr<T> objPt_{nullptr};
};

// Handle carrying the split sub-communicator grid built by the VirtualMachine
// and consumed by the scope-aware grid/object fetch. The owned GridCartesian
// frees its MPI sub-communicator on destruction.
struct SubGrids
{
    std::unique_ptr<GridCartesian> grid;      // base 4d sub-communicator grid (owned)
    Coordinate                     mpiSplit;   // per-subcomm processor layout
    int                            nrhs;       // number of sub-comms
    int                            me;         // this rank's sub-comm index [0, nrhs)
};

#define DEFINE_ENV_ALIAS \
inline Environment & env(void) const\
{\
    return Environment::getInstance();\
}

#define DEFINE_ENV_LAMBDA \
auto env = [](void)->Environment &{return Environment::getInstance();}

class Environment
{
    SINGLETON(Environment);
public:
    typedef SITE_SIZE_TYPE                         Size;
    typedef std::unique_ptr<GridCartesian>         GridPt;
    typedef std::unique_ptr<GridRedBlackCartesian> GridRbPt;
    typedef std::unique_ptr<GridParallelRNG>       RngPt;
    typedef std::unique_ptr<GridSerialRNG>         SerialRngPt;
    GRID_SERIALIZABLE_ENUM(Storage, undef, standard, 0, cache, 1, temporary, 2);
private:
    struct ObjInfo
    {
        Size                      size{0};
        Storage                   storage{Storage::standard};
        unsigned int              Ls{0};
        const std::type_info      *type{nullptr}, *derivedType{nullptr};
        std::string               name;
        int                       module{-1};
        std::unique_ptr<Object>   data{nullptr};
        std::set<unsigned int>    dependency;
    };
    typedef std::pair<size_t, unsigned int>     FineGridKey;
    typedef std::pair<size_t, std::vector<int>> CoarseGridKey;
public:
    // grids
    Coordinate simdDecomposition(const unsigned int nd, const unsigned int nSimd);
    template <typename VType = vComplex>
    void                    createGrid(const unsigned int Ls);
    template <typename VType = vComplex>
    void                    createCoarseGrid(const std::vector<int> &blockSize,
                                             const unsigned int Ls);
    template <typename VType = vComplex>
    void                    createSliceGrid(const unsigned int orthDim);
    template <typename VType = vComplex>
    void                    createSubGrid(void);
    template <typename VType = vComplex>
    GridCartesian *         getGrid(void);
    template <typename VType = vComplex>
    GridRedBlackCartesian * getRbGrid(void);
    template <typename VType = vComplex>
    GridCartesian *         getCoarseGrid(const std::vector<int> &blockSize);
    template <typename VType = vComplex>
    GridCartesian *         getSliceGrid(const unsigned int orthDir);
    template <typename VType = vComplex>
    GridCartesian *         getGrid(const unsigned int Ls);
    template <typename VType = vComplex>
    GridRedBlackCartesian * getRbGrid(const unsigned int Ls);
    template <typename VType = vComplex>
    GridCartesian *         getCoarseGrid(const std::vector<int> &blockSize,
                                          const unsigned int Ls);
    std::vector<int>        getDim(void) const;
    int                     getDim(const unsigned int mu) const;
    unsigned int            getNd(void) const;
    double                  getVolume(void) const;
    // random number generator
    GridParallelRNG *       get4dRng(void);
    GridSerialRNG *         getSerialRng(void);
    // subgrid scope management
    void                    setActiveSubGrid(GridCartesian *subGrid,
                                             const int splitKey);
    void                    clearActiveSubGrid(void);
    bool                    isSubGridActive(void) const;
    int                     getActiveSplitKey(void) const { return activeSplitKey_; }
    // toggle the subgrid scope on/off for the current module without clearing
    // the subgrid caches or shadow store (which must persist across the split
    // phase). setActiveSubGrid/clearActiveSubGrid bracket the whole split phase.
    void                    setSubGridScope(const bool on);
    // subgrid shadow store (scattered/rebuilt global objects, populated by VM)
    void                    addShadowObject(const unsigned int address,
                                            const int splitKey,
                                            std::unique_ptr<Object> obj);
    bool                    hasShadowObject(const unsigned int address) const;
    // scatter a lattice global object onto the active subgrid; returns true if
    // the object was scattered (it is a Lattice), false otherwise (non-lattice
    // objects must be rebuilt by the VM via setup() in shadow-create mode).
    bool                    scatterObject(const unsigned int address,
                                          const int splitKey);
    // clone a grid-independent non-lattice metadata object (e.g.
    // std::vector<Integer>) into the shadow store. Returns true if the object
    // was cloned, false if it is not cloneable (the VM then rebuilds the
    // producer via setup()+execute() in shadow-create mode).
    bool                    cloneObject(const unsigned int address,
                                        const int splitKey);
    // lightweight, non-mutating test: true iff the (created) object is a
    // scatterable Lattice field. Used by the GC scheduler to mirror
    // ensureShadowed's lattice-leaf behaviour when computing which global
    // objects must survive for the subgrid rebuild.
    bool                    isLatticeObject(const unsigned int address) const;
    // route subsequent createDerivedObject calls into the shadow store so a
    // non-lattice object can be rebuilt on the subgrid leaving the world-grid
    // global copy intact for global consumers.
    void                    enterShadowCreate(const int splitKey);
    void                    exitShadowCreate(void);
    // general memory management
    void                    addObject(const std::string name,
                                      const int moduleAddress = -1);
    template <typename B, typename T, typename ... Ts>
    void                    createDerivedObject(const std::string name,
                                                const Environment::Storage storage,
                                                const unsigned int Ls,
                                                Ts && ... args);
    template <typename T, typename ... Ts>
    void                    createObject(const std::string name,
                                         const Environment::Storage storage,
                                         const unsigned int Ls,
                                         Ts && ... args);
    void                    setObjectStorage(const unsigned int objAddress,
                                             const Environment::Storage storage);
    void                    setObjectModule(const unsigned int objAddress,
                                            const int modAddress);
    void                    addObjectDependency(const unsigned int objAddress,
                                                const unsigned int depAddress);
    void                    removeObjectDependency(const unsigned int objAddress,
                                                   const unsigned int depAddress);
    const std::set<unsigned int> &getObjectDependencies(const unsigned int objAddress) const;
    bool                    hasDependency(const unsigned int objAddress,
                                          const unsigned int depAddress) const;
    template <typename B, typename T>
    T *                     getDerivedObject(const unsigned int address) const;
    template <typename B, typename T>
    T *                     getDerivedObject(const std::string name) const;
    template <typename T>
    T *                     getObject(const unsigned int address) const;
    template <typename T>
    T *                     getObject(const std::string name) const;
    unsigned int            getMaxAddress(void) const;
    unsigned int            getObjectAddress(const std::string name) const;
    std::string             getObjectName(const unsigned int address) const;
    std::string             getObjectType(const unsigned int address) const;
    std::string             getObjectType(const std::string name) const;
    std::string             getObjectDerivedType(const unsigned int address) const;
    std::string             getObjectDerivedType(const std::string name) const;
    Size                    getObjectSize(const unsigned int address) const;
    Size                    getObjectSize(const std::string name) const;
    Storage                 getObjectStorage(const unsigned int address) const;
    Storage                 getObjectStorage(const std::string name) const;
    int                     getObjectModule(const unsigned int address) const;
    int                     getObjectModule(const std::string name) const;
    unsigned int            getObjectLs(const unsigned int address) const;
    unsigned int            getObjectLs(const std::string name) const;
    bool                    hasObject(const unsigned int address) const;
    bool                    hasObject(const std::string name) const;
    bool                    hasCreatedObject(const unsigned int address) const;
    bool                    hasCreatedObject(const std::string name) const;
    bool                    isObject5d(const unsigned int address) const;
    bool                    isObject5d(const std::string name) const;
    template <typename T>
    bool                    isObjectOfType(const unsigned int address) const;
    template <typename T>
    bool                    isObjectOfType(const std::string name) const;
    template <typename B, typename T>
    bool                    isObjectOfDerivedType(const unsigned int address) const;
    template <typename B, typename T>
    bool                    isObjectOfDerivedType(const std::string name) const;
    Environment::Size       getTotalSize(void) const;
    void                    freeObject(const unsigned int address, const bool recursive = false);
    void                    freeObject(const std::string name, const bool recursive = false);
    void                    freeSet(const std::set<unsigned int> &objects);
    void                    freeAll(void);
    void                    protectObjects(const bool protect);
    bool                    objectsProtected(void) const;
    // print environment content
    void                    printContent(void) const;
private:
    // general
    double                              vol_;
    bool                                protect_{true}, simdReverse_{false};
    // grids
    std::vector<int>                    dim_;
    std::vector<bool>                   simdMask_;
    std::map<FineGridKey, GridPt>       grid3d_;
    std::map<FineGridKey, GridPt>       grid4d_;
    std::map<FineGridKey, GridPt>       grid5d_;
    std::map<FineGridKey, GridRbPt>     gridRb4d_;
    std::map<FineGridKey, GridRbPt>     gridRb5d_;
    std::map<CoarseGridKey, GridPt>     gridCoarse4d_;
    std::map<CoarseGridKey, GridPt>     gridCoarse5d_;
    unsigned int                        nd_;
    // subgrid scope state
    GridCartesian                          *activeSubGrid_{nullptr};
    // subGridScopeOn_ toggles per-module: it is true only while a split-phase
    // module executes. activeSubGrid_ stays set for the whole split phase (so
    // the subgrid caches / shadow store persist), but getGrid/getObject only
    // serve the subgrid/shadow when subGridScopeOn_ is true — global modules
    // interleaved between split-phase modules thus run on the world grid.
    bool                                    subGridScopeOn_{false};
    // (no activeSubRbGrid_ — createSubGrid builds per-VType RB grids from per-VType 4d subgrids)
    int                                     activeSplitKey_{-1};
    std::map<FineGridKey, GridPt>           gridSub4d_;
    std::map<FineGridKey, GridRbPt>         gridSubRb4d_;
    std::map<std::pair<unsigned int, int>, std::unique_ptr<Object>> shadowStore_;
    // shadow-create mode: when true, createDerivedObject routes the newly built
    // object into shadowStore_ (keyed by shadowCreateKey_) instead of object_.
    // Used by the VM to rebuild non-lattice global objects (e.g. fermion
    // actions) onto the subgrid without overwriting the world-grid original.
    bool                                    shadowCreateMode_{false};
    int                                     shadowCreateKey_{-1};
    // random number generator
    RngPt                               rng4d_{nullptr};
    SerialRngPt                         rngSerial_{nullptr};
    // object store
    std::vector<ObjInfo>                object_;
    std::map<std::string, unsigned int> objectAddress_;
};

/******************************************************************************
 *                       Holder template implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T>
Holder<T>::Holder(T *pt)
: objPt_(pt)
{}

// access //////////////////////////////////////////////////////////////////////
template <typename T>
T & Holder<T>::get(void) const
{
    return *objPt_.get();
}

template <typename T>
T * Holder<T>::getPt(void) const
{
    return objPt_.get();
}

template <typename T>
void Holder<T>::reset(T *pt)
{
    objPt_.reset(pt);
}

// subgrid scatter ////////////////////////////////////////////////////////////
template <typename T>
std::unique_ptr<Object> Holder<T>::splitTo(void) const
{
    if constexpr (IsLattice<T>::value)
    {
        using VType = typename T::vector_type;
        auto  &env  = Environment::getInstance();
        auto  *sg   = env.template getGrid<VType>();
        auto   split = std::make_unique<T>(sg);
        Grid_split(*objPt_, *split);
        return std::make_unique<Holder<T>>(split.release());
    }
    else if constexpr (IsLatticeVector<T>::value)
    {
        // std::vector<Lattice<Vobj>>: scatter each element onto the subgrid.
        // All elements share the same lattice grid, so derive VType from the
        // element type and split element-by-element.
        using ElemT  = typename T::value_type;          // Lattice<Vobj>
        using VType  = typename ElemT::vector_type;
        auto  &env   = Environment::getInstance();
        auto  *sg    = env.template getGrid<VType>();
        auto   split = std::make_unique<T>();
        split->reserve(objPt_->size());
        for (auto &e : *objPt_)
        {
            split->emplace_back(sg);
            Grid_split(e, split->back());
        }
        return std::make_unique<Holder<T>>(split.release());
    }
    else
    {
        return nullptr;
    }
}

/******************************************************************************
 *                     Environment template implementation                    *
 ******************************************************************************/
// grids ///////////////////////////////////////////////////////////////////////
#define HADRONS_DUMP_GRID(...)\
LOG(Debug) << "New grid " << (__VA_ARGS__) << std::endl;\
LOG(Debug) << " - cb  : " << (__VA_ARGS__)->_isCheckerBoarded << std::endl;\
LOG(Debug) << " - fdim: " << (__VA_ARGS__)->_fdimensions << std::endl;\
LOG(Debug) << " - gdim: " << (__VA_ARGS__)->_gdimensions << std::endl;\
LOG(Debug) << " - ldim: " << (__VA_ARGS__)->_ldimensions << std::endl;\
LOG(Debug) << " - rdim: " << (__VA_ARGS__)->_rdimensions << std::endl;\
LOG(Debug) << " - SIMD: " << (__VA_ARGS__)->_simd_layout << std::endl;

template <typename VType>
void Environment::createGrid(const unsigned int Ls)
{
    size_t hash = typeHash<VType>();

    if (grid4d_.find({hash, 1}) == grid4d_.end())
    {
        grid4d_[{hash, 1}].reset(
            SpaceTimeGrid::makeFourDimGrid(getDim(), 
                                        simdDecomposition(getNd(), VType::Nsimd()),
                                        GridDefaultMpi()));
        HADRONS_DUMP_GRID(grid4d_[{hash, 1}].get());
        gridRb4d_[{hash, 1}].reset(
            SpaceTimeGrid::makeFourDimRedBlackGrid(grid4d_[{hash, 1}].get()));
        HADRONS_DUMP_GRID(gridRb4d_[{hash, 1}].get());
    }
    if (grid5d_.find({hash, Ls}) == grid5d_.end())
    {
        auto g = grid4d_[{hash, 1}].get();
        
        grid5d_[{hash, Ls}].reset(SpaceTimeGrid::makeFiveDimGrid(Ls, g));
        HADRONS_DUMP_GRID(grid5d_[{hash, Ls}].get());
        gridRb5d_[{hash, Ls}].reset(SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, g));
        HADRONS_DUMP_GRID(gridRb5d_[{hash, Ls}].get());
    }
}

template <typename VType>
void Environment::createCoarseGrid(const std::vector<int> &blockSize,
                                   const unsigned int Ls)
{
    int              nd      = getNd();
    std::vector<int> fineDim = getDim(), coarseDim(nd);
    unsigned int     cLs;
    auto             key4d = blockSize, key5d = blockSize;
    size_t           hash  = typeHash<VType>();

    createGrid(Ls);
    for (int d = 0; d < coarseDim.size(); d++)
    {
        coarseDim[d] = fineDim[d]/blockSize[d];
        if (coarseDim[d]*blockSize[d] != fineDim[d])
        {
            HADRONS_ERROR(Size, "Fine dimension " + std::to_string(d) 
                         + " (" + std::to_string(fineDim[d]) 
                         + ") not divisible by coarse dimension ("
                         + std::to_string(coarseDim[d]) + ")"); 
        }
    }
    if (blockSize.size() > nd)
    {
        cLs = Ls/blockSize[nd];
        if (cLs*blockSize[nd] != Ls)
        {
            HADRONS_ERROR(Size, "Fine Ls (" + std::to_string(Ls) 
                         + ") not divisible by coarse Ls ("
                         + std::to_string(cLs) + ")");
        }
    }
    else
    {
        cLs = Ls;
    }
    key4d.resize(nd);
    key5d.push_back(Ls);

    CoarseGridKey hkey4d = {hash, key4d}, hkey5d = {hash, key5d};

    if (gridCoarse4d_.find(hkey4d) == gridCoarse4d_.end())
    {
        gridCoarse4d_[hkey4d].reset(
            SpaceTimeGrid::makeFourDimGrid(coarseDim, 
                simdDecomposition(nd, VType::Nsimd()), GridDefaultMpi()));
        HADRONS_DUMP_GRID(gridCoarse4d_[hkey4d].get());
    }
    if (gridCoarse5d_.find(hkey5d) == gridCoarse5d_.end())
    {
        gridCoarse5d_[hkey5d].reset(
            SpaceTimeGrid::makeFiveDimGrid(cLs, gridCoarse4d_[hkey4d].get()));
        HADRONS_DUMP_GRID(gridCoarse5d_[hkey5d].get());
    }
}

template <typename VType>
void Environment::createSliceGrid(const unsigned int orthDim)
{
    size_t hash = typeHash<VType>();

    if (grid3d_.find({hash, orthDim}) == grid3d_.end())
    {
        GridCartesian *g         = getGrid<VType>();
        int           nd         = static_cast<int>(g->_ndimension);
        unsigned int  hd         = 0;
        Coordinate    latt_size  = g->_gdimensions;
        Coordinate    simd3      = simdDecomposition(nd - 1, VType::Nsimd());
        Coordinate    simd;
        Coordinate    mpi        = g->_processors;

        latt_size[orthDim] = 1;
        for (unsigned int d = 0; d < nd; d++)
        {
            if (d == orthDim)
            {
                simd.push_back(1);
            }
            else
            {
                simd.push_back(simd3[hd]);
                hd++;
            }
        }
        mpi[orthDim] = 1;
        grid3d_[{hash, orthDim}].reset( 
            new GridCartesian(latt_size, simd, mpi, *g));
        HADRONS_DUMP_GRID(grid3d_[{hash, orthDim}].get());
    }
}

template <typename VType>
void Environment::createSubGrid(void)
{
    size_t hash = typeHash<VType>();
    FineGridKey key = {hash, 1};

    if (gridSub4d_.find(key) == gridSub4d_.end())
    {
        Coordinate simd = simdDecomposition(activeSubGrid_->_ndimension,
                                            VType::Nsimd());
        gridSub4d_[key].reset(
            new GridCartesian(activeSubGrid_->_gdimensions, simd,
                              activeSubGrid_->_processors, *activeSubGrid_));
        HADRONS_DUMP_GRID(gridSub4d_[key].get());
        gridSubRb4d_[key].reset(
            SpaceTimeGrid::makeFourDimRedBlackGrid(gridSub4d_[key].get()));
        HADRONS_DUMP_GRID(gridSubRb4d_[key].get());
    }
}

#undef HADRONS_DUMP_GRID

template <typename VType>
GridCartesian * Environment::getGrid(void)
{
    if (activeSubGrid_ && subGridScopeOn_)
    {
        FineGridKey key = {typeHash<VType>(), 1};
        auto it = gridSub4d_.find(key);
        if (it != gridSub4d_.end())
            return it->second.get();
        createSubGrid<VType>();
        return gridSub4d_.at(key).get();
    }
    FineGridKey key = {typeHash<VType>(), 1};

    auto it = grid4d_.find(key);

    if (it != grid4d_.end())
    {
        return it->second.get();
    }
    else
    {
        createGrid<VType>(1);

        return grid4d_.at(key).get();
    }
}

template <typename VType>
GridRedBlackCartesian * Environment::getRbGrid(void)
{
    if (activeSubGrid_ && subGridScopeOn_)
    {
        FineGridKey key = {typeHash<VType>(), 1};
        auto it = gridSubRb4d_.find(key);
        if (it != gridSubRb4d_.end())
            return it->second.get();
        createSubGrid<VType>();
        return gridSubRb4d_.at(key).get();
    }
    FineGridKey key = {typeHash<VType>(), 1};
    auto        it  = gridRb4d_.find(key);

    if (it != gridRb4d_.end())
    {
        return it->second.get();
    }
    else
    {
        createGrid<VType>(1);

        return gridRb4d_.at(key).get();
    }
}

template <typename VType>
GridCartesian * Environment::getCoarseGrid(const std::vector<int> &blockSize)
{
    std::vector<int> s = blockSize;

    s.resize(getNd());

    CoarseGridKey key = {typeHash<VType>(), s};
    auto          it  = gridCoarse4d_.find(key);

    if (it != gridCoarse4d_.end())
    {
        return it->second.get();
    }
    else
    {
        createCoarseGrid<VType>(blockSize, 1);
        
        return gridCoarse4d_.at(key).get();
    }
}

template <typename VType>
GridCartesian * Environment::getSliceGrid(const unsigned int orthDir)
{
    FineGridKey key = {typeHash<VType>(), orthDir};

    auto it = grid3d_.find(key);

    if (it != grid3d_.end())
    {
        return it->second.get();
    }
    else
    {
        createSliceGrid<VType>(orthDir);

        return grid3d_.at(key).get();
    }
}

template <typename VType>
GridCartesian * Environment::getGrid(const unsigned int Ls)
{
    if (activeSubGrid_ && subGridScopeOn_)
    {
        HADRONS_ERROR(Logic, "5d subgrid not yet supported in split scope");
    }
    FineGridKey key = {typeHash<VType>(), Ls};
    auto        it  = grid5d_.find(key);

    if (it != grid5d_.end())
    {
        return it->second.get();
    }
    else
    {
        createGrid<VType>(Ls);

        return grid5d_.at(key).get();
    }
}

template <typename VType>
GridRedBlackCartesian * Environment::getRbGrid(const unsigned int Ls)
{
    if (activeSubGrid_ && subGridScopeOn_)
    {
        HADRONS_ERROR(Logic, "5d subgrid not yet supported in split scope");
    }
    FineGridKey key = {typeHash<VType>(), Ls};
    auto        it  = gridRb5d_.find(key);

    if (it != gridRb5d_.end())
    {
        return it->second.get();
    }
    else
    {
        createGrid<VType>(Ls);

        return gridRb5d_.at(key).get();
    }
}

template <typename VType>
GridCartesian * Environment::getCoarseGrid(const std::vector<int> &blockSize,
                                           const unsigned int Ls)
{
    std::vector<int> s = blockSize;

    s.push_back(Ls);

    CoarseGridKey key = {typeHash<VType>(), s};

    auto it = gridCoarse5d_.find(key);
    if (it != gridCoarse5d_.end())
    {
        return it->second.get();
    }
    else
    {
        createCoarseGrid<VType>(blockSize, Ls);

        return gridCoarse5d_.at(key).get();
    }
}


// general memory management ///////////////////////////////////////////////////
template <typename B, typename T, typename ... Ts>
void Environment::createDerivedObject(const std::string name,
                                      const Environment::Storage storage,
                                      const unsigned int Ls,
                                      Ts && ... args)
{
    if (!hasObject(name))
    {
        addObject(name);
    }
    
    unsigned int address = getObjectAddress(name);
    
    // shadow-create mode: a non-lattice global object is being rebuilt onto the
    // subgrid by the VM. Route the new object into the shadow store (keyed by
    // splitKey) so the world-grid object_ entry is left intact for global
    // consumers; getDerivedObject under scope returns this shadow copy.
    if (shadowCreateMode_)
    {
        MemoryStats memStats;
        if (!MemoryProfiler::stats)
        {
            MemoryProfiler::stats = &memStats;
        }
        shadowStore_[{address, shadowCreateKey_}].reset(
            new Holder<B>(new T(std::forward<Ts>(args)...)));
        if (MemoryProfiler::stats == &memStats)
        {
            MemoryProfiler::stats = nullptr;
        }
        return;
    }

    if (!object_[address].data or !objectsProtected())
    {
        MemoryStats memStats;
    
        if (!MemoryProfiler::stats)
        {
            MemoryProfiler::stats = &memStats;
        }
        size_t initMem               = MemoryProfiler::stats->currentlyAllocated;
        object_[address].storage     = storage;
        object_[address].Ls          = Ls;
        object_[address].data.reset(new Holder<B>(new T(std::forward<Ts>(args)...)));
        object_[address].size        = MemoryProfiler::stats->currentlyAllocated - initMem;
        object_[address].type        = typeIdPt<B>();
        object_[address].derivedType = typeIdPt<T>();
        if (MemoryProfiler::stats == &memStats)
        {
            MemoryProfiler::stats = nullptr;
        }
    }
    // object already exists, no error if it is a cache, error otherwise
    else if ((object_[address].storage               != Storage::cache) or 
             (object_[address].storage               != storage)        or
             (object_[address].name                  != name)           or
             (typeHash(object_[address].type)        != typeHash<B>())  or
             (typeHash(object_[address].derivedType) != typeHash<T>()))
    {
        HADRONS_ERROR_REF(ObjectDefinition, "object '" + name + "' already allocated", address);
    }
}

template <typename T, typename ... Ts>
void Environment::createObject(const std::string name, 
                               const Environment::Storage storage,
                               const unsigned int Ls,
                               Ts && ... args)
{
    createDerivedObject<T, T>(name, storage, Ls, std::forward<Ts>(args)...);
}

template <typename B, typename T>
T * Environment::getDerivedObject(const unsigned int address) const
{
    // subgrid scope: check shadow store first (only while scope is on for the
    // current module — a global module must see the world-grid global store)
    if (activeSubGrid_ && subGridScopeOn_)
    {
        auto key = std::make_pair(address, activeSplitKey_);
        auto it  = shadowStore_.find(key);
        if (it != shadowStore_.end())
        {
            if (auto h = dynamic_cast<Holder<B> *>(it->second.get()))
            {
                if (&typeid(T) == &typeid(B))
                {
                    return dynamic_cast<T *>(h->getPt());
                }
                else
                {
                    if (auto hder = dynamic_cast<T *>(h->getPt()))
                    {
                        return hder;
                    }
                    else
                    {
                        HADRONS_ERROR_REF(ObjectType, "object with address " +
                            std::to_string(address) +
                            " cannot be casted to '" + typeName(&typeid(T)) +
                            "' (has type '" + typeName(&typeid(h->get())) + "')", address);
                    }
                }
            }
            else
            {
                HADRONS_ERROR_REF(ObjectType, "object with address " +
                            std::to_string(address) +
                            " does not have type '" + typeName(&typeid(B)) +
                            "' (has type '" + getObjectType(address) + "')", address);
            }
        }
        // no shadow entry: fall through to the global store below
    }
    if (hasObject(address))
    {
        if (hasCreatedObject(address))
        {
            if (auto h = dynamic_cast<Holder<B> *>(object_[address].data.get()))
            {
                if (&typeid(T) == &typeid(B))
                {
                    return dynamic_cast<T *>(h->getPt());
                }
                else
                {
                    if (auto hder = dynamic_cast<T *>(h->getPt()))
                    {
                        return hder;
                    }
                    else
                    {
                        HADRONS_ERROR_REF(ObjectType, "object with address " +
                            std::to_string(address) +
                            " cannot be casted to '" + typeName(&typeid(T)) +
                            "' (has type '" + typeName(&typeid(h->get())) + "')", address);
                    }
                }
            }
            else
            {
                HADRONS_ERROR_REF(ObjectType, "object with address " + 
                            std::to_string(address) +
                            " does not have type '" + typeName(&typeid(B)) +
                            "' (has type '" + getObjectType(address) + "')", address);
            }
        }
        else
        {
            HADRONS_ERROR_REF(ObjectDefinition, "object with address " + 
                              std::to_string(address) + " is empty", address);
        }
    }
    else
    {
        HADRONS_ERROR_REF(ObjectDefinition, "no object with address " + 
                          std::to_string(address), address);
    }
}

template <typename B, typename T>
T * Environment::getDerivedObject(const std::string name) const
{
    return getDerivedObject<B, T>(getObjectAddress(name));
}

template <typename T>
T * Environment::getObject(const unsigned int address) const
{
    return getDerivedObject<T, T>(address);
}

template <typename T>
T * Environment::getObject(const std::string name) const
{
    return getObject<T>(getObjectAddress(name));
}

template <typename T>
bool Environment::isObjectOfType(const unsigned int address) const
{
    if (hasCreatedObject(address))
    {
        if (auto h = dynamic_cast<Holder<T> *>(object_[address].data.get()))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        HADRONS_ERROR_REF(ObjectDefinition, "no initialised object with address " 
                          + std::to_string(address), address);
    }
}

template <typename T>
bool Environment::isObjectOfType(const std::string name) const
{
    return isObjectOfType<T>(getObjectAddress(name));
}

template <typename B, typename T>
bool Environment::isObjectOfDerivedType(const unsigned int address) const
{
    try
    {
        auto o = getDerivedObject<B,T>(address);
    }
    catch(Exceptions::ObjectType)
    {
        return false;
    }
    return true;
}

template <typename B, typename T>
bool Environment::isObjectOfDerivedType(const std::string name) const
{
    return isObjectOfDerivedType<B, T>(getObjectAddress(name));
}

END_HADRONS_NAMESPACE

#endif // Hadrons_Environment_hpp_
