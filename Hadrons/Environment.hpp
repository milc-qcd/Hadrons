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
// lattice element type of a scatterable leaf: Lattice<Vobj> itself, or the
// element of a std::vector<Lattice<Vobj>> (lazy, unlike a conditional_t on
// T::value_type, which Lattice does not define)
template <typename T>
struct LatticeElem { using type = T; };
template <typename Vobj>
struct LatticeElem<std::vector<Lattice<Vobj>>> { using type = Lattice<Vobj>; };
// map-like containers whose MAPPED type is itself scatterable — e.g.
// TGammaMap<PropagatorField> solver guesses or
// TGammaMap<std::vector<PropagatorField>> (n_noise > 1). Detected through
// the nested mapped_type (inherited by std::map subclasses such as
// TGammaMap, so exact-type specializations cannot catch them); the mapped
// type must pass the leaf test, so metadata maps (TGammaMap<StagGamma>)
// correctly stay non-scatterable and keep their cheap rebuild path. Like
// the vector case these are scatter-leaves: the container is
// default-constructed on the subgrid and every mapped value Grid_split
// into place, the producer is NOT rebuilt.
template <typename T, typename = void>
struct IsLatticeMap : std::false_type {};
template <typename T>
struct IsLatticeMap<T, std::void_t<typename T::mapped_type>>
: std::integral_constant<bool,
    IsLattice<typename T::mapped_type>::value
    || IsLatticeVector<typename T::mapped_type>::value> {};
// anything scatterable onto a subgrid (scalar lattice, vector of lattices,
// or map-like container of lattices)
template <typename T>
struct IsScatterable
: std::integral_constant<bool, IsLattice<T>::value
                             || IsLatticeVector<T>::value
                             || IsLatticeMap<T>::value> {};

// Grid-independent metadata that can be cheaply deep-copied into the shadow
// store instead of triggering a full producer rebuild. E.g. std::vector<Integer>
// (time-dilution source shifts), scalar parameters. These contain no Grid
// lattices and no RNG state, so a copy is valid on any subgrid.
template <typename T>
struct IsCloneableMetadata : std::false_type {};
template <typename T>
struct IsCloneableMetadata<std::vector<T>> : std::is_arithmetic<T> {};

// internal-slot -> sub-comm permutation of the batched split /////////////////
// The collapse loop below (mirroring Grid's vector Grid_split) digit-encodes
// the INTERNAL slot index over the SPLIT dimensions with the first-collapsed
// (largest-d) dimension MOST significant; the sub-comm index (splitRank,
// CartesianCommunicator(processors, parent, srank)) digit-encodes the same
// coordinates with the SMALLEST-d dimension most significant
// (IndexFromCoorReversed). One split dimension: the orders coincide and
// slotToSub[p] == p. Two or more: a digit-reversal permutation — e.g. a
// (1,1,2,2) world split into four (1,1,1,1) subgrids delivers slot 2t+z to
// the sub-comm indexed 2z+t (slots 1 and 2 swapped). Rank-independent
// (computed from the two processor layouts alone): every rank derives the
// same map, so permuting the fill cannot itself diverge across ranks.
inline std::vector<unsigned int> gridSplitBatchSlotToSub(
    const GridBase &fullGrid, const GridBase &splitGrid)
{
    const int                 ndim = fullGrid._ndimension;
    std::vector<unsigned int> subExt(ndim), slotStride(ndim), subStride(ndim);
    for (int d = 0; d < ndim; ++d)
    {
        subExt[d] = static_cast<unsigned int>(fullGrid._processors[d])
                    / static_cast<unsigned int>(splitGrid._processors[d]);
    }
    unsigned int acc = 1;
    for (int d = 0; d < ndim; ++d)       // slot digits: larger d MORE
    {                                    // significant (collapse order)
        slotStride[d] = acc;
        acc           *= subExt[d];
    }
    acc = 1;
    for (int d = ndim - 1; d >= 0; --d)  // sub-comm digits: smaller d MORE
    {                                    // significant (reversed lex)
        subStride[d] = acc;
        acc          *= subExt[d];
    }
    std::vector<unsigned int> slotToSub(acc);   // acc == nvector
    for (unsigned int p = 0; p < slotToSub.size(); ++p)
    {
        // decode the internal slot index into per-dim digits (most
        // significant first), then re-encode with the sub-comm significance
        unsigned int              rem = p;
        std::vector<unsigned int> digit(ndim);
        for (int d = ndim - 1; d >= 0; --d)
        {
            digit[d] = rem / slotStride[d];
            rem     %= slotStride[d];
        }
        unsigned int k = 0;
        for (int d = 0; d < ndim; ++d)
        {
            k += digit[d] * subStride[d];
        }
        slotToSub[p] = k;
    }
    return slotToSub;
}

// batched subgrid scatter (R4) //////////////////////////////////////////////
// Pointer-taking reimplementation of Grid's vector Grid_split
// (Grid/Grid/lattice/Lattice_transfer.h:1598-1711). The raw collapse
// address arithmetic (lex_c = c+chunk*m+chunk*M*s) pairs the internal slot
// position with the sub-comm index ONLY for single-dim splits (the two
// digit orders coincide); for multi-dim splits the initial fill below
// permutes the source objects through gridSplitBatchSlotToSub so the
// documented contract — slot k of `fulls` lands on sub-comm k — holds for
// every layout. Unlike the scalar wrapper (:1713-1725, which builds full_v
// from nrhs world-grid COPIES of the same field) the slots here hold
// DISTINCT objects and no world-grid copies are made. All building blocks
// (unvectorizeToLexOrdArray, AllToAll, vectorizeFromLexOrdArray,
// Lexicographic) are public Grid APIs.
template <typename Vobj>
void gridSplitBatch(const std::vector<const Lattice<Vobj> *> &fulls,
                    Lattice<Vobj>                             &split)
{
    typedef typename Vobj::scalar_object Sobj;

    int full_vecs   = fulls.size();

    GRID_ASSERT(full_vecs >= 1);

    GridBase * full_grid = fulls[0]->Grid();
    GridBase *split_grid = split.Grid();

    int       ndim  = full_grid->_ndimension;
    int  full_nproc = full_grid->_Nprocessors;
    int split_nproc =split_grid->_Nprocessors;

    // checkerboard management
    int cb = fulls[0]->Checkerboard();
    split.Checkerboard() = cb;

    // checks (mirror Grid's vector overload)
    GRID_ASSERT(full_grid->_ndimension == split_grid->_ndimension);
    for (int n = 0; n < full_vecs; n++)
    {
        GRID_ASSERT(fulls[n]->Checkerboard() == cb);
        for (int d = 0; d < ndim; d++)
        {
            GRID_ASSERT(fulls[n]->Grid()->_gdimensions[d]
                        == split.Grid()->_gdimensions[d]);
            GRID_ASSERT(fulls[n]->Grid()->_fdimensions[d]
                        == split.Grid()->_fdimensions[d]);
        }
    }

    int   nvector = full_nproc / split_nproc;
    GRID_ASSERT(nvector * split_nproc == full_nproc);
    GRID_ASSERT(nvector == full_vecs);

    Coordinate ratio(ndim);
    for (int d = 0; d < ndim; d++)
    {
        ratio[d] = full_grid->_processors[d] / split_grid->_processors[d];
    }

    uint64_t lsites = full_grid->lSites();
    uint64_t     sz = lsites * nvector;
    std::vector<Sobj> tmpdata(sz);
    std::vector<Sobj> alldata(sz);
    std::vector<Sobj> scalardata(lsites);

    // slot k of `fulls` must land on sub-comm k (the documented contract
    // the driver relies on). The collapse arithmetic routes internal slot
    // position p to sub-comm slotToSub[p], so object k rides position
    // subToSlot[k]. One split dimension: slotToSub is the identity and this
    // loop is byte-for-byte the old fill (single-dim runs stay
    // bitwise-identical).
    const auto slotToSub = gridSplitBatchSlotToSub(*full_grid, *split_grid);
    std::vector<unsigned int> subToSlot(slotToSub.size());
    for (unsigned int p = 0; p < slotToSub.size(); ++p)
    {
        subToSlot[slotToSub[p]] = p;
    }
    for (unsigned int k = 0; k < static_cast<unsigned int>(nvector); ++k)
    {
        unvectorizeToLexOrdArray(scalardata, *fulls[k]);
        thread_for(site, lsites, {
            alldata[subToSlot[k] * lsites + site] = scalardata[site];
        });
    }

    int      nvec  = nvector;   // counts down to 1 as we collapse dims
    Coordinate ldims = full_grid->_ldimensions;

    for (int d = ndim - 1; d >= 0; d--)
    {
        if (ratio[d] != 1)
        {
            full_grid->AllToAll(d, alldata, tmpdata);
            if (split_grid->_processors[d] > 1)
            {
                alldata = tmpdata;
                split_grid->AllToAll(d, alldata, tmpdata);
            }

            auto  rdims = ldims;
            auto      M = ratio[d];
            auto rsites = lsites * M;   // increases rsites by M
            nvec       /= M;            // reduce nvec by subdivision factor
            rdims[d]   *= M;            // increase local dim by same factor

            int sP = split_grid->_processors[d];
            int fP =  full_grid->_processors[d];

            int fvol  = lsites;
            int chunk = (nvec * fvol) / sP;
            GRID_ASSERT(chunk * sP == nvec * fvol);

            thread_for(c, chunk, {
                Coordinate coor(ndim);
                for (int m = 0; m < M; m++)
                {
                    for (int s = 0; s < sP; s++)
                    {
                        int lex_r;
                        uint64_t lex_c        = c + chunk * m + chunk * M * s;
                        uint64_t lex_fvol_vec = c + chunk * s;
                        uint64_t lex_fvol     = lex_fvol_vec % fvol;
                        uint64_t lex_vec      = lex_fvol_vec / fvol;

                        Lexicographic::CoorFromIndex(coor, lex_fvol, ldims);
                        coor[d] += m * ldims[d];
                        Lexicographic::IndexFromCoor(coor, lex_r, rdims);
                        lex_r += lex_vec * rsites;

                        alldata[lex_r] = tmpdata[lex_c];
                    }
                }
            });
            ldims[d] *= ratio[d];
            lsites   *= ratio[d];
        }
    }
    vectorizeFromLexOrdArray(alldata, split);
}

// gather context: this rank's sub-comm index, the producing sub-comm index,
// and the WORLD rank of the producing sub-comm's first rank (root for the
// shape/checkerboard descriptor broadcast — non-owning ranks hold no data
// at the address yet must join Grid_unsplit with rank-uniform leaf counts,
// keys and checkerboard).
struct GatherCtx
{
    int me{-1};        // this rank's sub-comm index
    int sourceSub{-1}; // sub-comm whose piece is the gather source
    int bossRank{-1};  // world rank of sourceSub's first rank
};

// broadcast the variable-size gather descriptor (int32 payload) from the
// producing sub-comm's boss rank to every world rank. The boss passes a
// packed buffer; all other ranks pass an empty one and receive. Collective
// — every rank must call it (the phase-exit gather driver calls it in
// lockstep, like every world collective in the split machinery).
inline void broadcastGatherDescriptor(std::vector<int32_t> &desc,
                                      const GatherCtx      &ctx)
{
    if (ctx.bossRank < 0)
    {
        HADRONS_ERROR(Logic, "broadcastGatherDescriptor: invalid boss rank");
    }
    uint64_t n = desc.size();
    CartesianCommunicator::BroadcastWorld(ctx.bossRank, &n, sizeof(n));
    desc.resize(n);
    if (n > 0)
    {
        CartesianCommunicator::BroadcastWorld(ctx.bossRank, desc.data(),
                                              n * sizeof(int32_t));
    }
}

// leaf-level gather: assemble the world-grid twin of ONE lattice leaf from
// this rank's subgrid piece via Grid_unsplit (the exact inverse of the
// Grid_split behind splitTo). `fulls` is the caller-owned scratch buffer of
// nrhs world-grid lattices, REUSED across the leaves of one object; slot
// subToSlot[sourceSub] receives the producing sub-comm's assembly, is moved
// out into `twin`, and is re-constructed for the next leaf. `piece` is the
// real leaf on ranks of the producing sub-comm and a zero-filled dummy
// elsewhere: Grid_unsplit never reads destination contents, so non-owning
// contributions only ever land in discarded slots.
template <typename LeafT>
void gridGatherLeaf(std::vector<LeafT>              &fulls,
                    LeafT                           *piece,
                    LeafT                           &twin,
                    GridCartesian                   *world,
                    GridCartesian                   *sub,
                    const std::vector<unsigned int> &subToSlot,
                    const unsigned int               sourceSub,
                    const int                        cb)
{
    LeafT dummy(sub);
    dummy = Zero();   // discarded slots only; zeroed for determinism
    LeafT &src = (piece != nullptr) ? *piece : dummy;
    if (piece != nullptr && static_cast<int>(piece->Checkerboard()) != cb)
    {
        HADRONS_ERROR(Definition, "subgrid gather: piece checkerboard "
                      "mismatch against the broadcast descriptor");
    }
    for (auto &f : fulls)
    {
        f.Checkerboard() = cb;   // Grid_unsplit propagates full[0]'s cb
    }
    Grid_unsplit(fulls, src);
    const unsigned int slot = subToSlot[sourceSub];
    twin = std::move(fulls[slot]);
    fulls[slot] = LeafT(world);   // re-arm the moved-from scratch slot
}

// local byte footprint of one subgrid lattice on this rank: oSites * sizeof
// (vobj) is exactly what the aligned allocator reserved (mirrors the
// profiler-delta accounting of world objects closely enough for runtime
// visibility; rebuilt non-lattice shadows record their exact profiler
// delta at creation instead — see createDerivedObject)
template <typename Vobj>
SITE_SIZE_TYPE latticeFootprint(const Lattice<Vobj> &lat)
{
    return static_cast<SITE_SIZE_TYPE>(lat.Grid()->oSites()) * sizeof(Vobj);
}

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
    // lightweight (non-mutating) lattice test — true iff the held object is
    // a scatterable lattice value or a container of lattice values (vector
    // or map-like). Used by the GC scheduler to tell scatter-leaf objects
    // (lattices and lattice containers) from rebuildable objects
    // (solvers/actions) without performing a Grid_split.
    virtual bool isLattice(void) const
    {
        return false;
    }
    // lightweight (non-mutating) cloneability test — true iff the held
    // object is grid-independent metadata safe to deep-copy into the
    // shadow store (IsCloneableMetadata). Mirrors isLattice(): lets the
    // VM's consumer-set walk classify objects without performing the
    // clone.
    virtual bool isCloneable(void) const
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
    // Batched variant of splitTo for sole-consumer objects (R4): `fulls`
    // holds nrhs type-erased objects of the SAME concrete scatterable type
    // (the dispatcher guarantees this); slot k must land on sub-comm k.
    // Returns this rank's piece (of the slot matching its sub-comm) as a
    // new Object, or nullptr when shapes are incompatible (the caller falls
    // back to per-object scatters). Default: not supported.
    virtual std::unique_ptr<Object> splitBatchTo(
        const std::vector<const Object *> &fulls) const
    {
        return nullptr;
    }
    // type-erased gather back to the world grid (the reverse of splitTo):
    // assembles this rank's subgrid piece of the held object into a
    // world-grid twin via Grid_unsplit. Called on the REAL holder where this
    // rank holds data, and on a default-constructed ghost Holder otherwise
    // (non-owning ranks hold no data but must join the collective) — see
    // gatherDispatch<B> and Environment::gatherObject. Returns nullptr for
    // non-scatterable types.
    virtual std::unique_ptr<Object> unsplitTo(const GatherCtx &ctx) const
    {
        return nullptr;
    }
    // approximate local byte footprint of the held object on this rank
    // (lattices: oSites*sizeof(vobj); summed element/key-wise for vector
    // and map containers; 0 for non-lattice holders — rebuilt shadows
    // record their exact profiler-measured size at creation instead).
    // Used by addShadowObject so getTotalSize/memPeak/printMemory see
    // shadow memory (R9, runtime-only visibility).
    virtual SITE_SIZE_TYPE footprint(void) const
    {
        return 0;
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
    std::unique_ptr<Object> splitBatchTo(
        const std::vector<const Object *> &fulls) const override;
    std::unique_ptr<Object> unsplitTo(const GatherCtx &ctx) const override;
    SITE_SIZE_TYPE footprint(void) const override;
    bool      isLattice(void) const override
    {
        return IsScatterable<T>::value;
    }
    bool      isCloneable(void) const override
    {
        return IsCloneableMetadata<T>::value;
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

// type-erased gather dispatch: routes through the holder's virtual unsplitTo
// when this rank holds the object; otherwise a default-constructed ghost
// Holder<B> carries the static type for the dummy-piece path. Registered per
// object at createDerivedObject time (see below) — the memory-profile pass
// creates every program object world-scope on every rank, so registration
// is rank-uniform and survives frees exactly like type/derivedType do.
template <typename B>
std::unique_ptr<Object> gatherDispatch(Object *piece, const GatherCtx &ctx)
{
    if (piece != nullptr)
    {
        return piece->unsplitTo(ctx);
    }
    Holder<B> ghost;

    return ghost.unsplitTo(ctx);
}

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
    typedef std::unique_ptr<Object> (*GatherDispatch)(Object *,
                                                      const GatherCtx &);
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
        // type-erased gather dispatch (subgrid gather, the reverse of
        // splitTo): registered at createDerivedObject time on every rank
        // that ever created the object. Survives frees exactly like
        // type/derivedType above. Null when never registered.
        GatherDispatch gather{nullptr};
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
    GridCartesian *         getWorldGrid(void);
    template <typename VType = vComplex>
    GridCartesian *         getSubGrid(void);
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
    // free ONE shadow-store entry (local to the owning sub-comm; no
    // collectives -- consumer sets are rank-independent, so every owning
    // rank derives the same free step). Returns true if an entry existed.
    bool                    freeShadowObject(const unsigned int address,
                                             const int splitKey);
    // erase every shadow-store entry with the given split key (end-of-phase
    // sweep / safety net). Local; returns the number of entries erased.
    unsigned int            sweepShadows(const int splitKey);
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
    // non-storing variant of scatterObject: perform the split and return
    // the scattered copy WITHOUT adding it to the shadow store (the VM
    // decides per-rank whether to store — consumer-gated shadow storage).
    // Returns nullptr for non-scatterable objects.
    std::unique_ptr<Object>  splitObject(const unsigned int address);
    // batched sole-consumer scatter (R4): `addresses` holds nrhs objects of
    // the same concrete scatterable type, slot k consumed only by sub-comm
    // k. Dispatches on the first element's concrete type; returns this
    // rank's piece WITHOUT storing (the VM gates the store per consumer),
    // or nullptr when the batch is shape-incompatible.
    std::unique_ptr<Object> splitObjectBatch(
        const std::vector<unsigned int> &addresses);
    // subgrid gather (the reverse of scatterObject): assemble a
    // tagged-phase-produced object back into a world-grid twin at the SAME
    // address, on every rank, via the registered gather dispatch and
    // Grid_unsplit. COLLECTIVE — every world rank must call it in lockstep
    // (the VM does, at phase exit). Returns true when the twin was
    // installed; false when the address has no registered dispatch (the
    // caller hard-errors: the object's type cannot be gathered).
    bool                    gatherObject(const unsigned int address,
                                         const GatherCtx   &ctx);
    // lightweight cloneability probe (see Object::isCloneable)
    bool                    isCloneableObject(const unsigned int address) const;
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
    // shadow memory accounting (R9, runtime-only): local byte size of each
    // shadow entry — recorded at store time (container footprint for
    // scattered pieces; profiler-measured delta for rebuilt objects in the
    // shadow-create branch of createDerivedObject), erased on free/sweep
    // and at clearActiveSubGrid. Per-rank local, like object_[].size.
    std::map<std::pair<unsigned int, int>, Size> shadowSize_;
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
    else if constexpr (IsLatticeMap<T>::value)
    {
        // map-like container of scatterable lattice values (e.g. a
        // TGammaMap<PropagatorField> guess map): default-construct the
        // container (TGammaMap is `= default`) and Grid_split every mapped
        // value onto the subgrid, mirroring the element-wise vector branch.
        // Both value shapes are covered: Lattice<Vobj> mapped values and
        // std::vector<Lattice<Vobj>> mapped values (n_noise > 1 guesses).
        using MappedT = typename T::mapped_type;
        using ElemT   = typename LatticeElem<MappedT>::type;
        using VType   = typename ElemT::vector_type;
        auto  &env    = Environment::getInstance();
        auto  *sg     = env.template getGrid<VType>();
        auto   split  = std::make_unique<T>();
        for (auto &p : *objPt_)
        {
            if constexpr (IsLattice<MappedT>::value)
            {
                auto ret = split->emplace(p.first, sg);
                Grid_split(p.second, ret.first->second);
            }
            else
            {
                auto ret = split->try_emplace(p.first);
                ret.first->second.reserve(p.second.size());
                for (auto &e : p.second)
                {
                    ret.first->second.emplace_back(sg);
                    Grid_split(e, ret.first->second.back());
                }
            }
        }
        return std::make_unique<Holder<T>>(split.release());
    }
    else
    {
        return nullptr;
    }
}

// subgrid gather (the reverse of splitTo) /////////////////////////////////////
template <typename T>
std::unique_ptr<Object> Holder<T>::unsplitTo(const GatherCtx &ctx) const
{
    if constexpr (IsLattice<T>::value)
    {
        using VType = typename T::vector_type;
        auto  &env   = Environment::getInstance();
        auto  *world = env.template getWorldGrid<VType>();
        auto  *sub   = env.template getSubGrid<VType>();
        const auto  slotToSub = gridSplitBatchSlotToSub(*world, *sub);
        std::vector<unsigned int> subToSlot(slotToSub.size());
        for (unsigned int p = 0; p < slotToSub.size(); ++p)
        {
            subToSlot[slotToSub[p]] = p;
        }
        // shape/checkerboard descriptor: [cb] (a default-constructed
        // lattice's cb is 0)
        std::vector<int32_t> desc;
        if (CartesianCommunicator::RankWorld() == ctx.bossRank)
        {
            desc.push_back(static_cast<int32_t>((*objPt_).Checkerboard()));
        }
        broadcastGatherDescriptor(desc, ctx);
        auto twin = std::make_unique<T>(world);
        std::vector<T> fulls(slotToSub.size(), T(world));
        gridGatherLeaf(fulls, objPt_.get(), *twin, world, sub, subToSlot,
                       static_cast<unsigned int>(ctx.sourceSub),
                       static_cast<int>(desc[0]));
        return std::make_unique<Holder<T>>(twin.release());
    }
    else if constexpr (IsLatticeVector<T>::value)
    {
        using ElemT  = typename T::value_type;
        using VType  = typename ElemT::vector_type;
        auto  &env   = Environment::getInstance();
        auto  *world = env.template getWorldGrid<VType>();
        auto  *sub   = env.template getSubGrid<VType>();
        const auto  slotToSub = gridSplitBatchSlotToSub(*world, *sub);
        std::vector<unsigned int> subToSlot(slotToSub.size());
        for (unsigned int p = 0; p < slotToSub.size(); ++p)
        {
            subToSlot[slotToSub[p]] = p;
        }
        // descriptor: [cb, element count] (empty vectors: cb 0, count 0)
        std::vector<int32_t> desc;
        if (CartesianCommunicator::RankWorld() == ctx.bossRank)
        {
            desc.push_back(objPt_->empty()
                           ? 0
                           : static_cast<int32_t>((*objPt_)[0].Checkerboard()));
            desc.push_back(static_cast<int32_t>(objPt_->size()));
        }
        broadcastGatherDescriptor(desc, ctx);
        const auto n = static_cast<typename T::size_type>(desc[1]);
        if (objPt_ && objPt_->size() != n)
        {
            HADRONS_ERROR(Definition, "subgrid gather: vector length "
                          "divergence between producing ranks");
        }
        auto twin = std::make_unique<T>();
        twin->reserve(n);
        std::vector<ElemT> fulls(slotToSub.size(), ElemT(world));
        for (typename T::size_type e = 0; e < n; ++e)
        {
            twin->emplace_back(world);
            gridGatherLeaf(fulls, objPt_ ? &(*objPt_)[e] : nullptr,
                           twin->back(), world, sub, subToSlot,
                           static_cast<unsigned int>(ctx.sourceSub),
                           static_cast<int>(desc[0]));
        }
        return std::make_unique<Holder<T>>(twin.release());
    }
    else if constexpr (IsLatticeMap<T>::value)
    {
        using MappedT = typename T::mapped_type;
        using KeyT    = typename T::key_type;
        if constexpr (!(std::is_enum<KeyT>::value
                        || std::is_integral<KeyT>::value))
        {
            return nullptr;   // map keys not wire-broadcastable: unsupported
        }
        else
        {
            using ElemT  = typename LatticeElem<MappedT>::type;
            using VType  = typename ElemT::vector_type;
            auto  &env   = Environment::getInstance();
            auto  *world = env.template getWorldGrid<VType>();
            auto  *sub   = env.template getSubGrid<VType>();
            const auto  slotToSub = gridSplitBatchSlotToSub(*world, *sub);
            std::vector<unsigned int> subToSlot(slotToSub.size());
            for (unsigned int p = 0; p < slotToSub.size(); ++p)
            {
                subToSlot[slotToSub[p]] = p;
            }
            // descriptor: [cb, key count, (key, per-key vector length)*]
            // (length 0/unused for lattice-valued mapped types; std::map
            // iterates keys in sorted order, identical on every producing
            // rank because the object was produced by one module run)
            std::vector<int32_t> desc;
            if (CartesianCommunicator::RankWorld() == ctx.bossRank)
            {
                if constexpr (IsLattice<MappedT>::value)
                {
                    desc.push_back(objPt_->empty() ? 0
                        : static_cast<int32_t>(
                            objPt_->begin()->second.Checkerboard()));
                }
                else
                {
                    desc.push_back(objPt_->empty() ? 0
                        : static_cast<int32_t>(
                            objPt_->begin()->second[0].Checkerboard()));
                }
                desc.push_back(static_cast<int32_t>(objPt_->size()));
                for (auto &p : *objPt_)
                {
                    desc.push_back(static_cast<int32_t>(p.first));
                    if constexpr (IsLattice<MappedT>::value)
                    {
                        desc.push_back(0);
                    }
                    else
                    {
                        desc.push_back(static_cast<int32_t>(p.second.size()));
                    }
                }
            }
            broadcastGatherDescriptor(desc, ctx);
            const int      cb = static_cast<int>(desc[0]);
            const auto     nk = static_cast<std::size_t>(desc[1]);
            if (objPt_ && objPt_->size() != nk)
            {
                HADRONS_ERROR(Definition, "subgrid gather: map key set "
                              "divergence between producing ranks");
            }
            auto twin = std::make_unique<T>();
            std::vector<ElemT> fulls(slotToSub.size(), ElemT(world));
            for (std::size_t i = 0; i < nk; ++i)
            {
                const auto key = static_cast<KeyT>(desc[2 + 2*i]);
                const auto len = static_cast<std::size_t>(desc[3 + 2*i]);
                if constexpr (IsLattice<MappedT>::value)
                {
                    auto ret = twin->emplace(key, world);
                    gridGatherLeaf(fulls, objPt_ ? &(objPt_->at(key)) : nullptr,
                                   ret.first->second, world, sub, subToSlot,
                                   static_cast<unsigned int>(ctx.sourceSub),
                                   cb);
                }
                else
                {
                    auto ret = twin->try_emplace(key);
                    if (objPt_ && objPt_->at(key).size() != len)
                    {
                        HADRONS_ERROR(Definition, "subgrid gather: mapped "
                                      "vector length divergence between "
                                      "producing ranks");
                    }
                    ret.first->second.reserve(len);
                    for (std::size_t e = 0; e < len; ++e)
                    {
                        ret.first->second.emplace_back(world);
                        gridGatherLeaf(fulls,
                                       objPt_ ? &(objPt_->at(key)[e]) : nullptr,
                                       ret.first->second.back(), world, sub,
                                       subToSlot,
                                       static_cast<unsigned int>(ctx.sourceSub),
                                       cb);
                    }
                }
            }
            return std::make_unique<Holder<T>>(twin.release());
        }
    }
    else
    {
        return nullptr;
    }
}

// approximate local byte footprint (R9): lattice element/key-wise
// oSites*sizeof(vobj); 0 for non-lattice holders
template <typename T>
SITE_SIZE_TYPE Holder<T>::footprint(void) const
{
    if constexpr (IsLattice<T>::value)
    {
        return latticeFootprint(*objPt_);
    }
    else if constexpr (IsLatticeVector<T>::value)
    {
        SITE_SIZE_TYPE sz = 0;
        for (auto &e : *objPt_)
        {
            sz += latticeFootprint(e);
        }
        return sz;
    }
    else if constexpr (IsLatticeMap<T>::value)
    {
        SITE_SIZE_TYPE sz = 0;
        for (auto &p : *objPt_)
        {
            if constexpr (IsLattice<typename T::mapped_type>::value)
            {
                sz += latticeFootprint(p.second);
            }
            else
            {
                for (auto &e : p.second)
                {
                    sz += latticeFootprint(e);
                }
            }
        }
        return sz;
    }
    else
    {
        return 0;
    }
}

// batched subgrid scatter (R4) //////////////////////////////////////////////
template <typename T>
std::unique_ptr<Object> Holder<T>::splitBatchTo(
    const std::vector<const Object *> &fulls) const
{
    if constexpr (IsScatterable<T>::value)
    {
        auto &env = Environment::getInstance();
        // every slot must hold the same concrete Holder<T>; any mismatch
        // (or shape incompatibility below) aborts the whole batch and the
        // caller falls back to per-object scatters
        std::vector<const Holder<T> *> slots;
        slots.reserve(fulls.size());
        for (auto *o : fulls)
        {
            auto *h = dynamic_cast<const Holder<T> *>(o);
            if (h == nullptr)
            {
                return nullptr;
            }
            slots.push_back(h);
        }
        if constexpr (IsLattice<T>::value)
        {
            using VType = typename T::vector_type;
            auto  *sg   = env.template getGrid<VType>();
            auto   split = std::make_unique<T>(sg);
            std::vector<const T *> ptrs;
            for (auto *h : slots)
            {
                ptrs.push_back(h->getPt());
            }
            gridSplitBatch(ptrs, *split);
            return std::make_unique<Holder<T>>(split.release());
        }
        else if constexpr (IsLatticeVector<T>::value)
        {
            using ElemT  = typename T::value_type;         // Lattice<Vobj>
            using VType  = typename ElemT::vector_type;
            const auto   n = slots[0]->getPt()->size();
            for (auto *h : slots)
            {
                if (h->getPt()->size() != n)
                {
                    return nullptr;   // unequal vector sizes: no batch
                }
            }
            auto  *sg    = env.template getGrid<VType>();
            auto   split = std::make_unique<T>();
            split->reserve(n);
            std::vector<const ElemT *> ptrs(slots.size(), nullptr);
            for (typename T::size_type e = 0; e < n; ++e)
            {
                for (std::size_t k = 0; k < slots.size(); ++k)
                {
                    ptrs[k] = &(*slots[k]->getPt())[e];
                }
                split->emplace_back(sg);
                gridSplitBatch(ptrs, split->back());
            }
            return std::make_unique<Holder<T>>(split.release());
        }
        else if constexpr (IsLatticeMap<T>::value)
        {
            // identical key sequences required (std::map subclasses
            // iterate in key order; compare slot 0 against every slot)
            using MappedT = typename T::mapped_type;
            using ElemT   = typename LatticeElem<MappedT>::type;
            using VType   = typename ElemT::vector_type;
            auto &m0      = *slots[0]->getPt();
            for (auto *h : slots)
            {
                auto &m = *h->getPt();
                if (m.size() != m0.size())
                {
                    return nullptr;
                }
                auto it0 = m0.begin();
                auto it  = m.begin();
                for (; it0 != m0.end(); ++it0, ++it)
                {
                    if (it->first != it0->first)
                    {
                        return nullptr;
                    }
                }
            }
            auto *sg    = env.template getGrid<VType>();
            auto  split = std::make_unique<T>();
            std::vector<const ElemT *> ptrs(slots.size(), nullptr);
            for (auto &p0 : m0)
            {
                if constexpr (IsLattice<MappedT>::value)
                {
                    auto ret = split->emplace(p0.first, sg);
                    for (std::size_t k = 0; k < slots.size(); ++k)
                    {
                        ptrs[k] = &slots[k]->getPt()->at(p0.first);
                    }
                    gridSplitBatch(ptrs, ret.first->second);
                }
                else
                {
                    auto ret = split->try_emplace(p0.first);
                    const auto n = p0.second.size();
                    for (auto *h : slots)
                    {
                        if (h->getPt()->at(p0.first).size() != n)
                        {
                            return nullptr;   // mapped vector sizes differ
                        }
                    }
                    ret.first->second.reserve(n);
                    for (typename MappedT::size_type e = 0; e < n; ++e)
                    {
                        for (std::size_t k = 0; k < slots.size(); ++k)
                        {
                            ptrs[k] = &slots[k]->getPt()->at(p0.first)[e];
                        }
                        ret.first->second.emplace_back(sg);
                        gridSplitBatch(ptrs, ret.first->second.back());
                    }
                }
            }
            return std::make_unique<Holder<T>>(split.release());
        }
        else
        {
            return nullptr;
        }
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

// world grid REGARDLESS of sub-grid scope: getGrid diverts to the sub-grid
// cache when scope is on, but world-lattice builders (e.g. the batch
// slot-delivery probe's markers, Phase 12 — built at pass 2 of the batched
// pre-scatter where scope IS on) need the world layout.
template <typename VType>
GridCartesian * Environment::getWorldGrid(void)
{
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

// per-VType subgrid cache REGARDLESS of scope: the gather runs at phase
// exit with scope OFF, but the pieces live on gridSub4d_ keyed by VType
// (the cache survives until the split key changes / clearActiveSubGrid —
// both are strictly after the gather). Create-on-miss requires an active
// subgrid, which spans the whole phase. Mirrors getWorldGrid's
// scope-freeness (the dcbe257a lesson: never build world/subgrid lattices
// through the scope-aware accessor at transitional points).
template <typename VType>
GridCartesian * Environment::getSubGrid(void)
{
    FineGridKey key = {typeHash<VType>(), 1};

    auto it = gridSub4d_.find(key);
    if (it != gridSub4d_.end())
    {
        return it->second.get();
    }
    if (activeSubGrid_ == nullptr)
    {
        HADRONS_ERROR(Logic, "getSubGrid: no active subgrid (gather "
                      "outside a split phase?)");
    }
    createSubGrid<VType>();
    return gridSub4d_.at(key).get();
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
        const size_t initMem = MemoryProfiler::stats->currentlyAllocated;
        shadowStore_[{address, shadowCreateKey_}].reset(
            new Holder<B>(new T(std::forward<Ts>(args)...)));
        // R9 (Slice 7): record the newly computed allocation delta —
        // `initMem` is NEW content (HEAD's branch computes no baseline) —
        // the natural exact size hook for rebuilt (non-lattice) shadows
        shadowSize_[{address, shadowCreateKey_}]
            = MemoryProfiler::stats->currentlyAllocated - initMem;
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
        // subgrid gather dispatch: records how to gather this object type back
        // to the world grid even on ranks where the object is later absent
        // (tagged-phase outputs exist only on their owning sub-comm's ranks;
        // the profiling pass registers this on every rank)
        object_[address].gather = &gatherDispatch<B>;
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
                        HADRONS_ERROR_REF(ObjectType, "object '" + object_[address].name +
                            "' (address " + std::to_string(address) +
                            ") cannot be casted to '" + typeName(&typeid(T)) +
                            "' (has type '" + typeName(&typeid(h->get())) + "')", address);
                    }
                }
            }
            else
            {
                HADRONS_ERROR_REF(ObjectType, "object '" + object_[address].name +
                            "' (address " + std::to_string(address) +
                            ") does not have type '" + typeName(&typeid(B)) +
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
                        HADRONS_ERROR_REF(ObjectType, "object '" + object_[address].name +
                            "' (address " + std::to_string(address) +
                            ") cannot be casted to '" + typeName(&typeid(T)) +
                            "' (has type '" + typeName(&typeid(h->get())) + "')", address);
                    }
                }
            }
            else
            {
                HADRONS_ERROR_REF(ObjectType, "object '" + object_[address].name + "' (address " + 
                            std::to_string(address) +
                            ") does not have type '" + typeName(&typeid(B)) +
                            "' (has type '" + getObjectType(address) + "')", address);
            }
        }
        else
        {
            HADRONS_ERROR_REF(ObjectDefinition, "object '" + object_[address].name + "' (address " + 
                              std::to_string(address) + ") is empty", address);
        }
    }
    else
    {
        HADRONS_ERROR_REF(ObjectDefinition, "no object with address " + 
                          std::to_string(address) + " (out of range)", address);
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
        HADRONS_ERROR_REF(ObjectDefinition, "no initialised object '" + (hasObject(address) ? object_[address].name : std::string("(unknown)")) + "' (address " + std::to_string(address) + ")", address);
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
