/*
 * VirtualMachine.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Andrew Zhen Ning Yong <andrew.yong@ed.ac.uk>
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fabian Joswig <fabian.joswig@ed.ac.uk>
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

#ifndef Hadrons_VirtualMachine_hpp_
#define Hadrons_VirtualMachine_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Database.hpp>
#include <Hadrons/Graph.hpp>
#include <Hadrons/Environment.hpp>
#include <limits>
#include <set>

BEGIN_HADRONS_NAMESPACE

#define DEFINE_VM_ALIAS \
inline VirtualMachine & vm(void) const\
{\
    return VirtualMachine::getInstance();\
}

/******************************************************************************
 *                   Virtual machine for module execution                     *
 ******************************************************************************/
// forward declaration of Module
class ModuleBase;

class VirtualMachine
{
    SINGLETON_DEFCTOR(VirtualMachine);
public:
    typedef SITE_SIZE_TYPE                      Size;
    typedef std::unique_ptr<ModuleBase>         ModPt;
    typedef std::vector<std::set<unsigned int>> GarbageSchedule;
    typedef std::vector<unsigned int>           Program;
    // inclusive step range of one split phase within a program; `phase` is the
    // raw <splitPhase> value (sparse indices allowed; sequence must be
    // non-decreasing over tagged modules — see validateSplitPhases)
    struct PhaseRange
    {
        unsigned int phase{0};
        unsigned int first{0}, last{0};
    };
    struct MemoryPrint
    {
        Size                 size;
        Environment::Storage storage;
        int                  module;
    };
    struct MemoryProfile
    {
        std::vector<std::map<unsigned int, Size>> module;
        std::vector<MemoryPrint>                  object;
    };
    class GeneticPar: Serializable
    {
    public:
        GeneticPar(void):
            popSize{20}, maxGen{1000}, maxCstGen{100}, mutationRate{.1} {};
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GeneticPar,
                                        unsigned int, popSize,
                                        unsigned int, maxGen,
                                        unsigned int, maxCstGen,
                                        double      , mutationRate);
    };

    // serializable classes for database entries
    struct GlobalEntry: SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlUnique<SqlNotNull<std::string>>, name,
                           SqlNotNull<std::string>           , value);
    };

    struct ModuleEntry: SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlUnique<SqlNotNull<unsigned int>>, moduleId,
                           SqlUnique<SqlNotNull<std::string>> , name,
                           SqlNotNull<unsigned int>           , moduleTypeId,
                           std::string                        , parameters);
    };

    struct ModuleTypeEntry: SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlUnique<unsigned int>           , moduleTypeId,
                           SqlUnique<SqlNotNull<std::string>>, type);
                           
    };

    struct ObjectEntry: SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlUnique<SqlNotNull<unsigned int>>, objectId,
                           SqlUnique<SqlNotNull<std::string>> , name,
                           SqlNotNull<unsigned int>           , objectTypeId,
                           SqlNotNull<SITE_SIZE_TYPE>         , size,
                           SqlNotNull<Environment::Storage>   , storageType,
                           SqlNotNull<unsigned int>           , moduleId);
    };

    struct ObjectTypeEntry: SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlUnique<unsigned int>           , objectTypeId,
                           SqlUnique<SqlNotNull<std::string>>, type,
                           SqlNotNull<std::string>, baseType);
    };

    struct ScheduleEntry: SqlEntry
    {
        HADRONS_SQL_FIELDS(SqlUnique<SqlNotNull<unsigned int>>, step,
                           SqlUnique<SqlNotNull<unsigned int>>, moduleId);
    };
private:
    struct ModuleInfo
    {
        const std::type_info      *type{nullptr};
        std::string               name;
        ModPt                     data{nullptr};
        std::vector<unsigned int> input, output;
        size_t                    maxAllocated;
        int                       subgrid{-1};   // subgrid index, -1 = global
        int                       splitPhase{0}; // split phase index, 0 = first
                                                 // phase (meaningful only when
                                                 // subgrid >= 0)
    };
public:
    // trajectory counter
    void                setTrajectory(const unsigned int traj);
    unsigned int        getTrajectory(void) const;
    // run tag
    void                setRunId(const std::string id);
    std::string         getRunId(void) const;
    // database
    void                setDatabase(Database &db);
    void                dbRestoreMemoryProfile(void);
    void                dbRestoreModules(void);
    Program             dbRestoreSchedule(void);
    // module management
    void                pushModule(ModPt &pt, const int subgrid = -1,
                                   const int splitPhase = 0);
    template <typename M>
    void                createModule(const std::string name);
    template <typename M>
    void                createModule(const std::string name,
                                         const typename M::Par &par);
    void                createModule(const std::string name,
                                     const std::string type,
                                     XmlReader &reader,
                                     const std::string blockName = "options",
                                     const int subgrid = -1,
                                     const int splitPhase = 0);
    unsigned int        getNModule(void) const;
    ModuleBase *        getModule(const unsigned int address) const;
    ModuleBase *        getModule(const std::string name) const;
    template <typename M>
    M *                 getModule(const unsigned int address) const;
    template <typename M>
    M *                 getModule(const std::string name) const;
    unsigned int        getModuleAddress(const std::string name) const;
    std::string         getModuleName(const unsigned int address) const;
    std::string         getModuleType(const unsigned int address) const;
    std::string         getModuleType(const std::string name) const;
    std::string         getModuleNamespace(const unsigned int address) const;
    std::string         getModuleNamespace(const std::string name) const;
    int                 getModuleSubgrid(const unsigned int address) const;
    int                 getModuleSplitPhase(const unsigned int address) const;
    void                buildSubGrids(const std::vector<int> &mpiSplit);
    void                computeSplitPhase(void);
    int                 getMe(void) const { return me_; }
    int                 getCurrentModule(void) const;
    bool                hasModule(const unsigned int address) const;
    bool                hasModule(const std::string name) const;
    // print VM content
    void                printContent(void) const;
    // module graph (could be a const reference if topoSort was const)
    Graph<unsigned int> getModuleGraph(void);
    // dump GraphViz graph
    void                dumpModuleGraph(std::ostream &out);
    void                dumpModuleGraph(void);
    void                dumpModuleGraph(const std::string filename);
    // memory profile
    const MemoryProfile &getMemoryProfile(void);
    const MemoryProfile &getMemoryProfile(const Program &p);
    void                printMemoryProfile(void) const;
    // garbage collection
    GarbageSchedule     makeGarbageSchedule(const Program &p,
                                            const bool pin = true) const;
    // high-water memory function
    Size                memoryNeeded(const Program &p);
    // genetic scheduler
    Program             schedule(const GeneticPar &par);
    // naive scheduler
    Program             naiveSchedule(void);
    // program validation (completeness, duplicates, dependency order)
    void                validateProgram(const Program &p) const;
    // split-phase boundary computation + validation: each phase's tagged
    // modules must form a contiguous ordered range of the schedule (phase
    // index non-decreasing over tagged modules in schedule order; untagged
    // modules may interleave). computePhaseRanges validates then returns the
    // ordered per-phase ranges; empty when split is unconfigured.
    void                validateSplitPhases(const Program &p) const;
    std::vector<PhaseRange> computePhaseRanges(const Program &p) const;
    // general execution
    void                executeProgram(const Program &p);
    void                executeProgram(const std::vector<std::string> &p);
    // generate result DB
    void                generateResultDb(void);
private:
    // environment shortcut
    DEFINE_ENV_ALIAS;
    // module graph
    void makeModuleGraph(void);
    void findModuleCycle(const Graph<unsigned int> &graph) const;
    // memory profile
    void                makeMemoryProfile(const Program &p);
    Program             defaultProfileProgram(void);
    [[noreturn]] void   rethrowWithContext(const unsigned int moduleAddress,
                                           const Exceptions::ObjectDefinition &exc) const;
    void                resetProfile(void);
    void                resizeProfile(void);
    void                updateProfile(const unsigned int address);
    void                memoryProfile(const unsigned int address);
    // database handling
    bool         hasDatabase(void) const;
    void         initDatabase(void);
    unsigned int dbInsertModuleType(const std::string type);
    unsigned int dbInsertObjectType(const std::string type, const std::string baseType);
private:
    // general
    std::string                         runId_;
    unsigned int                        traj_;
    // database
    Database                            *db_{nullptr};
    bool                                makeModuleDb_{true}, makeObjectDb_{true}, makeScheduleDb_{true};
    // module and related maps
    std::vector<ModuleInfo>             module_;
    std::map<std::string, unsigned int> moduleAddress_;
    int                                 currentModule_{-1};
    // module graph
    bool                                graphOutdated_{true};
    Graph<unsigned int>                 graph_;
    // memory profile
    bool                                memoryProfileOutdated_{true};
    MemoryProfile                       profile_;     
    // cache of which objects are Lattice fields (scatter-leaves), populated
    // during the memory-profile pre-pass when objects briefly exist. Used by
    // makeGarbageSchedule to mirror ensureShadowed's lattice-leaf behaviour so
    // the GC pins only lattice inputs actually scattered for the subgrid
    // rebuild (e.g. gauge links) and not, e.g., eigenvector packs reachable
    // only through a scattered lattice producer.
    std::map<unsigned int, bool>        objectIsLattice_;

    // per-phase addresses whose world-grid originals free at the phase's
    // FIRST step, once the pre-scatter has scattered their shadows (the
    // tagged modules read only the shadows from then on, so the world
    // original is dead). Filled by the pinned makeGarbageSchedule pass and
    // read at each phase entry for the free itself (freeSet) and its
    // report; mutable because
    // makeGarbageSchedule is const. Cleared at the top of every
    // makeGarbageSchedule call, so genetic-scheduler memoryNeeded probes
    // cannot leak stale content into executeProgram (which always rebuilds
    // its schedule first).
    mutable std::map<unsigned int, std::set<unsigned int>>
                                      earlyFreedByPhase_;

    // time profile
    GridTime                            totalTime_;
    std::map<std::string, GridTime>     moduleTimeProfile_, moduleTypeTimeProfile_;               
    // subgrid split state
    SubGrids                              subGrid_;            // this rank's subgrid
    int                                   me_{-1};             // this rank's subcomm index
    bool                                  splitConfigured_{false};
    std::set<unsigned int>                splitPhaseModules_;  // modules tagged into a subgrid (subgrid >= 0)
    // objects already scattered/rebuilt onto the subgrid during the current
    // split phase (avoids re-scattering; rank-independent, like the schedule).
    std::set<unsigned int>                shadowedObjects_;
    // per-object set of subgrid indices whose tagged modules (transitively)
    // consume the object -- computed by computeShadowConsumers before the
    // pre-scatter pass so store decisions see the COMPLETE consumer set.
    // Rank-independent. Empty/absent entry = store everywhere (per-step
    // late path, or nrhs == 1).
    std::map<unsigned int, std::set<int>> shadowConsumers_;
    // last step (within the active phase) whose tagged closure consumes
    // {object, subgrid} -- filled by computeShadowConsumers alongside
    // shadowConsumers_; drives the local per-step shadow frees (R2).
    // Rank-independent, like shadowConsumers_.
    std::map<std::pair<unsigned int, int>, unsigned int> shadowLastUse_;

// ============================================================================
// INSERT A: private state, after the shadowLastUse_ map declaration
// (VirtualMachine.hpp:297) — subgrid gather (reverse split) state
// ============================================================================

    // ---- subgrid gather (reverse split) state ----
    // per-phase gather set: raw <splitPhase> value -> ADDRESSES (ascending)
    // of that phase's tagged-produced objects with post-phase untagged
    // consumers. The phase-exit driver gathers them back to the world grid
    // right after the exit Barrier. Rank-independent (schedule-derived), so
    // every rank issues the identical collective sequence.
    std::map<unsigned int, std::vector<unsigned int>> gatherSchedule_;
    // producing phase (raw <splitPhase> value) of every gathered object;
    // consumed by the scatter-refinement (computeShadowConsumers /
    // ensureShadowed) so earlier-phase gathered twins are treated as world
    // objects
    std::map<unsigned int, unsigned int> producerPhaseOf_;
    // every address in any phase's gather set (fast membership for the
    // scatter-refinement)
    std::set<unsigned int> gatherAddresses_;
    // first world rank of each sub-comm (the gather shape/checkerboard
    // descriptor broadcast root, GatherCtx::bossRank), tabulated once in
    // buildSubGrids via one world allgather of me_
    std::vector<int> subCommFirstRank_;
    // once-per-program gather-side slot-delivery self-check flag (the mirror
    // of batchSlotDeliveryVerified_): before trusting any subToSlot
    // extraction, every rank proves unsplit slot p carries sub-comm
    // slotToSub[p]'s marker piece on this grid pair
    bool                                  gatherSlotDeliveryVerified_{false};    // lazily scatter (lattices) or rebuild (non-lattice objects) the global
    // inputs of a split-phase module onto the subgrid. Collective: the schedule,
    // splitPhase set and module inputs are identical on every rank, so every
    // rank performs the same scatter/rebuild at the same schedule step.
    void                                  ensureShadowed(const std::vector<unsigned int> &inputs);
    // phase 1 of consumer-gated shadow storage: fill shadowConsumers_ by
    // walking every tagged module's input closure with ensureShadowed's
    // classification rules (no side effects, no collectives).
    void                                  computeShadowConsumers(const Program &p,
                                                                const PhaseRange &range);
    // ---- batched sole-consumer scatter state (R4, Slice 6) ----
    // collection mode: when true, ensureShadowed defers sole-consumer
    // scatterable objects to batchQueue_ and rebuild-path objects to
    // rebuildQueue_ instead of scattering/rebuilding inline; the driver
    // then runs batches and deferred rebuilds (collect -> batch -> rebuild).
    // All four members are per-phase (cleared by the driver at its end —
    // which runs at each phase entry before any use) and rank-independent
    // (derived from the schedule + consumer sets), so batch group/slot
    // assignment is identical on every rank.
    bool                                  collectingBatches_{false};
    std::vector<unsigned int>             batchQueue_;
    std::set<unsigned int>                batchQueued_;
    std::vector<unsigned int>             rebuildQueue_;
    std::set<unsigned int>                rebuildQueued_;
    // one-time per-program self-check that gridSplitBatch really delivers
    // slot k to sub-comm k on this rank's grid pair (Phase 11): guards the
    // analytic permutation against layouts not hand-verified
    bool                                  batchSlotDeliveryVerified_{false};
    void                                  verifyBatchSlotDelivery(void);
    // three-pass pre-scatter (R4, D-pre-scatter-three-pass): collect
    // (defer sole-consumer scatterables + rebuilds), batch (one vector
    // split per type/shape group; slot k -> subgrid k), rebuild (deferred
    // producers in reverse collection order = valid topological order).
    void                                  executeBatchedPreScatter(
                                            const Program &p,
                                            const PhaseRange &range);

// ============================================================================
// INSERT B: private methods, after the executeBatchedPreScatter declaration
// (VirtualMachine.hpp:331-333) — first two declarations; the third
// (verifyGatherSlotDelivery) is authored in Phase 3 per the slice contract
// ============================================================================

    // schedule pre-pass (inside executeProgram, after computePhaseRanges,
    // before any step executes): fill gatherSchedule_/producerPhaseOf_/
    // gatherAddresses_ and hard-error the two illegal shapes — (1) an
    // untagged consumer INSIDE the producing phase's [first,last] span
    // (would need a mid-phase world collective: the forbidden deadlock
    // class), (2) a post-phase untagged consumer of a non-scatterable
    // tagged output. Scatterability comes from the objectIsLattice_ cache
    // (populated by the profiling pass, which always precedes execution).
    // Rank-independent.
    void                                  computeGatherSchedule(
                                            const Program &p,
                                            const std::vector<PhaseRange>
                                                &phaseRanges);
    // phase-exit driver: gather the exiting phase's set back to the world
    // grid — called from the phase-exit block right after the world
    // Barrier, before the per-phase state clears, with phaseIdx still
    // naming the EXITING phase. COLLECTIVE: every rank calls it in
    // lockstep; ascending-address order keeps the collective sequence
    // rank-uniform. Empty set: returns immediately, no collective.
    void                                  executePhaseGather(
                                            const unsigned int phase);
    // gather-side slot-delivery probe (the mirror of the Phase 11 scatter
    // probe): fill this rank's marker piece with its sub-comm index,
    // Grid_unsplit nrhs world-grid markers, and check every slot p carries
    // the constant slotToSub[p] — exactly what gridGatherLeaf's subToSlot
    // extraction relies on (gridSplitBatchSlotToSub is NOT an involution for
    // non-palindromic mixed radixes). Markers are vComplex-typed through
    // the same vComplex grid pair the production gather runs on (NOT the
    // scatter probe's vInteger — see the implementation's marker-type
    // rationale). COLLECTIVE: called by every rank in lockstep, once per
    // program, at the first non-empty phase gather.
    void                                  verifyGatherSlotDelivery(void);
};

/******************************************************************************
 *                   VirtualMachine template implementation                   *
 ******************************************************************************/
// module management ///////////////////////////////////////////////////////////
template <typename M>
void VirtualMachine::createModule(const std::string name)
{
    ModPt pt(new M(name));
    
    pushModule(pt);
}

template <typename M>
void VirtualMachine::createModule(const std::string name,
                               const typename M::Par &par)
{
    ModPt pt(new M(name));
    
    static_cast<M *>(pt.get())->setPar(par);
    pushModule(pt);
}

template <typename M>
M * VirtualMachine::getModule(const unsigned int address) const
{
    if (auto *pt = dynamic_cast<M *>(getModule(address)))
    {
        return pt;
    }
    else
    {
        HADRONS_ERROR(Definition, "module '" + module_[address].name
                     + "' does not have type " + typeid(M).name()
                     + "(has type: " + getModuleType(address) + ")");
    }
}

template <typename M>
M * VirtualMachine::getModule(const std::string name) const
{
    return getModule<M>(getModuleAddress(name));
}

END_HADRONS_NAMESPACE

#endif // Hadrons_VirtualMachine_hpp_
