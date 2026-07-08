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
    void                pushModule(ModPt &pt, const int subgrid = -1);
    template <typename M>
    void                createModule(const std::string name);
    template <typename M>
    void                createModule(const std::string name,
                                         const typename M::Par &par);
    void                createModule(const std::string name,
                                     const std::string type,
                                     XmlReader &reader,
                                     const std::string blockName = "options",
                                     const int subgrid = -1);
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
    void                printMemoryProfile(void) const;
    // garbage collection
    GarbageSchedule     makeGarbageSchedule(const Program &p) const;
    // high-water memory function
    Size                memoryNeeded(const Program &p);
    // genetic scheduler
    Program             schedule(const GeneticPar &par);
    // naive scheduler
    Program             naiveSchedule(void);
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
    // memory profile
    void makeMemoryProfile(void);
    void resetProfile(void);
    void resizeProfile(void);
    void updateProfile(const unsigned int address);
    void cleanEnvironment(void);
    void memoryProfile(const std::string name);
    void memoryProfile(const unsigned int address);
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
    // lazily scatter (lattices) or rebuild (non-lattice objects) the global
    // inputs of a split-phase module onto the subgrid. Collective: the schedule,
    // splitPhase set and module inputs are identical on every rank, so every
    // rank performs the same scatter/rebuild at the same schedule step.
    void                                  ensureShadowed(const std::vector<unsigned int> &inputs);
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
