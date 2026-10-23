/*
 * VirtualMachine.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
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

#include <Hadrons/VirtualMachine.hpp>
#include <functional>
#include <Hadrons/GeneticScheduler.hpp>
#include <Hadrons/StatLogger.hpp>
#include <Hadrons/ModuleFactory.hpp>

using namespace Grid;
 
using namespace Hadrons;

/******************************************************************************
 *                      VirtualMachine implementation                         *
 ******************************************************************************/
// trajectory counter //////////////////////////////////////////////////////////
void VirtualMachine::setTrajectory(const unsigned int traj)
{
    traj_ = traj;
}

unsigned int VirtualMachine::getTrajectory(void) const
{
    return traj_;
}

// run tag /////////////////////////////////////////////////////////////////////
void VirtualMachine::setRunId(const std::string id)
{
    runId_ = id;
}

std::string VirtualMachine::getRunId(void) const
{
    return runId_;
}

// database ////////////////////////////////////////////////////////////////////
void VirtualMachine::setDatabase(Database &db)
{
    db_ = &db;
    initDatabase();
}

void VirtualMachine::dbRestoreMemoryProfile(void)
{
    if (hasDatabase())
    {
        if (db_->tableExists("objects"))
        {
            auto         table = db_->getTable<ObjectEntry>("objects", "ORDER BY objectId");
            unsigned int nMod;

            auto comp  = [](const ObjectEntry &a, const ObjectEntry &b)
            {
                return (a.moduleId < b.moduleId);
            };

            if (env().getMaxAddress() > 0)
            {
                HADRONS_ERROR(Database, "environment is not empty");
            }
            if (table.size() == 0)
            {
                HADRONS_ERROR(Database, "object table is empty");
            }
            resetProfile();
            nMod = std::max_element(table.begin(), table.end(), comp)->moduleId + 1;
            profile_.module.resize(nMod);
            profile_.object.resize(table.size());
            for (auto &e: table)
            {
                profile_.object[e.objectId].module      = e.moduleId;
                profile_.object[e.objectId].size        = e.size;
                profile_.object[e.objectId].storage     = e.storageType;
                profile_.module[e.moduleId][e.objectId] = e.size;
                env().addObject(e.name, e.moduleId);
                assert(env().getObjectAddress(e.name) == e.objectId);
                env().setObjectStorage(e.objectId, e.storageType);
            }
            memoryProfileOutdated_ = false;
        }
    }
    else
    {
        HADRONS_ERROR(Database, "no database connected");
    }
}

void VirtualMachine::dbRestoreModules(void)
{
    if (hasDatabase())
    {
        if (db_->tableExists("modules"))
        {
            std::string prefix    = "HADRONS_NAMESPACE::";
            auto        modTable  = db_->getTable<ModuleEntry>("modules", "ORDER BY moduleId");
           
            if (getNModule() > 0)
            {
                HADRONS_ERROR(Database, "module graph is not empty");
            }
            if (modTable.size() == 0)
            {
                HADRONS_ERROR(Database, "module table is empty");
            }
            for (auto &e: modTable)
            {
                auto typeTable = db_->getTable<ModuleTypeEntry>("moduleTypes", 
                    "WHERE moduleTypeId = " + std::to_string(e.moduleTypeId));
                std::string type = typeTable.front().type;

                if ((type.size() > prefix.size()) 
                    and (type.substr(0, prefix.size()) == prefix))
                {
                    type = type.substr(prefix.size(), type.size() - prefix.size());
                }

                std::string par = "<top>" + e.parameters + "</top>";
                XmlReader   reader(par, true, "top");
                auto        &factory = ModuleFactory::getInstance();
                auto        pt       = factory.create(type, e.name);

                pt->parseParameters(reader, pt->parClassName());
                pushModule(pt);
            }
        }
    }
    else
    {
        HADRONS_ERROR(Database, "no database connected");
    }
}

VirtualMachine::Program VirtualMachine::dbRestoreSchedule(void)
{
    Program program;

    if (hasDatabase())
    {
        if (db_->tableExists("schedule"))
        {
            auto table = db_->getTable<ScheduleEntry>("schedule");

            if (table.size() == 0)
            {
                HADRONS_ERROR(Database, "schedule table is empty");
            }
            program.resize(table.size());
            for (auto &e: table)
            {
                program[e.step] = e.moduleId;
            }
        }
    }
    else
    {
        HADRONS_ERROR(Database, "no database connected");
    }

    return program;
}

bool VirtualMachine::hasDatabase(void) const
{
    return ((db_ != nullptr) and db_->isConnected());
}

void VirtualMachine::initDatabase(void)
{
    db_->execute("PRAGMA foreign_keys = ON;");
    if (!db_->tableExists("global"))
    {
        db_->createTable<GlobalEntry>("global");
    }
    if (!db_->tableExists("moduleTypes"))
    {
        db_->createTable<ModuleTypeEntry>("moduleTypes", "PRIMARY KEY(moduleTypeId)");
    }
    if (!db_->tableExists("modules"))
    {
        db_->createTable<ModuleEntry>("modules", "PRIMARY KEY(moduleId)"
            "FOREIGN KEY(moduleTypeId) REFERENCES moduleTypes(moduleTypeId)");
    }
    else if (!db_->tableEmpty("modules"))
    {
        LOG(Message) << "The module table in '" << db_->getFilename() << "' is not empty, it will not be altered" << std::endl;
        makeModuleDb_ = false;
    }
    if (!db_->tableExists("objectTypes"))
    {
        db_->createTable<ObjectTypeEntry>("objectTypes", "PRIMARY KEY(objectTypeId)");
    }
    if (!db_->tableExists("objects"))
    {
        db_->createTable<ObjectEntry>("objects", "PRIMARY KEY(objectId)," 
            "FOREIGN KEY(moduleId) REFERENCES modules(moduleId),"
            "FOREIGN KEY(objectTypeId) REFERENCES objectTypes(objectTypeId)");
    }
    else if (!db_->tableEmpty("objects"))
    {
        LOG(Message) << "The object table in '" << db_->getFilename() << "' is not empty, it will not be altered" << std::endl;
        makeObjectDb_ = false;
    }
    if (!db_->tableExists("schedule"))
    {
        db_->createTable<ScheduleEntry>("schedule", "PRIMARY KEY(step)," 
            "FOREIGN KEY(moduleId) REFERENCES modules(moduleId)");
    }
    else if (!db_->tableEmpty("schedule"))
    {
        LOG(Message) << "The schedule table in '" << db_->getFilename() << "' is not empty, it will not be altered" << std::endl;
        makeScheduleDb_ = false;
    }
    db_->execute(
        "CREATE VIEW IF NOT EXISTS vModules AS                                     "
        "SELECT moduleId,                                                          "
        "       modules.name,                                                      "
        "       moduleTypes.type AS type,                                          "
        "       modules.parameters                                                 "
        "FROM modules                                                              "
        "INNER JOIN moduleTypes ON modules.moduleTypeId = moduleTypes.moduleTypeId "
        "ORDER BY moduleId;                                                        "
    );
    db_->execute(
        "CREATE VIEW IF NOT EXISTS vObjects AS                                     "
        "SELECT objectId,                                                          "
        "       objects.name,                                                      "
        "       objectTypes.type AS type,                                          "
        "       objectTypes.baseType AS baseType,                                  "
        "       objects.size*1.0/1024/1024 AS sizeMB,                              "
        "       objects.storageType,                                               "
        "       modules.name AS module                                             "
        "FROM objects                                                              "
        "INNER JOIN objectTypes ON objects.objectTypeId = objectTypes.objectTypeId "
        "INNER JOIN modules     ON objects.moduleId = modules.moduleId             "
        "ORDER BY objectId;                                                        "
    );
    db_->execute(
        "CREATE VIEW IF NOT EXISTS vSchedule AS                                    "
        "SELECT step,                                                              "
        "       modules.name AS module                                             "
        "FROM schedule                                                             "
        "INNER JOIN modules ON schedule.moduleId = modules.moduleId                "
        "ORDER BY step;                                                            "
    );
}

unsigned int VirtualMachine::dbInsertModuleType(const std::string type)
{
    QueryResult r = db_->execute("SELECT moduleTypeId FROM moduleTypes "
                                 "WHERE type = '" + type + "';");

    if (r.rows() == 0)
    {
        ModuleTypeEntry e;

        r = db_->execute("SELECT COUNT(*) FROM moduleTypes;");
        e.moduleTypeId = std::stoi(r[0][0]);
        e.type         = type;
        db_->insert("moduleTypes", e);

        return e.moduleTypeId;
    }
    else
    {
        return std::stoi(r[0][0]);
    }
}

unsigned int VirtualMachine::dbInsertObjectType(const std::string type, 
                                             const std::string baseType)
{
    QueryResult r = db_->execute("SELECT objectTypeId FROM objectTypes "
                                 "WHERE type = '" + type + "' "
                                 "AND baseType = '" + baseType + "';");

    if (r.rows() == 0)
    {
        ObjectTypeEntry e;

        r = db_->execute("SELECT COUNT(*) FROM objectTypes;");
        e.objectTypeId = std::stoi(r[0][0]);
        e.type         = type;
        e.baseType     = baseType;
        db_->insert("objectTypes", e);

        return e.objectTypeId;
    }
    else
    {
        return std::stoi(r[0][0]);
    }
}

// module management ///////////////////////////////////////////////////////////
void VirtualMachine::pushModule(VirtualMachine::ModPt &pt, const int subgrid,
                                const int splitPhase)
{
    std::string name = pt->getName();
    
    if (!hasModule(name))
    {   
        // module registration -------------------------------------------------
        unsigned int address;
        ModuleInfo   mtmp;
        ModuleBase   *m;

        mtmp.data = std::move(pt);
        mtmp.type = typeIdPt(*mtmp.data.get());
        mtmp.name = name;
        module_.push_back(std::move(mtmp));
        address              = static_cast<unsigned int>(module_.size() - 1);
        moduleAddress_[name] = address;
        m                    = getModule(address);

        // input & output scan -------------------------------------------------
        ModuleInfo &mInfo = module_[address];
        mInfo.subgrid = subgrid;        // store subgrid tag (-1 = global)
        mInfo.splitPhase = splitPhase;  // store split-phase tag (0 = first phase)

        // scan inputs and add objects to the environment if necessary
        for (auto &in: m->getInput())
        {
            if (!env().hasObject(in))
            {
                // if object does not exist, add it with no creator module
                env().addObject(in);
                memoryProfileOutdated_ = true;
            }
            mInfo.input.push_back(env().getObjectAddress(in));
        }

        // scan outputs and add objects to the environment if necessary
        for (auto &out: m->getOutput())
        {
            if (!env().hasObject(out))
            {
                // output does not exists, add it
                env().addObject(out, address);
                memoryProfileOutdated_ = true;
            }
            else
            {
                if (env().getObjectModule(env().getObjectAddress(out)) < 0)
                {
                    // output exists but without creator, correct it
                    env().setObjectModule(env().getObjectAddress(out), address);
                }
                else if (env().getObjectModule(env().getObjectAddress(out)) != address)
                {
                    // output already produced by another module, error
                    HADRONS_ERROR_REF(ObjectDefinition, "object '" + out
                                 + "' is already produced by module '"
                                 + module_[env().getObjectModule(out)].name
                                 + "' (while pushing module '" + name + "')",
                                 env().getObjectAddress(out));
                }
            }
            mInfo.output.push_back(env().getObjectAddress(out));
        }

        // extra user-specified dependencies -----------------------------------
 #define VEC_HAS_ELEMENT(v, x) (std::find(v.begin(), v.end(), x) != v.end())

        for (auto &dep: m->getObjectDependencies())
        {
            unsigned int o = env().getObjectAddress(dep.first);
            unsigned int d = env().getObjectAddress(dep.second);

            if (!VEC_HAS_ELEMENT(mInfo.input, o)
                and !VEC_HAS_ELEMENT(mInfo.input, d)
                and !VEC_HAS_ELEMENT(mInfo.output, o)
                and !VEC_HAS_ELEMENT(mInfo.output, d))
            {
                HADRONS_ERROR(Definition, "Module '" + name + "' has a dependency with"
                              + " an object which is neither an input or an output"
                              + " (dependency between objects '" + dep.first
                              + "' and '" + dep.second + "')");
            }
            env().addObjectDependency(o, d);
        }

#undef VEC_HAS_ELEMENT

        // creating entry in database ------------------------------------------
        if (hasDatabase() and makeModuleDb_)
        {
            ModuleEntry e;

            e.moduleId     = address;
            e.name         = name;
            e.moduleTypeId = dbInsertModuleType(getModuleType(address));
            e.parameters   = getModule(address)->parString();
            db_->insert("modules", e);
        }
        graphOutdated_ = true;
    }
    else
    {
        HADRONS_ERROR(Definition, "module '" + name + "' already exists");
    }
}

unsigned int VirtualMachine::getNModule(void) const
{
    return module_.size();
}

void VirtualMachine::createModule(const std::string name, const std::string type,
                                  XmlReader &reader, const std::string blockName,
                                  const int subgrid, const int splitPhase)
{
    auto &factory = ModuleFactory::getInstance();
    auto pt       = factory.create(type, name);
    
    pt->parseParameters(reader, blockName);
    pushModule(pt, subgrid, splitPhase);
}

int VirtualMachine::getModuleSubgrid(const unsigned int address) const
{
    if (hasModule(address))
    {
        return module_[address].subgrid;
    }
    return -1;
}

int VirtualMachine::getModuleSplitPhase(const unsigned int address) const
{
    if (hasModule(address))
    {
        return module_[address].splitPhase;
    }
    return 0;
}

void VirtualMachine::buildSubGrids(const std::vector<int> &mpiSplit)
{
    if (mpiSplit.empty())
    {
        splitConfigured_ = false;
        return;
    }
    GridCartesian *world = env().getGrid();          // default 4d world grid
    int            childsize = 1;
    for (auto &p: mpiSplit) childsize *= p;
    int            nproc = world->_Nprocessors;
    int            nrhs  = nproc / childsize;
    if (childsize * nrhs != nproc)
    {
        HADRONS_ERROR(Definition, "mpiSplit product " + std::to_string(childsize)
                      + " does not divide world size " + std::to_string(nproc));
    }
    for (unsigned int d = 0; d < mpiSplit.size(); ++d)
    {
        if (world->_processors[d] % mpiSplit[d] != 0)
        {
            HADRONS_ERROR(Definition, "mpiSplit[" + std::to_string(d) + "]="
                          + std::to_string(mpiSplit[d])
                          + " does not divide world procs["
                          + std::to_string(d) + "]="
                          + std::to_string(world->_processors[d]));
        }
    }
    Coordinate dims  = world->_gdimensions;
    Coordinate simd  = world->_simd_layout;
    Coordinate procSplit(mpiSplit);   // AcceleratorVector<int> from std::vector<int>
    int        splitRank;
    subGrid_.grid.reset(new GridCartesian(dims, simd, procSplit, *world, splitRank));
    subGrid_.mpiSplit = procSplit;
    subGrid_.nrhs     = nrhs;
    subGrid_.me       = splitRank;
    me_               = splitRank;
    splitConfigured_  = true;
    computeSplitPhase();
    LOG(Message) << "Split grid: " << nrhs << " sub-comms, this rank is subcomm "
                 << me_ << std::endl;
}

void VirtualMachine::computeSplitPhase(void)
{
    // The split phase is exactly the set of modules explicitly TAGGED into a
    // subgrid (subgrid >= 0). Untagged producers (solvers / fermion actions)
    // are deliberately NOT pulled in by a transitive closure: they execute once
    // on the world grid at their own schedule step and are rebuilt on the
    // subgrid lazily by ensureShadowed when a tagged consumer first needs them
    // (see ensureShadowed and makeGarbageSchedule). Pulling them in caused them
    // to execute under the subgrid scope at their original step, producing a
    // grid-mismatched object and skipping the subgrid rebuild.
    // This cached tagged set is reused by makeGarbageSchedule to pin the lattice
    // inputs (e.g. gauge links) needed for those lazy subgrid rebuilds.
    splitPhaseModules_.clear();
    for (unsigned int m = 0; m < module_.size(); ++m)
    {
        if (module_[m].subgrid >= 0)
        {
            splitPhaseModules_.insert(m);
        }
    }
}

ModuleBase * VirtualMachine::getModule(const unsigned int address) const
{
    if (hasModule(address))
    {
        return module_[address].data.get();
    }
    else
    {
        HADRONS_ERROR(Definition, "no module with address " + std::to_string(address));
    }
}

ModuleBase * VirtualMachine::getModule(const std::string name) const
{
    return getModule(getModuleAddress(name));
}

unsigned int VirtualMachine::getModuleAddress(const std::string name) const
{
    if (hasModule(name))
    {
        return moduleAddress_.at(name);
    }
    else
    {
        HADRONS_ERROR(Definition, "no module with name '" + name + "'");
    }
}

std::string VirtualMachine::getModuleName(const unsigned int address) const
{
    if (hasModule(address))
    {
        return module_[address].name;
    }
    else
    {
        HADRONS_ERROR(Definition, "no module with address " + std::to_string(address));
    }
}

std::string VirtualMachine::getModuleType(const unsigned int address) const
{
    if (hasModule(address))
    {
        return typeName(module_[address].type);
    }
    else
    {
        HADRONS_ERROR(Definition, "no module with address " + std::to_string(address));
    }
}

std::string VirtualMachine::getModuleType(const std::string name) const
{
    return getModuleType(getModuleAddress(name));
}

std::string VirtualMachine::getModuleNamespace(const unsigned int address) const
{
    std::string type = getModuleType(address), ns;
    
    auto pos2 = type.rfind("::");
    auto pos1 = type.rfind("::", pos2 - 2);
    
    return type.substr(pos1 + 2, pos2 - pos1 - 2);
}

std::string VirtualMachine::getModuleNamespace(const std::string name) const
{
    return getModuleNamespace(getModuleAddress(name));
}

int VirtualMachine::getCurrentModule(void) const
{
    return currentModule_;
}

bool VirtualMachine::hasModule(const unsigned int address) const
{
    return (address < module_.size());
}

bool VirtualMachine::hasModule(const std::string name) const
{
    return (moduleAddress_.find(name) != moduleAddress_.end());
}

// print VM content ////////////////////////////////////////////////////////////
void VirtualMachine::printContent(void) const
{
    LOG(Debug) << "Modules: " << std::endl;
    for (unsigned int i = 0; i < module_.size(); ++i)
    {
        LOG(Debug) << std::setw(4) << i << ": "
                   << getModuleName(i) << std::endl;
    }
}

// module graph ////////////////////////////////////////////////////////////////
Graph<unsigned int> VirtualMachine::getModuleGraph(void)
{
    if (graphOutdated_)
    {
        makeModuleGraph();
        graphOutdated_ = false;
    }

    return graph_;
}

void VirtualMachine::makeModuleGraph(void)
{
    Graph<unsigned int> graph;
    
    // create vertices
    for (unsigned int m = 0; m < module_.size(); ++m)
    {
        graph.addVertex(m);
    }
    // create edges
    for (unsigned int m = 0; m < module_.size(); ++m)
    {
        for (auto &in: module_[m].input)
        {
            int min = env().getObjectModule(in);

            if (min < 0)
            {
                HADRONS_ERROR_REF(ObjectDefinition, "dependency '" 
                             + env().getObjectName(in) + "' (address " 
                             + std::to_string(in)
                             + ", of module '" + getModuleName(m) 
                             + "') is not produced by any module", in);
            }
            else
            {
                graph.addEdge(min, m);
            }
        }
    }
    // name the cycle (if any) before assignment — Graph::topoSort's generic
    // error carries no module names
    findModuleCycle(graph);
    graph_ = graph;
}

void VirtualMachine::findModuleCycle(const Graph<unsigned int> &graph) const
{
    // three-colour DFS with the current path: a back-edge to a grey vertex
    // closes a cycle, reported as a module-name chain
    std::map<unsigned int, bool> grey, done;
    std::vector<unsigned int>    path;

    std::function<void(const unsigned int)> visit = [&](const unsigned int v)
    {
        grey[v] = true;
        path.push_back(v);
        for (auto &c: graph.getChildren(v))
        {
            if (grey.count(c) and grey[c])
            {
                auto it      = std::find(path.begin(), path.end(), c);
                auto cycle   = std::vector<unsigned int>(it, path.end());
                std::string msg;

                for (auto &m: cycle)
                {
                    msg += "module '" + getModuleName(m) + "' -> ";
                }
                msg += "module '" + getModuleName(c) + "'";
                HADRONS_ERROR(Definition, "module graph has a dependency"
                              " cycle: " + msg);
            }
            if (!done.count(c) or !done[c])
            {
                visit(c);
            }
        }
        grey[v] = false;
        done[v] = true;
        path.pop_back();
    };

    for (auto &v: graph.getVertices())
    {
        if (!done.count(v) or !done[v])
        {
            visit(v);
        }
    }
}

// dump GraphViz graph /////////////////////////////////////////////////////////
void VirtualMachine::dumpModuleGraph(std::ostream &out)
{
    makeModuleGraph();
    out << "digraph hadrons {" << std::endl;
    out << "node [shape=record, fontname=\"Courier\", fontsize=\"11\"];" << std::endl;
    out << "graph [fontname = \"Courier\", fontsize=\"11\"];" << std::endl;
    out << "edge [fontname = \"Courier\", fontsize=\"11\"];"<< std::endl;
    for (unsigned int m = 0; m < module_.size(); ++m)
    {
        for (auto &in: module_[m].input)
        {
            int min = env().getObjectModule(in);

            out << min << " -> " << m << " [ label = \""
                << env().getObjectName(in) << "\" ];" << std::endl;
        }
    }
    for (unsigned int m = 0; m < module_.size(); ++m)
    {
        out <<  m << " [ label = \"{<f0> " << getModule(m)->getRegisteredName()
            << " |<f1> " << getModuleName(m) << "}\" ];" << std::endl;
    }
    out << "}\n" << std::endl;
}

void VirtualMachine::dumpModuleGraph(void)
{
    dumpModuleGraph(std::cout);
}

void VirtualMachine::dumpModuleGraph(const std::string filename)
{
    std::ofstream f(filename);

    dumpModuleGraph(f);
}

// memory profile //////////////////////////////////////////////////////////////
const VirtualMachine::MemoryProfile & VirtualMachine::getMemoryProfile(void)
{
    if (memoryProfileOutdated_)
    {
        makeMemoryProfile(defaultProfileProgram());
        memoryProfileOutdated_ = false;
    }

    return profile_;
}

const VirtualMachine::MemoryProfile & VirtualMachine::getMemoryProfile(const Program &p)
{
    if (memoryProfileOutdated_)
    {
        makeMemoryProfile(p);
        memoryProfileOutdated_ = false;
    }

    return profile_;
}

void VirtualMachine::makeMemoryProfile(const Program &p)
{
    bool protect = env().objectsProtected();
    bool hmsg    = HadronsLogMessage.isActive();
    bool gmsg    = GridLogMessage.isActive();

    resetProfile();
    profile_.module.resize(getNModule());
    env().protectObjects(false);
    GridLogMessage.Active(false);
    HadronsLogMessage.Active(false);
    // Single forward pass over the program, mirroring executeProgram's loop:
    // each module's setup() runs exactly once, in a valid topological order,
    // and objects are freed exactly as the real run would free them. The
    // garbage schedule for the pass is built WITHOUT split-grid pinning:
    // pinning only protects lattice inputs from early freeing while
    // ensureShadowed scatters them at execution time, and a setup()-only
    // world-grid pass never scatters. With the objectIsLattice_ cache still
    // empty at this point, pinning would over-pin (walk through lattice
    // scatter-leaves into unrelated heavy producers). executeProgram rebuilds
    // its own schedule, pinned, after this pass fills the cache.
    GarbageSchedule freeProg;
    try
    {
        freeProg = makeGarbageSchedule(p, /*pin=*/false);
        for (unsigned int i = 0; i < p.size(); ++i)
        {
            auto a = p[i];

            LOG(Debug) << "Profiling memory for module '" << module_[a].name
                       << "' (" << a << ")" << std::endl;
            memoryProfile(a);
            env().freeSet(freeProg[i]);
            // mirror executeProgram's temporary sweep: envTmp objects created by
            // this module's setup() live only for this step
            for (unsigned int o = 0; o < env().getMaxAddress(); ++o)
            {
                if (env().getObjectStorage(o) == Environment::Storage::temporary)
                {
                    env().freeObject(o);
                }
            }
        }
        // objects invisible to the up-front schedule (undeclared outputs, caches)
        // survive the pass and are released here
        env().freeAll();
    }
    catch (...)
    {
        // restore log/protect state before the error propagates: every throw
        // path now ends in Exceptions::abort, whose output must not be
        // silenced (Step-5 triage; design follow-up:
        // .rpiv/artifacts/designs/2026-10-15_18-30-55_single-pass-memory-profiler.md)
        env().protectObjects(protect);
        GridLogMessage.Active(gmsg);
        HadronsLogMessage.Active(hmsg);
        throw;
    }
    env().protectObjects(protect);
    GridLogMessage.Active(gmsg);
    HadronsLogMessage.Active(hmsg);
    if (hasDatabase() and makeObjectDb_)
    {
        for (unsigned int i = 0; i < profile_.object.size(); ++i)
        {
            ObjectEntry o;

            o.objectId     = i;
            o.name         = env().getObjectName(i);
            o.objectTypeId = dbInsertObjectType(env().getObjectDerivedType(i),
                                                env().getObjectType(i));
            o.size         = profile_.object[i].size;
            o.moduleId     = profile_.object[i].module;
            o.storageType  = profile_.object[i].storage;
            db_->insert("objects", o);
        }
    }
}

VirtualMachine::Program VirtualMachine::defaultProfileProgram(void)
{
    Program p;

    for (unsigned int m = 0; m < getNModule(); ++m)
    {
        p.push_back(m);
    }
    // declaration order is a valid profiling seed unless an input's producer
    // is declared after its consumer (the naive scheduler's condition); fall
    // back to the module graph's deterministic topological order otherwise
    for (unsigned int m = 0; m < getNModule(); ++m)
    {
        for (auto &in: module_[m].input)
        {
            if (env().getObjectModule(in) > static_cast<int>(m))
            {
                return getModuleGraph().topoSort();
            }
        }
    }

    return p;
}

void VirtualMachine::rethrowWithContext(const unsigned int moduleAddress,
                                        const Exceptions::ObjectDefinition &exc) const
{
    std::string msg = exc.what();

    // keep the original message at the front (substring contracts, e.g.
    // test/check_setup_error.py) and append dependency context: calling
    // module, its outputs, failing object, that object's producer
    msg += " [in module '" + getModuleName(moduleAddress)
         + "' (address " + std::to_string(moduleAddress);
    if (!module_[moduleAddress].output.empty())
    {
        msg += ", building:";
        for (auto &o: module_[moduleAddress].output)
        {
            msg += " '" + env().getObjectName(o) + "'";
        }
    }
    msg += ")";
    auto address = exc.getAddress();
    if (address < env().getMaxAddress())
    {
        int producer = env().getObjectModule(address);

        msg += "; object '" + env().getObjectName(address)
             + "' (address " + std::to_string(address)
             + ") is produced by ";
        if (producer >= 0)
        {
            msg += "module '" + getModuleName(producer)
                 + "' (address " + std::to_string(producer) + ")";
        }
        else
        {
            msg += "no module";
        }
    }
    msg += "]";
    HADRONS_ERROR_REF(ObjectDefinition, msg, address);
}

void VirtualMachine::printMemoryProfile(void) const
{
    LOG(Debug) << "Memory profile:" << std::endl;
    LOG(Debug) << "----------------" << std::endl;
    for (unsigned int a = 0; a < profile_.module.size(); ++a)
    {
        LOG(Debug) << getModuleName(a) << " (" << a << ")" << std::endl;
        for (auto &o: profile_.module[a])
        {
            LOG(Debug) << "|__ " << env().getObjectName(o.first) << " ("
                       << profile_.object[o.first].storage << " "
                       << sizeString(o.second) << ")" << std::endl;
        }
        LOG(Debug) << std::endl;
    }
    LOG(Debug) << "----------------" << std::endl;
}

void VirtualMachine::resetProfile(void)
{
    profile_.module.clear();
    profile_.object.clear();
    objectIsLattice_.clear();
}

void VirtualMachine::resizeProfile(void)
{
    if (env().getMaxAddress() > profile_.object.size())
    {
        MemoryPrint empty;

        empty.size   = 0;
        empty.module = -1;
        profile_.object.resize(env().getMaxAddress(), empty);
    }
}

void VirtualMachine::updateProfile(const unsigned int address)
{
    resizeProfile();
    for (unsigned int a = 0; a < env().getMaxAddress(); ++a)
    {
        int envMod = env().getObjectModule(a);

        if ((env().hasCreatedObject(a) or (envMod == address)) and (profile_.object[a].module == -1))
        {
            profile_.object[a].module   = address;
            profile_.object[a].size     = env().getObjectSize(a);
            profile_.object[a].storage  = env().getObjectStorage(a);
            profile_.module[address][a] = profile_.object[a].size;
            if (envMod < 0)
            {
                env().setObjectModule(a, address);
            }
            // record lattice-ness while the object exists (it is freed again
            // by the profiling pass), for makeGarbageSchedule's GC pinning.
            if (env().hasCreatedObject(a))
            {
                objectIsLattice_[a] = env().isLatticeObject(a);
            }
        }
    }
}

void VirtualMachine::memoryProfile(const unsigned int address)
{
    auto m = getModule(address);

    LOG(Debug) << "Setting up module '" << m->getName() 
               << "' (" << address << ")" << std::endl;
    currentModule_ = address;
    try
    {
        m->setup();
    }
    catch (Exceptions::ObjectDefinition &exc)
    {
        // In the forward pass every producer in p runs before its consumers,
        // so an uncreated dependency can only mean the module uses memory it
        // did not declare through getInput()/getObjectDependencies(): a
        // declaration bug, reported with full context, no recursive retry.
        rethrowWithContext(address, exc);
    }
    currentModule_ = -1;
    updateProfile(address);
}

// garbage collector ///////////////////////////////////////////////////////////
VirtualMachine::GarbageSchedule 
VirtualMachine::makeGarbageSchedule(const Program &p, const bool pin) const
{
    GarbageSchedule freeProg;
    
    freeProg.resize(p.size());

    // earliest time to destroy object ignoring dependencies
    std::function<unsigned int(const unsigned int)> earliestTimeNoDep = 
    [&](const unsigned int a)
    {
        
        auto pred = [a, this](const unsigned int b)
        {
            auto &in = module_[b].input;
            auto it  = std::find(in.begin(), in.end(), a);
            
            return (it != in.end()) or (b == env().getObjectModule(a));
        };
        auto it = std::find_if(p.rbegin(), p.rend(), pred);
        if (it == p.rend())
        {
            auto owner = env().getObjectModule(a);

            HADRONS_ERROR(Logic, "object '" + env().getObjectName(a)
                          + "' (address " + std::to_string(a) + ", owned by "
                          + ((owner >= 0)
                             ? ("module '" + getModuleName(owner) + "'")
                             : std::string("no module"))
                          + ") has neither a consumer nor an owner in the program");
        }

        return std::distance(it, p.rend()) - 1;
    };

    // earliest time to destroy object (taking dependencies into account)
    std::function<unsigned int(const unsigned int)> earliestTime = 
    [&](const unsigned int a)
    {
        unsigned int t = 0;

        t = std::max(t, earliestTimeNoDep(a));
        for (auto &d: env().getObjectDependencies(a))
        {
            t = std::max(t, earliestTime(d));
        }

        return t;
    };

    // ---- split-grid: keep lattice inputs alive for the subgrid rebuild --------
    // Tagged modules run on a sub-comm. Their non-lattice inputs (solvers /
    // fermion actions) are rebuilt on the subgrid by ensureShadowed, which
    // recursively Grid_split-scatters the producers' lattice inputs (e.g.
    // gauge links) from the GLOBAL store. Those lattice inputs would
    // otherwise be freed by normal GC right after their last standard
    // consumer — a use-after-free. Pin every object in the transitive input
    // closure of the tagged modules (stopping at tagged producers, which
    // self-rebuild, and lattice scatter-leaves) to the PRE-SCATTER STEP OF
    // EVERY PHASE WHOSE CLOSURE NEEDS IT (the latest such phase): each
    // phase's pre-scatter re-scatters from the WORLD original (earlier
    // phases' shadow copies are keyed to their own phase and invisible), so
    // the original must survive every needing phase's pre-scatter. After
    // that scatter a subgrid copy lives in the shadow store, so the global
    // original is dead and may be freed at the end of that step (R3: the
    // pin ends at the pre-scatter, not at the first tagged consumer — this
    // removes the double-hold window; untagged later readers keep normal GC
    // via the earliestTime max below).
    std::map<unsigned int, unsigned int> shadowFirstNeed;
    if (pin && splitConfigured_ && !splitPhaseModules_.empty())
    {
        const auto    phaseRanges = computePhaseRanges(p);
        unsigned int  phaseIdx    = 0;
        for (unsigned int i = 0; i < p.size(); ++i)
        {
            if (!splitPhaseModules_.count(p[i]))
            {
                continue;
            }
            while (phaseIdx + 1 < phaseRanges.size()
                   and i > phaseRanges[phaseIdx].last)
            {
                phaseIdx++;
            }
            const unsigned int pinStep = phaseRanges[phaseIdx].first;
            // BFS this tagged module's shadow closure; pin every object to
            // the latest phase pre-scatter needing it. The visited set is
            // per tagged module (unlike the old global shadowFirstNeed.count
            // dedup) so later phases can RAISE an object's pin to their own
            // pre-scatter step.
            std::set<unsigned int> visited;
            std::vector<unsigned int> stack(module_[p[i]].input);
            while (!stack.empty())
            {
                unsigned int a = stack.back();
                stack.pop_back();
                if (visited.count(a))
                {
                    continue;
                }
                visited.insert(a);
                auto it = shadowFirstNeed.find(a);
                if (it == shadowFirstNeed.end())
                {
                    shadowFirstNeed[a] = pinStep;
                }
                else
                {
                    it->second = std::max(it->second, pinStep);
                }
                int m = env().getObjectModule(a);
                if (m < 0)
                {
                    continue;   // externally provided object
                }
                if (module_[static_cast<unsigned int>(m)].subgrid >= 0)
                {
                    continue;   // tagged producer: rebuilt on its subgrid
                }
                // Lattice objects are scatter-leaves in ensureShadowed: they
                // are Grid_split onto the subgrid but their producer is NOT
                // rebuilt, so do NOT descend into the producer's inputs here.
                // Without this, the closure would walk through a scattered
                // lattice (e.g. a guess propagator) into unrelated heavy
                // producers (e.g. an eigenvector pack) and pin them, defeating
                // the split's memory reclamation. Objects not seen during
                // the profiling pass are descended to avoid under-pinning a
                // genuine rebuild input.
                auto lat = objectIsLattice_.find(a);
                if (lat != objectIsLattice_.end() && lat->second)
                {
                    continue;   // lattice scatter-leaf: do not descend
                }
                for (auto &in: module_[static_cast<unsigned int>(m)].input)
                {
                    stack.push_back(in);
                }
            }
        }
    }

    for (unsigned int a = 0; a < env().getMaxAddress(); ++a)
    {
        if (env().getObjectStorage(a) == Environment::Storage::standard)
        {
            unsigned int t = earliestTime(a);
            auto         it = shadowFirstNeed.find(a);

            if (it != shadowFirstNeed.end())
            {
                t = std::max(t, it->second);
            }
            freeProg[t].insert(a);
        }
    }

    return freeProg;
}

// high-water memory function //////////////////////////////////////////////////
VirtualMachine::Size VirtualMachine::memoryNeeded(const Program &p)
{
    const MemoryProfile &profile = getMemoryProfile();
    GarbageSchedule     freep    = makeGarbageSchedule(p);
    Size                current = 0, max = 0;

    for (unsigned int i = 0; i < p.size(); ++i)
    {
        for (auto &o: profile.module[p[i]])
        {
            current += o.second;
        }
        max = std::max(current, max);
        for (auto &o: freep[i])
        {
            current -= profile.object[o].size;
        }
        // temporaries live only during their own module's step (the
        // execution loop sweeps them at the end of each step), so they do
        // not contribute to the running total of later steps
        for (auto &o: profile.module[p[i]])
        {
            if (profile.object[o.first].storage == Environment::Storage::temporary)
            {
                current -= o.second;
            }
        }
    }

    return max;
}

// genetic scheduler ///////////////////////////////////////////////////////////
VirtualMachine::Program VirtualMachine::schedule(const GeneticPar &par)
{
    typedef GeneticScheduler<Size, unsigned int> Scheduler;

    auto graph = getModuleGraph();

    //constrained topological sort using a genetic algorithm
    LOG(Message) << "Scheduling computation..." << std::endl;
    LOG(Message) << "               #module= " << graph.size() << std::endl;
    LOG(Message) << "       population size= " << par.popSize << std::endl;
    LOG(Message) << "       max. generation= " << par.maxGen << std::endl;
    LOG(Message) << "  max. cst. generation= " << par.maxCstGen << std::endl;
    LOG(Message) << "         mutation rate= " << par.mutationRate << std::endl;
    // Build the memory profile once, up front, with a deterministic program
    // (declaration order, topoSort fallback — see defaultProfileProgram), so
    // the profiling pass never walks a random candidate from initPopulation
    // and the first fitness evaluations only consult the cached profile.
    getMemoryProfile();

    unsigned int          gen, prevPeak, nCstPeak = 0;
    std::random_device    rd;
    Scheduler::Parameters gpar;
    
    gpar.popSize      = par.popSize;
    gpar.mutationRate = par.mutationRate;
    gpar.seed         = rd();
    CartesianCommunicator::BroadcastWorld(0, &(gpar.seed), sizeof(gpar.seed));
    Scheduler::ObjFunc memPeak = [this](const Program &p)->Size
    {
        return memoryNeeded(p);
    };
    Scheduler scheduler(graph, memPeak, gpar);
    gen = 0;
    scheduler.initPopulation();
    LOG(Message) << "Start: " << sizeString(scheduler.getMinValue()) 
                 << std::endl;
    do
    {
        scheduler.nextGeneration();
        if (gen != 0)
        {
            if (prevPeak == scheduler.getMinValue())
            {
                nCstPeak++;
            }
            else
            {
                nCstPeak = 0;
            }
        }
        
        prevPeak = scheduler.getMinValue();
        if (gen % 10 == 0)
        {
            LOG(Message) << "Generation " << gen << ": "
                         << sizeString(scheduler.getMinValue()) << std::endl;
        }
        
        gen++;
    } while ((gen < par.maxGen) and (nCstPeak < par.maxCstGen));
    if (hasDatabase() and makeScheduleDb_)
    {
        Program p = scheduler.getMinSchedule();

        for (unsigned int i = 0; i < p.size(); ++i)
        {
            ScheduleEntry s;

            s.step     = i;
            s.moduleId = p[i];
            db_->insert("schedule", s);
        }
    }
    
    return scheduler.getMinSchedule();
}

// naive scheduler ///////////////////////////////////////////////////////////
VirtualMachine::Program VirtualMachine::naiveSchedule(void)
{
    LOG(Message) << "Using naive scheduler." << std::endl;
    auto graph = getModuleGraph();

    Program p;

    for (unsigned int i = 0; i < graph.size(); ++i)
    {
        p.push_back(i);
    }
    // declaration order must be a valid dependency order; validateProgram
    // reports the dependency, consumer and producer with correct names
    // (the former inline check passed a module index where an object
    // address was expected)
    validateProgram(p);

    if (hasDatabase() and makeScheduleDb_)
    {
        for (unsigned int i = 0; i < p.size(); ++i)
        {
            ScheduleEntry s;

            s.step     = i;
            s.moduleId = p[i];
            db_->insert("schedule", s);
        }
    }

    return p;
}

// program validation //////////////////////////////////////////////////////////
void VirtualMachine::validateProgram(const Program &p) const
{
    std::vector<int> position(getNModule(), -1);

    // every step must reference a valid module, exactly once
    for (unsigned int i = 0; i < p.size(); ++i)
    {
        if (p[i] >= getNModule())
        {
            HADRONS_ERROR(Range, "program step " + std::to_string(i)
                          + " refers to module address "
                          + std::to_string(p[i])
                          + ", but only " + std::to_string(getNModule())
                          + " modules exist");
        }
        if (position[p[i]] >= 0)
        {
            HADRONS_ERROR(Definition, "module '" + getModuleName(p[i])
                          + "' appears more than once in the program (steps "
                          + std::to_string(position[p[i]] + 1) + " and "
                          + std::to_string(i + 1) + ")");
        }
        position[p[i]] = static_cast<int>(i);
    }
    // every registered module must be scheduled
    for (unsigned int m = 0; m < getNModule(); ++m)
    {
        if (position[m] < 0)
        {
            HADRONS_ERROR(Definition, "module '" + getModuleName(m)
                          + "' is missing from the program");
        }
    }
    // every input's producer must be scheduled before its consumer
    for (unsigned int i = 0; i < p.size(); ++i)
    {
        for (auto &in: module_[p[i]].input)
        {
            int min = env().getObjectModule(in);

            // deliberately similar to (but independent of) makeModuleGraph's
            // unowned-input error: this copy covers DB-restored and loaded
            // programs, which never re-run makeModuleGraph (Step-5 triage note)
            if (min < 0)
            {
                HADRONS_ERROR_REF(ObjectDefinition, "dependency '"
                              + env().getObjectName(in) + "' (address "
                              + std::to_string(in) + ", of module '"
                              + getModuleName(p[i]) + "') is not produced by"
                              + " any module", in);
            }
        }
    }
}

// split-phase boundaries ///////////////////////////////////////////////////////
void VirtualMachine::validateSplitPhases(const Program &p) const
{
    if (!splitConfigured_)
    {
        return;   // no split configured: <subgrid>/<splitPhase> tags are inert
    }
    int          prevPhase = -1;
    unsigned int prevStep  = 0;
    for (unsigned int i = 0; i < p.size(); ++i)
    {
        if (module_[p[i]].subgrid < 0)
        {
            continue;
        }
        const int phase = module_[p[i]].splitPhase;
        if (phase < prevPhase)
        {
            HADRONS_ERROR(Definition, "split phase " + std::to_string(phase)
                          + " of module '" + getModuleName(p[i]) + "' (step "
                          + std::to_string(i + 1) + ") is lower than split"
                          + " phase " + std::to_string(prevPhase)
                          + " of module '" + getModuleName(p[prevStep])
                          + "' (step " + std::to_string(prevStep + 1) + "):"
                          + " each phase's tagged modules must form a"
                          + " contiguous ordered range of the schedule"
                          + " (untagged modules may interleave)");
        }
        prevPhase = phase;
        prevStep  = i;
    }
}

std::vector<VirtualMachine::PhaseRange>
VirtualMachine::computePhaseRanges(const Program &p) const
{
    std::vector<PhaseRange> ranges;

    if (!splitConfigured_)
    {
        return ranges;   // no split configured: one implicit phase (legacy)
    }
    validateSplitPhases(p);
    for (unsigned int i = 0; i < p.size(); ++i)
    {
        if (module_[p[i]].subgrid < 0)
        {
            continue;
        }
        const auto phase = static_cast<unsigned int>(module_[p[i]].splitPhase);
        if (!ranges.empty() and ranges.back().phase == phase)
        {
            ranges.back().last = i;   // same phase extends; untagged interleave
        }
        else
        {
            ranges.push_back(PhaseRange{phase, i, i});
        }
    }
    return ranges;
}

// general execution ///////////////////////////////////////////////////////////
#define BIG_SEP   "================"
#define SEP       "----------------"
#define SMALL_SEP "................"

// subgrid shadowing of a split-phase module's global inputs ////////////////////
// For each input produced by a global (non-split-phase) module and not yet
// shadowed:
//   - lattice inputs are scattered onto the subgrid via Grid_split;
//   - non-lattice inputs (e.g. fermion actions) cannot be Grid_split, so the
//     producer module's setup() is re-run under the subgrid scope in
//     shadow-create mode, building a fresh subgrid copy into the shadow store
//     (the world-grid original is preserved for global consumers).
// Recursion guarantees the producer's own inputs are shadowed first. This is
// collective: the schedule, splitPhase set and module inputs are identical on
// every rank, so all ranks issue the same Grid_split / setup() at the same step.
// consumer-subgrid sets for gated shadow storage ///////////////////////////
// Per split phase (computeShadowConsumers(p, range) at phase entry): for
// every object transitively needed by a TAGGED module of the phase, record
// the set of subgrids whose tagged modules need it. Mirrors ensureShadowed's
// traversal rules exactly (tagged-producer skip, no-module skip,
// created-ness check, scatter-leaf non-descent, clone non-descent, rebuild
// descent) so the sets cannot diverge from the stores phase 2 performs. A
// Grid_split yields every rank its OWN sub-comm's piece of the full-lattice
// object, so "consumed by subgrid s" means only ranks with me_ == s keep
// their piece; empty set (never pre-computed) means store everywhere.
// The created-ness check is the pre-scatter contract (Q7): every
// untagged-produced closure object must exist at its phase's pre-scatter —
// the per-step late path is the memory trap being removed.
void VirtualMachine::computeShadowConsumers(const Program &p,
                                            const PhaseRange &range)
{
    shadowConsumers_.clear();
    shadowLastUse_.clear();
    for (unsigned int i = range.first; i <= range.last; ++i)
    {
        if (module_[p[i]].subgrid < 0)
        {
            continue;
        }
        const int s = module_[p[i]].subgrid;
        // per-(object, subgrid) visited set: record once, descend once
        std::set<std::pair<unsigned int, int>> visited;
        std::vector<unsigned int> stack(module_[p[i]].input);
        while (!stack.empty())
        {
            unsigned int a = stack.back();
            stack.pop_back();
            if (visited.count({a, s}))
            {
                continue;
            }
            visited.insert({a, s});
            int m = env().getObjectModule(a);
            if (m < 0)
            {
                continue;   // externally provided object
            }
            if (module_[static_cast<unsigned int>(m)].subgrid >= 0)
            {
                continue;   // tagged producer: self-rebuilds on its subgrid
            }
            if (!env().hasCreatedObject(a))
            {
                // pre-scatter contract (Q7): the per-step late shadow path is
                // the memory trap being removed. Every untagged-produced
                // closure object must exist at its phase's pre-scatter; the
                // producer must be scheduled ahead of the phase's first
                // tagged step (pyfm hoists per phase).
                HADRONS_ERROR_REF(ObjectDefinition, "object '"
                              + env().getObjectName(a) + "' (address "
                              + std::to_string(a) + "), needed by tagged"
                              + " module '" + getModuleName(p[i])
                              + "' (step " + std::to_string(i + 1)
                              + ") and produced by module '" + getModuleName(m)
                              + "', does not exist at its split phase's"
                              + " pre-scatter: schedule the producer ahead of"
                              + " the phase's tagged block (or check that the"
                              + " producer creates its output)", a);
            }
            shadowConsumers_[a].insert(s);
            shadowLastUse_[{a, s}] = i;
            if (env().isLatticeObject(a) || env().isCloneableObject(a))
            {
                continue;   // scatter-leaf / clone: no producer descent
            }
            for (auto &in : module_[static_cast<unsigned int>(m)].input)
            {
                stack.push_back(in);
            }
        }
    }
}

void VirtualMachine::ensureShadowed(const std::vector<unsigned int> &inputs)
{
    for (auto &a: inputs)
    {
        if (shadowedObjects_.count(a) or batchQueued_.count(a)
            or rebuildQueued_.count(a))
        {
            continue;
        }
        int m = env().getObjectModule(a);
        if (m < 0)
        {
            continue;   // no producing module (e.g. externally provided object)
        }
        // Only producers that are TAGGED into a subgrid (subgrid >= 0) execute
        // on their subgrid and produce their output there directly, so they
        // need no shadow copy. UNTAGGED producers (solvers / fermion actions)
        // execute exactly once, at their original schedule step, on the WORLD
        // grid (with shadowCreateMode_ off), landing in the global store bound
        // to world-grid inputs. A subgrid consumer would otherwise fetch that
        // grid-mismatched global object (getDerivedObject finds no shadow entry
        // and falls through). Such producers MUST be rebuilt here on the subgrid
        // in shadow-create mode — so do NOT skip them (fall through below).
        if (module_[static_cast<unsigned int>(m)].subgrid >= 0)
        {
            continue;   // tagged producer: self-rebuilds on its subgrid
        }
        if (!env().hasCreatedObject(a))
        {
            continue;   // not created yet — will be shadowed by a later consumer
        }
        // ---- consumer-gated storage (scatter design D2) ----
        // shadowConsumers_[a] holds the COMPLETE set of subgrids consuming
        // a (computed by computeShadowConsumers before the pre-scatter
        // pass); a rank keeps its shadow copy only if its own sub-comm is
        // a consumer. A Grid_split yields each rank its OWN sub-comm's
        // piece, so non-consumers simply drop theirs. Empty set = never
        // pre-computed (per-step late path, nrhs == 1): store everywhere,
        // exactly the pre-gating behaviour.
        const bool keepShadow = shadowConsumers_[a].empty()
                                || shadowConsumers_[a].count(me_);

        // ---- batched sole-consumer scatter (R4, Slice 6) ----
        // sole-consumer scatterables are deferred: one vector Grid_split
        // per type/shape group delivers slot k to subgrid k without the
        // scalar wrapper's nrhs world-grid copies. Deferred rebuilds run
        // after the batches (their setup() reads lattice inputs through the
        // scope-aware fetch).
        if (collectingBatches_ and subGrid_.nrhs > 1
            and shadowConsumers_[a].size() == 1
            and env().isLatticeObject(a))
        {
            batchQueued_.insert(a);
            batchQueue_.push_back(a);
            continue;
        }

        // scatterable lattice value or container: the Grid_split inside
        // splitTo() is world-collective and runs on EVERY rank in lockstep
        // (the branch below is rank-independent); only the STORE is gated.
        auto scattered = env().splitObject(a);
        if (scattered)
        {
            if (keepShadow)
            {
                env().addShadowObject(a, env().getActiveSplitKey(),
                                      std::move(scattered));
            }
            shadowedObjects_.insert(a);
            continue;   // lattice / lattice container: scattered (or dropped)
        }
        if (env().isCloneableObject(a))
        {
            // cloneable metadata (e.g. std::vector<Integer>): the deep copy
            // is purely local (no collectives), so non-consumers skip it
            if (keepShadow)
            {
                env().cloneObject(a, env().getActiveSplitKey());
            }
            shadowedObjects_.insert(a);
            continue;   // no producer rebuild needed
        }
        // non-lattice: rebuild the producer on the subgrid into the shadow
        // store. In collection mode defer the rebuild until after the
        // batches: its setup() reads lattice inputs through the scope-aware
        // fetch, and a deferred sole-consumer input must be stored first.
        if (collectingBatches_)
        {
            // visit mark BEFORE the recursion (cycle break / memo guard —
            // unchanged); queue AFTER it: rebuildQueue_ must end up in
            // POST-order (every dependency queued ahead of every dependent
            // across all walk roots) so pass 3 can rebuild forward in
            // dependency order (Phase 14, the G/g-s4eps-e1000 fix — a
            // pre-order fill iterated in reverse rebuilds shared
            // dependencies AFTER their later dependents: six of seven
            // solvers before the shared actions, each binding the
            // world-grid original at setup() time).
            rebuildQueued_.insert(a);
            ensureShadowed(module_[static_cast<unsigned int>(m)].input);
            rebuildQueue_.push_back(a);
            continue;
        }
        // Mark before recursing to break dependency cycles.
        shadowedObjects_.insert(a);
        ensureShadowed(module_[static_cast<unsigned int>(m)].input);
        if (!keepShadow)
        {
            // The recursion above (world-collective scatters) ran on every
            // rank; the rebuild's own collectives are sub-comm-scoped
            // (setup/execute run under the split-grid scope on this rank's
            // sub-comm grids), so non-consuming ranks skip the work
            // entirely. If a future module's setup() introduced a WORLD
            // collective this skip would deadlock -- such a module would
            // already break the split-grid scope contract.
            continue;
        }
        env().enterShadowCreate(env().getActiveSplitKey());
        try
        {
            // Call both setup() and execute(). setup() alone fully constructs
            // grid-bound objects (solvers, fermion actions) whose execute() is
            // empty. But some non-lattice objects — notably sink functions
            // (std::function produced by MSink modules) — are only populated
            // during execute(), not setup(). Calling execute() here ensures
            // those objects are usable on the subgrid; for solver/action
            // modules execute() is a harmless no-op.
            // attribute failures to the producer being rebuilt, not the
            // consumer that triggered the shadowing — set before
            // resetShadowState() so even its failures attribute correctly
            // (Step-5 triage; design follow-up:
            // .rpiv/artifacts/designs/2026-10-15_18-30-55_single-pass-memory-profiler.md)
            currentModule_ = m;
            module_[static_cast<unsigned int>(m)].data->resetShadowState();
            module_[static_cast<unsigned int>(m)].data->setup();
            module_[static_cast<unsigned int>(m)].data->execute();
        }
        catch (Exceptions::ObjectDefinition &exc)
        {
            env().exitShadowCreate();
            currentModule_ = -1;
            rethrowWithContext(static_cast<unsigned int>(m), exc);
        }
        catch (...)
        {
            env().exitShadowCreate();
            currentModule_ = -1;
            throw;
        }
        env().exitShadowCreate();
    }
}

// one-time slot-delivery self-check (Phase 11) ///////////////////////////////////
// Scatters nrhs constant-filled marker lattices (slot k <- k) through the
// SAME gridSplitBatch the batches use; after the Phase 11 fix every rank's
// piece must carry its own sub-comm index. Guards the analytic permutation
// against layouts not hand-traced (3+ split dims, multi-rank subgrids on
// split dims) at the cost of one small-type scatter per program.
void VirtualMachine::verifyBatchSlotDelivery(void)
{
    using Mark = Lattice<iScalar<vInteger>>;
    using Sobj = typename Mark::scalar_object;   // per-site iScalar<int>
    auto *world = env().getWorldGrid();          // 4d world grid (scope-free:
                                                 // scope is ON here — pass 2
                                                 // of the batched pre-scatter)
    auto *sub   = env().template getGrid<vInteger>();   // active split's
                                                        // sub-grid (cached)
    const auto  nrhs = static_cast<unsigned int>(subGrid_.nrhs);
    std::vector<Mark>         fulls;
    std::vector<const Mark *> ptrs;
    fulls.reserve(nrhs);   // markers must never move: ptrs holds their
                           // addresses across the fill loop (emplace_back on
                           // an unreserved vector reallocates and dangles
                           // every earlier pointer — the C4j/2x2 abort)
    ptrs.reserve(nrhs);
    for (unsigned int k = 0; k < nrhs; ++k)
    {
        // marker k = the constant k, filled through the same lex-array
        // plumbing gridSplitBatch uses (no reliance on broadcast assignment)
        std::vector<Sobj> data(world->lSites());
        for (auto &s : data)
        {
            s() = static_cast<int>(k);
        }
        fulls.emplace_back(world);
        vectorizeFromLexOrdArray(data, fulls.back());
        ptrs.push_back(&fulls.back());
    }
    Mark piece(sub);
    gridSplitBatch(ptrs, piece);
    std::vector<Sobj> data;
    unvectorizeToLexOrdArray(data, piece);
    int got = (data.empty() ? -1 : static_cast<int>(data[0]()));
    for (auto &s : data)
    {
        if (static_cast<int>(s()) != got)
        {
            got = -1;   // not a constant: report as mis-delivery below
            break;
        }
    }
    if (got != me_)
    {
        HADRONS_ERROR(Definition, "batched split slot mis-delivery: subcomm "
                      + std::to_string(me_) + " received slot "
                      + std::to_string(got) + " (expected "
                      + std::to_string(me_) + "); gridSplitBatchSlotToSub"
                      + " disagrees with the delivered pairing on this grid"
                      + " pair — re-derive the collapse-arithmetic trace for"
                      + " this layout");
    }
    LOG(Message) << "Split-phase batch slot delivery verified (subcomm "
                 << me_ << " received slot " << me_ << ")." << std::endl;
}

// three-pass pre-scatter (R4) //////////////////////////////////////////////////
void VirtualMachine::executeBatchedPreScatter(const Program &p,
                                              const PhaseRange &range)
{
    // ---- pass 1: collect --------------------------------------------------
    // multi-consumer scatterables scatter immediately (the scalar wrapper's
    // nrhs-fold replication is load-bearing for multi-consumer delivery);
    // cloneables clone immediately (local); sole-consumer scatterables and
    // rebuild-path objects are queued.
    collectingBatches_ = true;
    for (unsigned int j = range.first; j <= range.last; ++j)
    {
        if (module_[p[j]].subgrid >= 0)
        {
            ensureShadowed(module_[p[j]].input);
        }
    }
    collectingBatches_ = false;

    // ---- pass 2: batch sole-consumer scatterables by concrete type --------
    if (!batchQueue_.empty())
    {
        if (!batchSlotDeliveryVerified_)
        {
            // once per program, prove slot k lands on sub-comm k on this grid
            // pair before trusting any batch. Rank-uniform site (pass-2 state
            // is rank-independent; every rank reaches it in lockstep at phase
            // entry), so the world-collective marker scatter is safe here.
            verifyBatchSlotDelivery();
            batchSlotDeliveryVerified_ = true;
        }
        // group by concrete derived type: same type name <=> same Holder<T>
        // <=> batchable together (shape compatibility checked by
        // splitBatchTo; nullptr -> per-object fallback for that batch)
        std::map<std::string, std::vector<unsigned int>> groups;
        for (auto a : batchQueue_)
        {
            groups[env().getObjectDerivedType(a)].push_back(a);
        }
        for (auto &kv : groups)
        {
            auto                  &members = kv.second;   // deduped, walk order
            std::set<unsigned int> placed;
            while (placed.size() < members.size())
            {
                // greedy slot fill: slot k <- an unplaced member whose sole
                // consumer is subgrid k. Slots are FRESH each iteration
                // (constructed here), so the first unplaced member always
                // fits and every member rides some batch. Pad short batches
                // with an in-batch pointer (R8) — the receiving subgrid
                // drops the duplicate.
                std::vector<unsigned int> slots(
                    subGrid_.nrhs, std::numeric_limits<unsigned int>::max());
                unsigned int anchor = std::numeric_limits<unsigned int>::max();
                bool progress = false;
                for (auto a : members)
                {
                    if (placed.count(a))
                    {
                        continue;
                    }
                    const int s = *shadowConsumers_.at(a).begin();
                    if (slots[s] == std::numeric_limits<unsigned int>::max())
                    {
                        slots[s] = a;
                        placed.insert(a);
                        if (anchor == std::numeric_limits<unsigned int>::max())
                        {
                            anchor = a;
                        }
                        progress = true;
                    }
                }
                if (!progress)
                {
                    // defensive (unreachable per the fresh-slots argument):
                    // scatter any unplaced member per-object so a queued
                    // object can NEVER be stranded — the memo guard in
                    // ensureShadowed would otherwise skip it forever
                    for (auto a : members)
                    {
                        if (placed.count(a))
                        {
                            continue;
                        }
                        const bool keep = shadowConsumers_.at(a).empty()
                                          || shadowConsumers_.at(a).count(me_);
                        auto scattered = env().splitObject(a);
                        if (scattered and keep)
                        {
                            env().addShadowObject(
                                a, env().getActiveSplitKey(),
                                std::move(scattered));
                        }
                        shadowedObjects_.insert(a);
                    }
                    break;
                }
                for (auto &s : slots)
                {
                    if (s == std::numeric_limits<unsigned int>::max())
                    {
                        s = anchor;   // R8 padding: receiver drops it
                    }
                }
                std::set<unsigned int> batchMembers(slots.begin(), slots.end());
                LOG(Debug) << "Split phase batch: " << batchMembers.size()
                           << " sole-consumer object(s) of type '" << kv.first
                           << "' into " << subGrid_.nrhs << " slots" << std::endl;
                // one collective per element/key; every rank derives the
                // same groups/slots/padding (rank-invariant inputs)
                auto piece = env().splitObjectBatch(slots);
                if (piece)
                {
                    const auto me = static_cast<unsigned int>(me_);
                    if (shadowConsumers_.at(slots[me]).count(me_))
                    {
                        env().addShadowObject(slots[me],
                                              env().getActiveSplitKey(),
                                              std::move(piece));
                    }
                }
                else
                {
                    // shape mismatch: per-object scalar fallback for THIS
                    // batch's members (replicate semantics per object)
                    for (auto a : batchMembers)
                    {
                        const bool keep = shadowConsumers_.at(a).empty()
                                          || shadowConsumers_.at(a).count(me_);
                        auto scattered = env().splitObject(a);
                        if (scattered and keep)
                        {
                            env().addShadowObject(
                                a, env().getActiveSplitKey(),
                                std::move(scattered));
                        }
                    }
                }
                // memo on EVERY rank for every member (incl. padded
                // duplicates): the memo must not diverge across ranks
                for (auto a : slots)
                {
                    shadowedObjects_.insert(a);
                }
            }
        }
    }

    // ---- pass 3: deferred rebuilds (dependency order) ---------------------
    // rebuildQueue_ is filled in POST-order DFS (visit mark before the
    // producer-input recursion, push after it), so every dependency of an
    // object is queued ahead of it across ALL walk roots — shared
    // dependencies included (queued at first discovery, before any later
    // dependent that reaches them). Forward iteration is therefore a valid
    // topological order of the closure DAG, mirroring the inline path's
    // recurse-then-rebuild ordering. (Phase 14: the pre-order fill +
    // reverse iteration this replaced rebuilt shared dependencies AFTER
    // their later dependents — with N solvers sharing one action, N-1
    // rebuilt before the action and bound the world-grid original at
    // setup() time.) Sub-comm-scoped work; non-consuming ranks skip
    // (existing contract).
    for (auto it = rebuildQueue_.begin(); it != rebuildQueue_.end(); ++it)
    {
        const unsigned int a = *it;
        const int          m = env().getObjectModule(a);
        if (m < 0 or module_[static_cast<unsigned int>(m)].subgrid >= 0)
        {
            continue;   // defensive: classification happened at collection
        }
        const bool keepShadow = shadowConsumers_.at(a).empty()
                                || shadowConsumers_.at(a).count(me_);
        if (!keepShadow)
        {
            continue;
        }
        env().enterShadowCreate(env().getActiveSplitKey());
        try
        {
            // attribute failures to the producer being rebuilt (same
            // contract as the inline path)
            currentModule_ = m;
            module_[static_cast<unsigned int>(m)].data->resetShadowState();
            module_[static_cast<unsigned int>(m)].data->setup();
            module_[static_cast<unsigned int>(m)].data->execute();
        }
        catch (Exceptions::ObjectDefinition &exc)
        {
            env().exitShadowCreate();
            currentModule_ = -1;
            rethrowWithContext(static_cast<unsigned int>(m), exc);
        }
        catch (...)
        {
            env().exitShadowCreate();
            currentModule_ = -1;
            throw;
        }
        env().exitShadowCreate();
        shadowedObjects_.insert(a);
    }

    // per-phase transient state outlives its use only until here
    batchQueue_.clear();
    batchQueued_.clear();
    rebuildQueue_.clear();
    rebuildQueued_.clear();
}

void VirtualMachine::executeProgram(const Program &p)
{
    Size            memPeak = 0, sizeBefore, sizeAfter;
    GarbageSchedule freeProg;

    // Fail fast on invalid programs: loaded .sched files (and schedules
    // restored from the application database) are otherwise never checked,
    // and a mis-ordered or incomplete program would only surface as a
    // mysterious mid-run object error — or silently skip modules.
    validateProgram(p);
    // Ensure the memory profile (and the objectIsLattice_ cache it populates)
    // is up to date. On the naive path this builds the profile by walking
    // THIS program, so the profiling pass frees exactly as this run will; on
    // the genetic and loaded-schedule paths the profile is already cached
    // (schedule() pre-seed / printSchedule->memoryNeeded) and this call is a
    // no-op — the profile content (per-module sizes, objectIsLattice_) is
    // order-independent either way. makeGarbageSchedule relies on
    // objectIsLattice_ to avoid pinning unrelated heavy objects (e.g.
    // eigenvector packs) reached through a scattered lattice.
    getMemoryProfile(p);
    // build garbage collection schedule
    LOG(Debug) << "Building garbage collection schedule..." << std::endl;
    freeProg = makeGarbageSchedule(p);
    for (unsigned int i = 0; i < freeProg.size(); ++i)
    {
        std::string msg = "";

        for (auto &a: freeProg[i])
        {
            msg += env().getObjectName(a) + " ";
        }
        msg += "]";
        LOG(Debug) << std::setw(4) << i + 1 << ": [" << msg << std::endl;
    }

    // program execution
    LOG(Debug) << "Executing program..." << std::endl;
    totalTime_ = GridTime::zero();
    moduleTimeProfile_.clear();
    moduleTypeTimeProfile_.clear();

    // Split-phase boundaries: ordered per-phase step ranges (empty when no
    // tagged modules or split unconfigured — the legacy single-window case
    // is then simply one range). All world-collective Grid_split calls are
    // front-loaded into a pre-scatter pass at each phase's entry so per-step
    // ensureShadowed becomes a no-op within the phase, eliminating
    // world-collective barriers from the tagged-step loop. Phase indices are
    // validated non-decreasing over tagged modules (untagged modules may
    // interleave freely within and between phases).
    const auto phaseRanges = computePhaseRanges(p);
    unsigned int phaseIdx    = 0;   // active index into phaseRanges
    // per-phase shadow free schedule (R2, Slice 5): shadowFreeSteps_[k] lists
    // {address, subgrid} entries whose LAST consumer is step range.first + k;
    // rebuilt (assign -- no stale per-phase content) at each phase entry.
    std::vector<std::vector<std::pair<unsigned int, int>>> shadowFreeSteps_;

    for (unsigned int i = 0; i < p.size(); ++i)
    {
        // ---- Split-phase exit: join barrier + per-phase VM state reset ----
        // Synchronize all sub-comms before any following global step or the
        // next phase's pre-scatter. During the phase the sub-comms ran
        // concurrently (no per-step world barriers); this is the "join".
        // Per-phase VM state (shadow memo, consumer sets) is reset HERE so
        // the next phase's pre-scatter cannot short-circuit on stale marks
        // (which would silently fall through to grid-mismatched world
        // objects). shadowStore_ entries of this phase's key are swept by
        // the shadow GC extension (Slice 5); entries of old keys are
        // invisible to lookups (keyed by the active phase).
        if (!phaseRanges.empty() and phaseIdx < phaseRanges.size()
            and i == phaseRanges[phaseIdx].last + 1)
        {
            env().setSubGridScope(false);
            env().getGrid()->Barrier();
            shadowedObjects_.clear();
            shadowConsumers_.clear();
            // safety net for shadows whose last-use free never fired (B1)
            env().sweepShadows(static_cast<int>(phaseRanges[phaseIdx].phase));
            // stale prior-phase entries would otherwise survive into the
            // next phase's shadowFreeSteps_ build (K1)
            shadowLastUse_.clear();
            phaseIdx++;
        }

        // ---- Split-phase entry: pre-scatter THIS phase's tagged inputs ----
        if (!phaseRanges.empty() and phaseIdx < phaseRanges.size()
            and i == phaseRanges[phaseIdx].first)
        {
            const auto &range = phaseRanges[phaseIdx];
            const auto  key   = static_cast<int>(range.phase);
            // switch the shadow split key to this phase (flushes subgrid grid
            // caches on key change; does NOT wipe the shadow store)
            env().setActiveSubGrid(subGrid_.grid.get(), key);
            env().setSubGridScope(true);
            LOG(Message) << "Pre-scattering subgrid inputs for split phase "
                         << key << " (steps " << range.first + 1 << "-"
                         << range.last + 1 << ")..." << std::endl;
            // consumer-gated shadow storage phase 1: which sub-comm consumes
            // each tagged-closure object of THIS phase (must complete before
            // any store decision), enforcing the pre-scatter contract (Q7)
            computeShadowConsumers(p, range);
            // rebuild the per-step shadow free schedule for this phase (R2):
            // shadowFreeSteps_[i - range.first] lists {address, subgrid}
            // entries whose last consumer is step i.
            shadowFreeSteps_.assign(range.last - range.first + 1, {});
            for (auto &as : shadowLastUse_)
            {
                shadowFreeSteps_[as.second - range.first].push_back(as.first);
            }
            // three-pass pre-scatter (R4): collect -> batch -> rebuild
            executeBatchedPreScatter(p, range);
            LOG(Message) << "Pre-scatter complete; sub-comms now execute "
                         << "concurrently." << std::endl;
        }

        // ---- Per-step scope gate ----
        if (splitConfigured_)
        {
            if (module_[p[i]].subgrid >= 0)
            {
                // scope ON for this tagged module's execution
                env().setSubGridScope(true);
                // NO per-step ensureShadowed: all shadowing happened at this
                // phase's entry pre-scatter (executeBatchedPreScatter), where
                // every rank joins the world-collective Grid_split calls in
                // lockstep. Once the shadow GC starts freeing spent shadows
                // (owner-local memo erases), a per-step call would scatter on
                // SOME ranks only — ranks already waiting in the phase-exit
                // Barrier never match the collective: deadlock. Mid-phase
                // shadowing is forbidden; the pre-scatter is authoritative.
                if (module_[p[i]].subgrid >= 0 && module_[p[i]].subgrid != me_)
                {
                    env().setSubGridScope(false);
                    // world-object GC is rank-uniform (Phase 13, the C4h
                    // fix): freeProg is schedule-derived and freeSet is
                    // local — a non-owning rank that skips this free keeps
                    // its piece of these world fields until program end.
                    // Shadow frees below are owner-local no-ops here.
                    env().freeSet(freeProg[i]);
                    continue;   // skip non-owning subcomm — no world barrier!
                }
            }
            else
            {
                // global module: must execute on the world grid
                env().setSubGridScope(false);
            }
        }
        // execute module
        LOG(Message) << SEP << " Measurement step " << i + 1 << "/"
                     << p.size() << " (module '" << module_[p[i]].name
                     << "') " << SEP << std::endl;
        LOG(Message) << SMALL_SEP << " Module execution" << std::endl;
        currentModule_ = p[i];
        try
        {
            (*module_[p[i]].data)();
        }
        catch (Exceptions::ObjectDefinition &exc)
        {
            // dependency failures during execution get the same context as
            // during profiling: module, its outputs, failing object, producer
            rethrowWithContext(p[i], exc);
        }
        currentModule_ = -1;
        sizeBefore = env().getTotalSize();
        // print time profile after execution
        LOG(Message) << SMALL_SEP << " Timings" << std::endl;

        std::map<std::string, GridTime> ctiming, gtiming;
        GridTime                        total;

        ctiming  = module_[p[i]].data->getTimings();
        total    = ctiming.at("_total");
        gtiming["total"]     = ctiming["_total"];   ctiming.erase("_total");
        gtiming["setup"]     = ctiming["_setup"];   ctiming.erase("_setup");
        gtiming["execution"] = ctiming["_execute"]; ctiming.erase("_execute");
        LOG(Message) << "* GLOBAL TIMERS" << std::endl;
        printTimeProfile(gtiming, total);
        if (!ctiming.empty())
        {
            LOG(Message) << "* CUSTOM TIMERS" << std::endl;
            printTimeProfile(ctiming, total);
        }
        moduleTimeProfile_[module_[p[i]].name] = total;
        std::string moduleType = getModuleType(p[i]);
        if (moduleTypeTimeProfile_.find(moduleType) == moduleTypeTimeProfile_.end())
        {
            moduleTypeTimeProfile_[getModuleType(p[i])] = total;
        }
        else
        {
            moduleTypeTimeProfile_.at(getModuleType(p[i])) += total;
        }
        totalTime_ += total;
        // print used memory after execution
        LOG(Message) << SMALL_SEP << " Memory management" << std::endl;
        MemoryUtils::printMemory();
        if (sizeBefore > memPeak)
        {
            memPeak = sizeBefore;
        }
        // garbage collection for step i
        LOG(Message) << "Garbage collection..." << std::endl;
        env().freeSet(freeProg[i]);

        // ---- Shadow GC (R2, Slice 5): free shadows after last consumer ----
        // Local to the owning sub-comm: consumer sets are rank-independent,
        // so every rank derives the same free step; only ranks with
        // me_ == s hold a {address, key} entry. No collectives.
        if (!phaseRanges.empty() and phaseIdx < phaseRanges.size()
            and i >= phaseRanges[phaseIdx].first
            and i <= phaseRanges[phaseIdx].last)
        {
            const auto key = static_cast<int>(phaseRanges[phaseIdx].phase);
            for (auto &as : shadowFreeSteps_[i - phaseRanges[phaseIdx].first])
            {
                if (me_ == as.second)
                {
                    env().freeShadowObject(as.first, key);
                    shadowedObjects_.erase(as.first);
                }
            }
        }

        // Clean up remaining temporary objects
        for (unsigned int a = 0; a < env().getMaxAddress(); ++a)
        {
            if (env().getObjectStorage(a) == Environment::Storage::temporary)
            {
                env().freeObject(a);
            }
        }

        // print used memory after garbage collection if necessary
        sizeAfter = env().getTotalSize();
        if (sizeBefore != sizeAfter)
        {
            MemoryUtils::printMemory();
        }
        else
        {
            LOG(Message) << "Nothing to free" << std::endl;
        }
    }
    // Join barrier (final phase): synchronize before cleanup and any
    // post-executeProgram MPI operations, exactly as HEAD did, plus the
    // final wholesale teardown (scope, shadow store, subgrid state).
    if (!phaseRanges.empty())
    {
        env().setSubGridScope(false);
        env().getGrid()->Barrier();
    }
    shadowConsumers_.clear();   // (already per-phase; final safety)
    shadowLastUse_.clear();     // (already per-phase; final safety)
    // C2 fix: clear scope whenever it was set (not position-based)
    if (splitConfigured_ && env().isSubGridActive())
    {
        env().clearActiveSubGrid();
        shadowedObjects_.clear();
    }
    // print total time profile
    LOG(Message) << SEP << " Measurement time profile" << SEP << std::endl;
    LOG(Message) << "Total measurement time: " << timeString(totalTime_) << std::endl;
    LOG(Message) << SMALL_SEP << " Module breakdown" << std::endl;
    printTimeProfile(moduleTimeProfile_, totalTime_);
    LOG(Message) << SMALL_SEP << " Module type breakdown" << std::endl;
    printTimeProfile(moduleTypeTimeProfile_, totalTime_);
}

void VirtualMachine::executeProgram(const std::vector<std::string> &p)
{
    Program pAddress;
    
    for (auto &n: p)
    {
        pAddress.push_back(getModuleAddress(n));
    }
    executeProgram(pAddress);
}

// generate result DB //////////////////////////////////////////////////////////
void VirtualMachine::generateResultDb(void)
{
    for (auto &m: module_)
    {
        m.data->generateResultDb();
    }
}
