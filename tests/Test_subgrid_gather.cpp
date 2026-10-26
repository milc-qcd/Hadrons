/*
 * Test_subgrid_gather.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2026
 *
 * Author: Michael Lynch
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

// Subgrid gather regression test: a tagged (subgrid) module's output is
// gathered back to the world grid at the producing phase's exit and consumed
// by an untagged world module. The gathered field must be BITWISE equal to
// the same field produced by an untagged world module: the gather is pure
// data movement and the plane-wave source content is deterministic and
// site-local. Run with the production multi-dim split layout: world
// --mpi 1 1 2 2 with mpiSplit 1 1 1 1 -> 4 sub-comms of one rank each,
// ratio (1,1,2,2), for which the unsplit slots 1<->2 are permuted
// (non-trivial subToSlot inverse, exercised exactly as load-bearing in
// production).

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

// the round-trip compare below loads LatticePropagator while the save
// modules write FIMPL::PropagatorField: tie the two types at compile time
// so a precision-macro change (FIMPL = HADRONS_IMPL(WilsonImpl, R),
// Global.hpp:103-107) fails here, not silently at compare time
static_assert(std::is_same<LatticePropagator, FIMPL::PropagatorField>::value,
              "LatticePropagator must match FIMPL::PropagatorField for the "
              "gather round-trip compare");

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;

    // run setup ///////////////////////////////////////////////////////////////
    // trajectory number feeds the result file names: <stem>.<traj>.bin
    const unsigned int      traj = 1500;
    Application             application;
    Application::GlobalPar  globalPar;

    globalPar.trajCounter.start = traj;
    globalPar.trajCounter.end   = traj + 20;
    globalPar.trajCounter.step  = 20;
    globalPar.runId             = "test_subgrid_gather";
    // production multi-dim split layout (see file header comment)
    globalPar.split.mpiSplit    = "1 1 1 1";
    application.setPar(globalPar);

    // untagged world reference: the same deterministic plane-wave source,
    // saved straight from the world grid
    MSource::Momentum::Par momPar;
    momPar.mom = "0 0 0 1";
    application.createModule<MSource::Momentum>("mom_ref", momPar);
    MIO::SavePropagator::Par saveRefPar;
    saveRefPar.name     = "mom_ref";
    saveRefPar.fileStem = "subgrid_gather_world_reference";
    application.createModule<MIO::SavePropagator>("save_ref", saveRefPar);

    // tagged producer: the same source on subgrid 0, split phase 0. Created
    // through the string/XML overload because only it carries the
    // <subgrid>/<splitPhase> tags (the template createModule has no such
    // parameters).
    const std::string momTagXml
        = "<module><options><mom>0 0 0 1</mom></options></module>";
    XmlReader         momTagReader(momTagXml, true, "module");
    application.createModule("mom_tag", "MSource::Momentum", momTagReader,
                             /* subgrid    */ 0,
                             /* splitPhase */ 0);

    // untagged POST-PHASE world consumer of the tagged output: this
    // dependency is what the VM's gather schedule pre-pass detects, so the
    // phase-exit gather installs a world twin of "mom_tag" on every rank
    // before this module runs
    MIO::SavePropagator::Par saveGatPar;
    saveGatPar.name     = "mom_tag";
    saveGatPar.fileStem = "subgrid_gather_gathered";
    application.createModule<MIO::SavePropagator>("save_gathered",
                                                  saveGatPar);

    // execution ///////////////////////////////////////////////////////////////
    application.run();

    // check ///////////////////////////////////////////////////////////////////
    // Load both saved fields on the world grid and require bitwise equality:
    // the gather is redistribution-only, so any difference indicts the slot
    // permutation, the descriptor broadcast, or the twin installation. The
    // saved metadata must also report the original producer module for each
    // field (the gather replaces content in place -- it does not
    // re-attribute the object to a new module).
    auto                     &env  = Environment::getInstance();
    GridBase                 *grid = env.getGrid<vComplex>();
    LatticePropagator         gathered(grid), reference(grid);
    MIO::GenericFieldMetadata mdGathered, mdReference;
    const std::string         gatheredFile
        = ModuleBase::resultFilename("subgrid_gather_gathered", traj, "bin");
    const std::string         referenceFile
        = ModuleBase::resultFilename("subgrid_gather_world_reference", traj,
                                     "bin");

    FieldReader<LatticePropagator> rGathered(gatheredFile, grid);
    rGathered.readField(gathered, mdGathered);
    rGathered.close();
    FieldReader<LatticePropagator> rReference(referenceFile, grid);
    rReference.readField(reference, mdReference);
    rReference.close();

    gathered -= reference;
    auto diff = norm2(gathered);
    LOG(Message) << "Gather round-trip norm2(gathered - reference) = "
                 << diff << std::endl;
    LOG(Message) << "Gathered field producer module:  '"
                 << mdGathered.moduleName << "' (expected 'mom_tag')"
                 << std::endl;
    LOG(Message) << "Reference field producer module: '"
                 << mdReference.moduleName << "' (expected 'mom_ref')"
                 << std::endl;

    bool ok = (diff == 0.);
    ok = ok and (mdGathered.moduleName == "mom_tag");
    ok = ok and (mdReference.moduleName == "mom_ref");
    if (!ok)
    {
        LOG(Error) << "Subgrid gather regression test FAILED" << std::endl;
        Grid_finalize();
        return EXIT_FAILURE;
    }
    LOG(Message) << "Subgrid gather regression test passed" << std::endl;

    // epilogue ////////////////////////////////////////////////////////////////
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();

    return EXIT_SUCCESS;
}
