#!/usr/bin/env python3
"""DDG4 (DD4hep + Geant4) steering for the SoLID DDVCS configuration.

Full Geant4 simulation of DDVCS events (HepMC3 input from solid_ddvcs_gen) through
geometry/solid_ddvcs.xml, writing EDM4hep/podio output that SoLIDEDM4hep_source reads.

Usage (inside eic-shell or any environment with DD4hep, Geant4, podio, edm4hep and solid_dd4hep):
    export SOLID_DD4HEP=/path/to/solid_dd4hep        # compact/ and libsolid_dd4hep
    export DETECTOR_PATH=/path/to/SoLID-Jana-EICRecon/geometry
    solid_ddvcs_gen -n 10000 -o ddvcs.hepmc
    python3 scripts/run_ddsim_ddvcs.py --compact $DETECTOR_PATH/solid_ddvcs.xml -i ddvcs.hepmc -n 10000 -o ddvcs_sim
    jana -Pplugins=solid_ddvcs -Psolid:source=edm4hep -Psolid:luminosity=1e38 ddvcs_sim.edm4hep.root

The generator already places the vertex inside the target (z = -315 +- 7.5 cm), so no
additional vertex offset is applied; a small transverse smearing models the beam spot.
"""
import argparse
import os
import sys

parser = argparse.ArgumentParser(prog="run_ddsim_ddvcs", description="SoLID DDVCS Geant4 simulation")
parser.add_argument("--compact", default=os.path.join(os.environ.get("DETECTOR_PATH", "geometry"), "solid_ddvcs.xml"))
parser.add_argument("-i", "--input", required=True, help="HepMC3 ASCII input (solid_ddvcs_gen output)")
parser.add_argument("-o", "--output", default="ddvcs_sim", help="output file stem (podio/edm4hep)")
parser.add_argument("-n", "--n-events", type=int, default=1000)
parser.add_argument("-s", "--n-skip", type=int, default=0)
parser.add_argument("--seed", type=int, default=854321)
parser.add_argument("--beam-sigma-mm", type=float, default=0.5, help="transverse beam spot sigma [mm]")
parser.add_argument("--physics", default="QGSP_BERT")
parser.add_argument("--vis", action="store_true")
parser.add_argument("-v", "--verbose", type=int, default=0)
args = parser.parse_args()

import DDG4
from g4units import GeV, MeV, keV, mm, ns


def run():
    if not args.vis:
        os.environ["G4UI_USE_TCSH"] = "0"
    kernel = DDG4.Kernel()
    description = kernel.detectorDescription()
    kernel.loadGeometry("file:" + args.compact)
    DDG4.importConstants(description)

    geant4 = DDG4.Geant4(kernel)
    if args.verbose > 0:
        geant4.printDetectors()
    if args.vis:
        geant4.setupUI(typ="qt", vis=True, ui=True, macro="macro/vis.mac")
    else:
        geant4.setupUI(typ="tcsh", vis=False, ui=False, macro=False)
    kernel.NumEvents = args.n_events

    geant4.setupTrackingField(stepper="ClassicalRK4", equation="Mag_UsualEqRhs")

    rndm = DDG4.Action(kernel, "Geant4Random/Random")
    rndm.Seed = args.seed
    rndm.initialize()

    # --- output: podio / EDM4hep
    podio = DDG4.EventAction(kernel, "Geant4Output2Podio/RootOutput", True)
    podio.HandleMCTruth = True
    podio.Control = True
    podio.Output = args.output
    podio.enableUI()
    kernel.eventAction().adopt(podio)

    # --- generator: HepMC3 file from solid_ddvcs_gen
    gen = DDG4.GeneratorAction(kernel, "Geant4GeneratorActionInit/GenerationInit")
    kernel.generatorAction().adopt(gen)

    reader = DDG4.GeneratorAction(kernel, "Geant4InputAction/hepmc3")
    reader.Mask = 0
    reader.Input = "HEPMC3FileReader|" + args.input
    reader.Sync = args.n_skip
    reader.enableUI()
    kernel.generatorAction().adopt(reader)

    smear = DDG4.GeneratorAction(kernel, "Geant4InteractionVertexSmear/SmearVert")
    smear.Mask = 0
    smear.Offset = (0 * mm, 0 * mm, 0 * mm, 0 * ns)
    smear.Sigma = (args.beam_sigma_mm * mm, args.beam_sigma_mm * mm, 0 * mm, 0 * ns)
    kernel.generatorAction().adopt(smear)

    merger = DDG4.GeneratorAction(kernel, "Geant4InteractionMerger/InteractionMerger")
    kernel.generatorAction().adopt(merger)
    primary = DDG4.GeneratorAction(kernel, "Geant4PrimaryHandler/PrimaryHandler")
    kernel.generatorAction().adopt(primary)

    part = DDG4.GeneratorAction(kernel, "Geant4ParticleHandler/ParticleHandler")
    part.SaveProcesses = ["Decay"]
    part.MinimalKineticEnergy = 10 * MeV
    part.KeepAllParticles = False
    part.enableUI()
    kernel.generatorAction().adopt(part)

    # --- sensitive detectors
    f_optical_reject = DDG4.Filter(kernel, "ParticleRejectFilter/OpticalPhotonRejector")
    f_optical_reject.particle = "opticalphoton"
    f_optical_select = DDG4.Filter(kernel, "ParticleSelectFilter/OpticalPhotonSelector")
    f_optical_select.particle = "opticalphoton"
    kernel.registerGlobalFilter(f_optical_reject)
    kernel.registerGlobalFilter(f_optical_select)

    geant4.setupTracker("GEMTracker_DDVCS")
    geant4.setupTracker("SPDLargeAngle")
    geant4.setupTracker("SPDForwardAngle")
    geant4.setupTracker("MRPCForwardAngle")
    geant4.setupTracker("MuonFAScint")
    geant4.setupTracker("MuonLAScint")
    geant4.setupTracker("MuonBarrelScint")
    for calo in ("FAECPreShower", "FAECShower", "LAECPreShower", "LAECShower"):
        geant4.setupCalorimeter(calo)
    for cher in ("LightGasCherenkov", "HeavyGasCherenkov"):
        seq, act = geant4.setupDetector(cher, "PhotoMultiplierSDAction")
        act.adopt(f_optical_select)

    # --- physics
    phys = geant4.setupPhysics(args.physics)
    geant4.addPhysics("Geant4PhysicsList/Myphysics")
    ph = DDG4.PhysicsList(kernel, "Geant4OpticalPhotonPhysics/OpticalPhotonPhys")
    phys.adopt(ph)
    ph = DDG4.PhysicsList(kernel, "Geant4CerenkovPhysics/CerenkovPhys")
    ph.MaxNumPhotonsPerStep = 10
    ph.MaxBetaChangePerStep = 10.0
    ph.TrackSecondariesFirst = True
    phys.adopt(ph)
    rg = geant4.addPhysics("Geant4DefaultRangeCut/GlobalRangeCut")
    rg.RangeCut = 0.7 * mm
    if args.verbose > 1:
        phys.dump()

    kernel.configure()
    kernel.initialize()
    kernel.run()
    kernel.terminate()


if __name__ == "__main__":
    run()
