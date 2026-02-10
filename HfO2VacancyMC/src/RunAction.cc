#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

RunAction::RunAction(const DetectorConstruction* det) : fDet(det) {}

void RunAction::BeginOfRunAction(const G4Run*) {
    // Reset edep accumulators at start of run
    const_cast<DetectorConstruction*>(fDet)->GetVoxelGrid().ResetAccumulators();
}

void RunAction::EndOfRunAction(const G4Run*) {
    // Export energy deposition per voxel at end of run
    const auto& grid = fDet->GetVoxelGrid();
    grid.ExportCSV(fOutCsv);
}
