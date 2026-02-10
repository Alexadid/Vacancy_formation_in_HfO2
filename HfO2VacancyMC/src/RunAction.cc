#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "G4Run.hh"

RunAction::RunAction(DetectorConstruction* det) : fDet(det) {}

void RunAction::BeginOfRunAction(const G4Run*) {
    fDet->GetVoxelGrid().ResetRunAccumulators();
    fDet->GetVoxelGrid().ResetEventAccumulators();
    fDet->GetVacancyModel().ResetToSeedOnly();
}

void RunAction::EndOfRunAction(const G4Run* run) {
    const auto& grid = fDet->GetVoxelGrid();
    const auto& vac    = fDet->GetVacancyModel();

    grid.ExportEdepCSV("hfO2_edep_voxels.csv");
    vac.ExportVacancyCSV("hfO2_vacancy_map.csv", grid);
    vac.ExportSummaryCSV("hfO2_vacancy_summary.csv", run->GetNumberOfEvent());
}
