#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "VoxelGrid.hh"
#include "VacancyModel.hh"

EventAction::EventAction(DetectorConstruction* det) : fDet(det) {}

void EventAction::BeginOfEventAction(const G4Event*) {
    fDet->GetVoxelGrid().ResetEventAccumulators();
}

void EventAction::EndOfEventAction(const G4Event*) {
    // Process this event's deposition into vacancy state
    fDet->GetVacancyModel().ProcessEvent(fDet->GetVoxelGrid());
}
