#include "SteppingAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

SteppingAction::SteppingAction(DetectorConstruction* det) : fDet(det) {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    const auto edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.0) return;

    const auto pre = step->GetPreStepPoint();
    const auto vol = pre->GetTouchableHandle()->GetVolume();
    if (!vol) return;

    // Only score inside HfO2 volume
    if (vol->GetName() != "HfO2PV") return;

    // Use mid-step position for binning
    const auto p1 = pre->GetPosition();
    const auto p2 = step->GetPostStepPoint()->GetPosition();
    const auto pmid = 0.5*(p1 + p2);

    fDet->GetVoxelGrid().AddEdep(pmid, edep);
}
