#pragma once
#include "G4UserSteppingAction.hh"

class DetectorConstruction;

class SteppingAction : public G4UserSteppingAction {
public:
    explicit SteppingAction(DetectorConstruction* det);
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step* step) override;

private:
    DetectorConstruction* fDet = nullptr;
};
