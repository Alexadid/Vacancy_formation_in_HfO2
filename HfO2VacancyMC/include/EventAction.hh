#pragma once
#include "G4UserEventAction.hh"

class DetectorConstruction;

class EventAction : public G4UserEventAction {
public:
    explicit EventAction(DetectorConstruction* det);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

private:
    DetectorConstruction* fDet = nullptr;
};
