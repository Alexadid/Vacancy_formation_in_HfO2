#pragma once
#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

class ActionInitialization : public G4VUserActionInitialization {
public:
    explicit ActionInitialization(DetectorConstruction* det);
    ~ActionInitialization() override = default;

    void Build() const override;

private:
    DetectorConstruction* fDet = nullptr;
};