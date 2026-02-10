#pragma once
#include "G4UserRunAction.hh"
#include <string>

class DetectorConstruction;

class RunAction : public G4UserRunAction {
public:
    explicit RunAction(const DetectorConstruction* det);
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

private:
    const DetectorConstruction* fDet = nullptr;
    std::string fOutCsv = "hfO2_edep_voxels.csv";
};
