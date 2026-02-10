#pragma once
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GenericMessenger.hh"
#include "G4ParticleGun.hh"

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* anEvent) override;

private:
    G4ParticleGun* fGun = nullptr;
    G4GenericMessenger* fMessenger = nullptr;

    double fEnergyKeV = 10.0;       // 1..30 keV
    double fBeamSigmaNm = 0.0;      // 0 => pencil beam; set >0 for Gaussian spot
};
