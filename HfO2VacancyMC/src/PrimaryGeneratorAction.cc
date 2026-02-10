#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction() {
    fGun = new G4ParticleGun(1);

    auto electron = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    fGun->SetParticleDefinition(electron);

    // Direction: along -Z into the stack (vacuum -> HfO2 -> Si)
    fGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1));

    // Default energy
    fGun->SetParticleEnergy(fEnergyKeV * keV);

    // Start position slightly above top surface (z=0 is top of HfO2)
    fGun->SetParticlePosition(G4ThreeVector(0,0, +50*nm));

    fMessenger = new G4GenericMessenger(this, "/gun/", "Gun control");
    fMessenger->DeclareProperty("energyKeV", fEnergyKeV, "Primary electron energy in keV");
    fMessenger->DeclareProperty("sigmaNm", fBeamSigmaNm, "Gaussian beam sigma in nm (0 => pencil)");
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fGun;
    delete fMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    fGun->SetParticleEnergy(fEnergyKeV * keV);

    // Beam spot (optional Gaussian)
    double x = 0.0, y = 0.0;
    if (fBeamSigmaNm > 0.0) {
        const double sigma = fBeamSigmaNm * nm;
        x = G4RandGauss::shoot(0.0, sigma);
        y = G4RandGauss::shoot(0.0, sigma);
    }
    fGun->SetParticlePosition(G4ThreeVector(x, y, +50*nm));

    fGun->GeneratePrimaryVertex(anEvent);
}
