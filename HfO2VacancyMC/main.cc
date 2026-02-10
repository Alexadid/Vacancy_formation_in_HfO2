#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4VModularPhysicsList.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmParameters.hh"
#include "G4DecayPhysics.hh"

int main(int argc, char** argv) {
    auto runManager = new G4RunManager();

    auto det = new DetectorConstruction();
    runManager->SetUserInitialization(det);

    // Physics list: minimal modular + Livermore EM
    auto phys = new G4VModularPhysicsList();
    phys->RegisterPhysics(new G4DecayPhysics());
    phys->RegisterPhysics(new G4EmLivermorePhysics());

    // Optional: lower edge for cuts (important if you later want more low-energy secondaries)
    // You can also set via macro: /cuts/setLowEdge 100 eV
    auto emParams = G4EmParameters::Instance();
    emParams->SetMinEnergy(100.0*eV);
    emParams->SetMaxEnergy(100.0*MeV);

    runManager->SetUserInitialization(phys);
    runManager->SetUserInitialization(new ActionInitialization(det));

    // runManager->Initialize();

    // UI
    auto UImanager = G4UImanager::GetUIpointer();
    if (argc == 1) {
        auto ui = new G4UIExecutive(argc, argv);
        auto vis = new G4VisExecutive();
        vis->Initialize();
        UImanager->ApplyCommand("/control/execute macros/run.mac");
        ui->SessionStart();
        delete vis;
        delete ui;
    } else {
        G4String macro = argv[1];
        UImanager->ApplyCommand("/control/execute " + macro);
    }

    delete runManager;
    return 0;
}
