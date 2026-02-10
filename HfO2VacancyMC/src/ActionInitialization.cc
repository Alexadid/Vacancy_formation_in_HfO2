#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

#include "G4RunManager.hh"

ActionInitialization::ActionInitialization(DetectorConstruction* det) : fDet(det) {}

void ActionInitialization::Build() const {
    SetUserAction(new PrimaryGeneratorAction());
    SetUserAction(new RunAction(fDet));
    SetUserAction(new SteppingAction(fDet));
}
