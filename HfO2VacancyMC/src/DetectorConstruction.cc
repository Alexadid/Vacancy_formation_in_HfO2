#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Region.hh"
#include "G4ProductionCuts.hh"

DetectorConstruction::DetectorConstruction() {
    fMessenger = new G4GenericMessenger(this, "/det/", "Detector control");

    fMessenger->DeclareProperty("hfo2ThicknessNm", fHfO2ThicknessNm, "HfO2 thickness in nm");
    fMessenger->DeclareProperty("padSizeUm", fPadSizeUm, "Pad lateral size in um (square)");
    fMessenger->DeclareProperty("siThicknessUm", fSiThicknessUm, "Si thickness in um");

    fMessenger->DeclareProperty("voxelDxNm", fVoxelDxNm, "Voxel size X in nm (HfO2 scoring grid)");
    fMessenger->DeclareProperty("voxelDyNm", fVoxelDyNm, "Voxel size Y in nm (HfO2 scoring grid)");
    fMessenger->DeclareProperty("voxelDzNm", fVoxelDzNm, "Voxel size Z in nm (HfO2 scoring grid)");

    fMessenger->DeclareProperty("vacConcCm3", fVacancy.GetParams().initConc_cm3, "Initial oxygen vacancy concentration in cm^-3");
    fMessenger->DeclareProperty("vacSeed", fVacancy.GetParams().initSeed, "Seed for vacancy initialization");
    fMessenger->DeclareProperty("hfo2Rho_g_cm3", fVacancy.GetParams().rho_g_cm3, "HfO2 density in g/cm3 (affects max vacancy capacity)");
}

void DetectorConstruction::DefineMaterials() {
    auto nist = G4NistManager::Instance();

    // Vacuum
    nist->FindOrBuildMaterial("G4_Galactic");

    // Silicon
    nist->FindOrBuildMaterial("G4_Si");

    // HfO2 (build explicitly: Hf + 2*O)
    // Density: set a reasonable default, can be tuned; ~9.68 g/cm3 often used for crystalline HfO2.
    auto elHf = nist->FindOrBuildElement("Hf");
    auto elO    = nist->FindOrBuildElement("O");

    const G4double density = 9.68 * g/cm3;
    auto HfO2 = new G4Material("HfO2", density, 2);
    HfO2->AddElement(elHf, 1);
    HfO2->AddElement(elO,    2);
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
    DefineMaterials();
    auto nist = G4NistManager::Instance();

    const auto padXY = fPadSizeUm * um;
    const auto tSi     = fSiThicknessUm * um;
    const auto tHf     = fHfO2ThicknessNm * nm;

    // World: enough space above and below
    const auto worldXY = 8.0 * um;
    const auto worldZ    = (tSi + tHf + 6.0*um);

    auto solidWorld = new G4Box("World", worldXY/2, worldXY/2, worldZ/2);
    auto worldMat = nist->FindOrBuildMaterial("G4_Galactic");
    fLogicWorld = new G4LogicalVolume(solidWorld, worldMat, "WorldLV");
    auto physWorld = new G4PVPlacement(nullptr, {}, fLogicWorld, "WorldPV", nullptr, false, 0);

    // Coordinate convention:
    // Put HfO2 on top of Si, vacuum above HfO2.
    // Let the top surface of HfO2 be at z = 0.
    // Then HfO2 spans z in [-tHf, 0], Si spans z in [-(tHf+tSi), -tHf].

    // HfO2 volume
    auto solidHf = new G4Box("HfO2", padXY/2, padXY/2, tHf/2);
    auto HfMat = G4Material::GetMaterial("HfO2");
    fLogicHfO2 = new G4LogicalVolume(solidHf, HfMat, "HfO2LV");
    const G4ThreeVector posHf(0, 0, -tHf/2);
    new G4PVPlacement(nullptr, posHf, fLogicHfO2, "HfO2PV", fLogicWorld, false, 0);

    // Si volume
    auto solidSi = new G4Box("Si", padXY/2, padXY/2, tSi/2);
    auto siMat = nist->FindOrBuildMaterial("G4_Si");
    fLogicSi = new G4LogicalVolume(solidSi, siMat, "SiLV");
    const G4ThreeVector posSi(0, 0, -(tHf + tSi/2));
    new G4PVPlacement(nullptr, posSi, fLogicSi, "SiPV", fLogicWorld, false, 0);

    // Setup scoring voxel grid bounds exactly matching HfO2 volume in world coords:
    // HfO2 box spans x,y in [-padXY/2, +padXY/2] and z in [-tHf, 0]
    const G4ThreeVector minCorner(-padXY/2, -padXY/2, -tHf);
    const G4ThreeVector maxCorner(+padXY/2, +padXY/2,    0.0);

    fGrid.Configure(minCorner, maxCorner, fVoxelDxNm*nm, fVoxelDyNm*nm, fVoxelDzNm*nm);
    fVacancy.ConfigureFromGrid(fGrid);

    // Regions & cuts
    SetupRegionsAndCuts();

    return physWorld;
}

void DetectorConstruction::SetupRegionsAndCuts() {
    // Assign HfO2 LV to its own region, allowing smaller cuts
    fHfO2Region = new G4Region("HfO2Region");
    fLogicHfO2->SetRegion(fHfO2Region);
    fHfO2Region->AddRootLogicalVolume(fLogicHfO2);

    auto cuts = new G4ProductionCuts();

    // Range cuts (tune later for performance vs fidelity).
    // Start with something moderate: 5 nm for electrons/positrons, 10 nm for gammas.
    cuts->SetProductionCut(10*nm, "gamma");
    cuts->SetProductionCut(5*nm,    "e-");
    cuts->SetProductionCut(5*nm,    "e+");

    fHfO2Region->SetProductionCuts(cuts);
}

void DetectorConstruction::ConstructSDandField() {
    // No sensitive detector volumes needed because we score via SteppingAction into VoxelGrid.
}
