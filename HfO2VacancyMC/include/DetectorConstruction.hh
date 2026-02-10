#pragma once
#include "G4VUserDetectorConstruction.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4Region.hh"
#include "VoxelGrid.hh"

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
    DetectorConstruction();
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // Access to voxel grid for scoring/postprocessing
    VoxelGrid& GetVoxelGrid() { return fGrid; }
    const VoxelGrid& GetVoxelGrid() const { return fGrid; }

    // Geometry params
    double GetHfO2ThicknessNm() const { return fHfO2ThicknessNm; }

private:
    void DefineMaterials();
    void SetupRegionsAndCuts();

    // Parameters (settable by UI)
    double fHfO2ThicknessNm = 10.0;     // d (nm)
    double fPadSizeUm = 5.0;            // 5 x 5 um pad
    double fSiThicknessUm = 5.0;        // 5 um silicon thickness

    // Voxel grid resolution (settable)
    double fVoxelDxNm = 50.0;
    double fVoxelDyNm = 50.0;
    double fVoxelDzNm = 1.0;

    G4GenericMessenger* fMessenger = nullptr;

    // Pointers to volumes
    G4LogicalVolume* fLogicWorld = nullptr;
    G4LogicalVolume* fLogicSi = nullptr;
    G4LogicalVolume* fLogicHfO2 = nullptr;

    // Region for HfO2 (cuts)
    G4Region* fHfO2Region = nullptr;

    // Voxel grid in HfO2 volume
    VoxelGrid fGrid;
};
