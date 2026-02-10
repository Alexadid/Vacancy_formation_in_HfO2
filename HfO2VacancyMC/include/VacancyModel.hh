#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <fstream>
#include <algorithm>
#include <random>
#include <cmath>

class VoxelGrid;

class VacancyModel {
public:
    struct Params {
        // stage-1 parameters:
        double W_eV             = 15.0;
        double Ea_base_eV       = 2.0;
        double Ea_fast_eV       = 1.3;
        bool   fastOnlyNearSeed = true;

        // stage-2 parameters:
        double initConc_cm3     = 0.0;      // initial vacancy concentration (cm^-3)
        uint64_t initSeed       = 12345;    // RNG seed for initial distribution

        // Material parameters for capacity calculation:
        double rho_g_cm3        = 9.68;     // HfO2 density (adjustable)
        double molarMass_g_mol  = 210.49;   // HfO2 molar mass
    };

        VacancyModel() = default;

    void ConfigureFromGrid(const VoxelGrid& grid);
    void ResetAndInit(const VoxelGrid& grid);    // uses Params.initConc_cm3

    void ProcessEvent(const VoxelGrid& grid);

    // Getters
    long long TotalCreated() const { return fTotalCreated; }
    int SeedCapturedElectrons() const { return fSeedCapturedElectrons; }
    const Params& GetParams() const { return fP; }
    Params& GetParams() { return fP; }

    void ExportVacancyCSV(const std::string& path, const VoxelGrid& grid) const;
    void ExportSummaryCSV(const std::string& path, long long nPrimaries) const;

private:
    // internal helpers
    bool IsInBounds(int ix, int iy, int iz) const;
    size_t Flatten(int ix, int iy, int iz) const;
    void Unflatten(size_t flat, int& ix, int& iy, int& iz) const;

    bool HasVacancyNeighbor6(int ix, int iy, int iz) const;
    bool IsNeighborOfSeed6(int ix, int iy, int iz) const;

    double OxygenSiteDensity_cm3() const;
    uint32_t CapacityPerVoxel(const VoxelGrid& grid) const;

private:
    Params fP;

    int fNx{0}, fNy{0}, fNz{0};
    int fSeedIx{0}, fSeedIy{0}, fSeedIz{0};
    size_t fSeedFlat{0};

    // stage-2 state:
    std::vector<uint32_t> fVacCount; // number of vacancies in voxel (0..cap)
    uint32_t fCapPerVoxel{0};

    // energy bank:
    std::vector<float> fEbank_eV;

    int fSeedCapturedElectrons{0}; // 0..2
    long long fTotalCreated{0};

    std::mt19937_64 fRng;
};
