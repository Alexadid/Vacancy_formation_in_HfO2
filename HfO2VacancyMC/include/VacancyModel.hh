#pragma once
#include <vector>
#include <string>
#include <cstdint>
#include <fstream>
#include <algorithm>
#include <cmath>

class VoxelGrid;

class VacancyModel {
public:
    struct Params {
        double W_eV             = 15.0; // energy per excitation (start value)
        double Ea_base_eV       = 2.0;  // barrier without 2e- on seed
        double Ea_fast_eV       = 1.3;  // barrier near charged seed
        bool fastOnlyNearSeed   = true; // v1: only neighbors of seed get fast barrier
    };

    VacancyModel() = default;  

    void ConfigureFromGrid(const VoxelGrid& grid);
    void ResetToSeedOnly();

    void ProcessEvent(const VoxelGrid& grid);

    // Getters
    long long TotalCreated() const { return fTotalCreated; }
    int SeedCapturedElectrons() const { return fSeedCapturedElectrons; }
    const Params& GetParams() const { return fP; }
    Params& GetParams() { return fP; }

    // Export
    void ExportVacancyCSV(const std::string& path, const VoxelGrid& grid) const;
    void ExportSummaryCSV(const std::string& path, long long nPrimaries) const;

private:
    // internal helpers
    bool IsInBounds(int ix, int iy, int iz) const;
    size_t Flatten(int ix, int iy, int iz) const;
    void Unflatten(size_t flat, int& ix, int& iy, int& iz) const;

    bool HasVacancyNeighbor6(int ix, int iy, int iz) const;
    bool IsNeighborOfSeed6(int ix, int iy, int iz) const;

private:
    Params fP;

    int fNx{0}, fNy{0}, fNz{0};
    int fSeedIx{0}, fSeedIy{0}, fSeedIz{0};
    size_t fSeedFlat{0};

    // State arrays (size Nx*Ny*Nz)
    std::vector<uint8_t> fIsVacancy;    // 0/1
    std::vector<float>     fEbank_eV;   // accumulated energy in eV
    
    int fSeedCapturedElectrons{0};      // 0..2
    long long fTotalCreated{0};
};
