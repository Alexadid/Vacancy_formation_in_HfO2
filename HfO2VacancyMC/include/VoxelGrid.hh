#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

class VoxelGrid {
public:
    struct Index3 { int ix, iy, iz; };

    void Configure(const G4ThreeVector& minCorner,
                                 const G4ThreeVector& maxCorner,
                                 G4double dx, G4double dy, G4double dz)
    {
        fMin = minCorner;
        fMax = maxCorner;
        fDx = dx; fDy = dy; fDz = dz;

        const auto size = fMax - fMin;
        fNx = (int)std::ceil(size.x() / fDx);
        fNy = (int)std::ceil(size.y() / fDy);
        fNz = (int)std::ceil(size.z() / fDz);

        if (fNx <= 0 || fNy <= 0 || fNz <= 0) {
            throw std::runtime_error("VoxelGrid: invalid dimensions.");
        }

        const size_t n = (size_t)fNx * (size_t)fNy * (size_t)fNz;
        fEdepRun.assign(n, 0.0);
        fEdepEvent.assign(n, 0.0);
        fTouchedFlag.assign(n, 0);
        fHasSeedVacancy.assign(n, 0);

        SetSeedVacancyAtCenter();
    }

    // Event accumulators
    void ResetEventAccumulators() {
        for (size_t flat : fTouched) {
            fEdepEvent[flat] = 0.0;
            fTouchedFlag[flat] = 0;
        }
        fTouched.clear();
    }

    void ResetRunAccumulators() {
        std::fill(fEdepRun.begin(), fEdepRun.end(), 0.0);
    }

    inline bool Contains(const G4ThreeVector& p) const {
        return (p.x() >= fMin.x() && p.x() < fMax.x() &&
                        p.y() >= fMin.y() && p.y() < fMax.y() &&
                        p.z() >= fMin.z() && p.z() < fMax.z());
    }

    inline Index3 ToIndex(const G4ThreeVector& p) const {
        Index3 idx;
        idx.ix = (int)std::floor((p.x() - fMin.x()) / fDx);
        idx.iy = (int)std::floor((p.y() - fMin.y()) / fDy);
        idx.iz = (int)std::floor((p.z() - fMin.z()) / fDz);
        idx.ix = std::max(0, std::min(fNx - 1, idx.ix));
        idx.iy = std::max(0, std::min(fNy - 1, idx.iy));
        idx.iz = std::max(0, std::min(fNz - 1, idx.iz));
        return idx;
    }

    inline size_t Flatten(const Index3& idx) const {
        return (size_t)idx.iz + (size_t)fNz * ((size_t)idx.iy + (size_t)fNy * (size_t)idx.ix);
    }

    inline Index3 Unflatten(size_t flat) const {
        const size_t yz = (size_t)fNy * (size_t)fNz;
        Index3 idx;
        idx.ix = (int)(flat / yz);
        const size_t rem = flat - (size_t)idx.ix * yz;
        idx.iy = (int)(rem / (size_t)fNz);
        idx.iz = (int)(rem - (size_t)idx.iy * (size_t)fNz);
        return idx;
    }

    void AddEdep(const G4ThreeVector& p, G4double edep) {
        if (edep <= 0.0) return;
        if (!Contains(p)) return;

        const auto idx = ToIndex(p);
        const size_t flat = Flatten(idx);

        fEdepRun[flat] += edep;
        fEdepEvent[flat] += edep;

        if (!fTouchedFlag[flat]) {
            fTouchedFlag[flat] = 1;
            fTouched.push_back(flat);
        }
    }

    void SetSeedVacancyAtCenter() {
        if (fNx<=0 || fNy<=0 || fNz<=0) return;
        Index3 c{fNx/2, fNy/2, fNz/2};
        fSeed = c;
        std::fill(fHasSeedVacancy.begin(), fHasSeedVacancy.end(), 0);
        fHasSeedVacancy[Flatten(c)] = 1;
    }

    Index3 GetSeedIndex() const { return fSeed; }

    const std::vector<size_t>& GetTouchedVoxels() const { return fTouched; }

    double GetEdepEvent_eV(size_t flat) const { return fEdepEvent[flat] / eV; }
    double GetEdepRun_eV(size_t flat) const { return fEdepRun[flat] / eV; }

    void ExportEdepCSV(const std::string& path) const {
        std::ofstream out(path);
        out << "ix,iy,iz,edepRun_eV,seed\n";
        for (int ix=0; ix<fNx; ++ix) {
            for (int iy=0; iy<fNy; ++iy) {
                for (int iz=0; iz<fNz; ++iz) {
                    Index3 idx{ix,iy,iz};
                    const auto flat = Flatten(idx);
                    out << ix << "," << iy << "," << iz << ","
                            << GetEdepRun_eV(flat) << ","
                            << (int)fHasSeedVacancy[flat] << "\n";
                }
            }
        }
    }

    int Nx() const { return fNx; }
    int Ny() const { return fNy; }
    int Nz() const { return fNz; }
    G4double Dx() const { return fDx; }
    G4double Dy() const { return fDy; }
    G4double Dz() const { return fDz; }

    G4ThreeVector Min() const { return fMin; }
    G4ThreeVector Max() const { return fMax; }

private:
    G4ThreeVector fMin{0,0,0}, fMax{0,0,0};
    G4double fDx{50*nm}, fDy{50*nm}, fDz{1*nm};
    int fNx{0}, fNy{0}, fNz{0};

    std::vector<double> fEdepRun;   // Geant4 energy units
    std::vector<double> fEdepEvent; // per-event accum
    std::vector<uint8_t> fTouchedFlag;
    std::vector<size_t> fTouched;

    std::vector<uint8_t> fHasSeedVacancy;
    Index3 fSeed{0,0,0};
};
