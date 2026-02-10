#include "VacancyModel.hh"
#include "VoxelGrid.hh"

void VacancyModel::ConfigureFromGrid(const VoxelGrid& grid) {
    fNx = grid.Nx();
    fNy = grid.Ny();
    fNz = grid.Nz();

    const size_t n = (size_t)fNx * (size_t)fNy * (size_t)fNz;
    fIsVacancy.assign(n, 0);
    fEbank_eV.assign(n, 0.0f);

    auto seed = grid.GetSeedIndex();
    fSeedIx = seed.ix;
    fSeedIy = seed.iy;
    fSeedIz = seed.iz;
    fSeedFlat = Flatten(fSeedIx, fSeedIy, fSeedIz);

    ResetToSeedOnly();
}

void VacancyModel::ResetToSeedOnly() {
    std::fill(fIsVacancy.begin(), fIsVacancy.end(), 0);
    std::fill(fEbank_eV.begin(), fEbank_eV.end(), 0.0f);
    fIsVacancy[fSeedFlat] = 1;

    fSeedCapturedElectrons = 0;
    fTotalCreated = 0;
}

bool VacancyModel::IsInBounds(int ix, int iy, int iz) const {
    return (ix >= 0 && ix < fNx && iy >= 0 && iy < fNy && iz >= 0 && iz < fNz);
}

size_t VacancyModel::Flatten(int ix, int iy, int iz) const {
    return (size_t)iz + (size_t)fNz * ((size_t)iy + (size_t)fNy * (size_t)ix);
}

void VacancyModel::Unflatten(size_t flat, int& ix, int& iy, int& iz) const {
    const size_t yz = (size_t)fNy * (size_t)fNz;
    ix = (int)(flat / yz);
    const size_t rem = flat - (size_t)ix * yz;
    iy = (int)(rem / (size_t)fNz);
    iz = (int)(rem - (size_t)iy * (size_t)fNz);
}

bool VacancyModel::HasVacancyNeighbor6(int ix, int iy, int iz) const {
    // 6-neighborhood
    const int dx[6] = {+1,-1, 0, 0, 0, 0};
    const int dy[6] = { 0, 0,+1,-1, 0, 0};
    const int dz[6] = { 0, 0, 0, 0,+1,-1};

    for (int k=0;k<6;++k) {
        int nx = ix + dx[k], ny = iy + dy[k], nz = iz + dz[k];
        if (!IsInBounds(nx,ny,nz)) continue;
        if (fIsVacancy[Flatten(nx,ny,nz)]) return true;
    }
    return false;
}

bool VacancyModel::IsNeighborOfSeed6(int ix, int iy, int iz) const {
    const int md = std::abs(ix - fSeedIx) + std::abs(iy - fSeedIy) + std::abs(iz - fSeedIz);
    return (md == 1);
}

void VacancyModel::ProcessEvent(const VoxelGrid& grid) {
    // 1) Add event Edep to Ebank
    const auto& touched = grid.GetTouchedVoxels();
    for (size_t flat : touched) {
        const double edep_eV = grid.GetEdepEvent_eV(flat);
        if (edep_eV > 0.0) {
            fEbank_eV[flat] += (float)edep_eV;
        }
    }

    // 2) Update seed captured electrons from THIS event
    if (fSeedCapturedElectrons < 2) {
        const double edepSeed_eV = grid.GetEdepEvent_eV(fSeedFlat);
        if (edepSeed_eV > 0.0 && fP.W_eV > 0.0) {
            const int dn = (int)std::floor(edepSeed_eV / fP.W_eV);
            if (dn > 0) {
                fSeedCapturedElectrons = std::min(2, fSeedCapturedElectrons + dn);
            }
        }
    }

    // 3) Attempt to create vacancies in candidate voxels (touched + adjacent to existing vacancy)
    for (size_t flat : touched) {
        if (fIsVacancy[flat]) continue; // already vacancy

        int ix, iy, iz;
        Unflatten(flat, ix, iy, iz);

        if (!HasVacancyNeighbor6(ix,iy,iz)) continue;

        // Choose barrier
        double Ea = fP.Ea_base_eV;
        if (fSeedCapturedElectrons >= 2) {
            if (!fP.fastOnlyNearSeed) {
                Ea = fP.Ea_fast_eV;
            } else {
                if (IsNeighborOfSeed6(ix,iy,iz)) Ea = fP.Ea_fast_eV;
            }
        }

        if ((double)fEbank_eV[flat] >= Ea) {
            // Create vacancy
            fIsVacancy[flat] = 1;
            fEbank_eV[flat] = (float)((double)fEbank_eV[flat] - Ea);
            fTotalCreated += 1;
        }
    }
}

void VacancyModel::ExportVacancyCSV(const std::string& path, const VoxelGrid& grid) const {
    std::ofstream out(path);
    out << "ix,iy,iz,isVacancy,Ebank_eV,edepRun_eV,seed\n";

    for (int ix=0; ix<fNx; ++ix) {
        for (int iy=0; iy<fNy; ++iy) {
            for (int iz=0; iz<fNz; ++iz) {
                const size_t flat = Flatten(ix,iy,iz);
                out << ix << "," << iy << "," << iz << ","
                        << (int)fIsVacancy[flat] << ","
                        << (double)fEbank_eV[flat] << ","
                        << grid.GetEdepRun_eV(flat) << ","
                        << ((flat==fSeedFlat)?1:0) << "\n";
            }
        }
    }
}

void VacancyModel::ExportSummaryCSV(const std::string& path, long long nPrimaries) const {
    std::ofstream out(path);
    out << "key,value\n";
    out << "W_eV," << fP.W_eV << "\n";
    out << "Ea_base_eV," << fP.Ea_base_eV << "\n";
    out << "Ea_fast_eV," << fP.Ea_fast_eV << "\n";
    out << "seedCapturedElectrons," << fSeedCapturedElectrons << "\n";
    out << "totalCreated," << fTotalCreated << "\n";
    out << "nPrimaries," << nPrimaries << "\n";
    out << "createdPerPrimary," << (nPrimaries>0 ? (double)fTotalCreated/(double)nPrimaries : 0.0) << "\n";
}
