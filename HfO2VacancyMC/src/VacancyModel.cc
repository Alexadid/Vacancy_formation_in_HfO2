#include "VacancyModel.hh"
#include "VoxelGrid.hh"
#include "G4SystemOfUnits.hh" // for cm

static constexpr double kNA = 6.02214076e23; // Avogadro (mol^-1)

double VacancyModel::OxygenSiteDensity_cm3() const {
    // n_O = 2 * (rho/M) * NA
    const double n_formula = (fP.rho_g_cm3 / fP.molarMass_g_mol) * kNA;
    return 2.0 * n_formula;
}

uint32_t VacancyModel::CapacityPerVoxel(const VoxelGrid& grid) const {
    const double nO = OxygenSiteDensity_cm3(); // cm^-3
    const double Vvox_cm3 =
            (grid.Dx() / cm) * (grid.Dy() / cm) * (grid.Dz() / cm);
    const double cap = std::floor(nO * Vvox_cm3);
    return (cap < 1.0) ? 1u : (uint32_t)cap;
}

void VacancyModel::ConfigureFromGrid(const VoxelGrid& grid) {
    fNx = grid.Nx();
    fNy = grid.Ny();
    fNz = grid.Nz();

    const size_t n = (size_t)fNx * (size_t)fNy * (size_t)fNz;
    fVacCount.assign(n, 0);
    fEbank_eV.assign(n, 0.0f);

    auto seed = grid.GetSeedIndex();
    fSeedIx = seed.ix;
    fSeedIy = seed.iy;
    fSeedIz = seed.iz;
    fSeedFlat = Flatten(fSeedIx, fSeedIy, fSeedIz);

    fCapPerVoxel = CapacityPerVoxel(grid);

    ResetAndInit(grid);
}

void VacancyModel::ResetAndInit(const VoxelGrid& grid) {
    std::fill(fVacCount.begin(), fVacCount.end(), 0);
    std::fill(fEbank_eV.begin(), fEbank_eV.end(), 0.0f);

    fSeedCapturedElectrons = 0;
    fTotalCreated = 0;

    // Init RNG
    fRng.seed(fP.initSeed);

    // Clamp concentration by physical maximum nO
    const double nO = OxygenSiteDensity_cm3();
    double C0 = std::max(0.0, fP.initConc_cm3);
    if (C0 > nO) C0 = nO;

    // Fill vacancies per voxel using Poisson(C0 * Vvox), then cap by capacity
    const double Vvox_cm3 =
            (grid.Dx() / cm) * (grid.Dy() / cm) * (grid.Dz() / cm);
    const double lambda = C0 * Vvox_cm3;

    // For speed/stability: Poisson is fine for your sizes; if needed, add a normal approx for huge lambda.
    std::poisson_distribution<int> pois(lambda);

    for (size_t i = 0; i < fVacCount.size(); ++i) {
        uint32_t v = 0;
        if (lambda > 0.0) {
            int draw = pois(fRng);
            if (draw < 0) draw = 0;
            v = (uint32_t)draw;
            if (v > fCapPerVoxel) v = fCapPerVoxel;
        }
        fVacCount[i] = v;
    }

    // Ensure at least one seed vacancy in the center voxel
    if (fVacCount[fSeedFlat] == 0) fVacCount[fSeedFlat] = 1;
}

// --- helpers (same as before, but now "vacancy exists" means vacCount>0)

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
    const int dx[6] = {+1,-1, 0, 0, 0, 0};
    const int dy[6] = { 0, 0,+1,-1, 0, 0};
    const int dz[6] = { 0, 0, 0, 0,+1,-1};
    for (int k=0;k<6;++k) {
        int nx = ix + dx[k], ny = iy + dy[k], nz = iz + dz[k];
        if (!IsInBounds(nx,ny,nz)) continue;
        if (fVacCount[Flatten(nx,ny,nz)] > 0) return true;
    }
    return false;
}

bool VacancyModel::IsNeighborOfSeed6(int ix, int iy, int iz) const {
    const int md = std::abs(ix - fSeedIx) + std::abs(iy - fSeedIy) + std::abs(iz - fSeedIz);
    return (md == 1);
}

void VacancyModel::ProcessEvent(const VoxelGrid& grid) {
    // 1) add event edep to energy bank
    const auto& touched = grid.GetTouchedVoxels();
    for (size_t flat : touched) {
        const double edep_eV = grid.GetEdepEvent_eV(flat);
        if (edep_eV > 0.0) fEbank_eV[flat] += (float)edep_eV;
    }

    // 2) update seed captured electrons
    if (fSeedCapturedElectrons < 2) {
        const double edepSeed_eV = grid.GetEdepEvent_eV(fSeedFlat);
        if (edepSeed_eV > 0.0 && fP.W_eV > 0.0) {
            const int dn = (int)std::floor(edepSeed_eV / fP.W_eV);
            if (dn > 0) fSeedCapturedElectrons = std::min(2, fSeedCapturedElectrons + dn);
        }
    }

    // 3) create new vacancies in touched voxels adjacent to existing vacancies
    for (size_t flat : touched) {
        // if voxel already "full" of vacancies, skip
        if (fVacCount[flat] >= fCapPerVoxel) continue;

        int ix, iy, iz;
        Unflatten(flat, ix, iy, iz);

        if (!HasVacancyNeighbor6(ix,iy,iz)) continue;

        double Ea = fP.Ea_base_eV;
        if (fSeedCapturedElectrons >= 2) {
            if (!fP.fastOnlyNearSeed) Ea = fP.Ea_fast_eV;
            else if (IsNeighborOfSeed6(ix,iy,iz)) Ea = fP.Ea_fast_eV;
        }

        if ((double)fEbank_eV[flat] >= Ea) {
            // Create ONE vacancy (you can allow multiple by while-loop if you want)
            fVacCount[flat] += 1;
            fEbank_eV[flat] = (float)((double)fEbank_eV[flat] - Ea);
            fTotalCreated += 1;
        }
    }
}

void VacancyModel::ExportVacancyCSV(const std::string& path, const VoxelGrid& grid) const {
    std::ofstream out(path);
    out << "ix,iy,iz,vacCount,Ebank_eV,edepRun_eV,seed\n";
    for (int ix=0; ix<fNx; ++ix) {
        for (int iy=0; iy<fNy; ++iy) {
            for (int iz=0; iz<fNz; ++iz) {
                const size_t flat = Flatten(ix,iy,iz);
                out << ix << "," << iy << "," << iz << ","
                        << fVacCount[flat] << ","
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
    out << "initConc_cm3," << fP.initConc_cm3 << "\n";
    out << "rho_g_cm3," << fP.rho_g_cm3 << "\n";
    out << "capPerVoxel," << fCapPerVoxel << "\n";
    out << "W_eV," << fP.W_eV << "\n";
    out << "Ea_base_eV," << fP.Ea_base_eV << "\n";
    out << "Ea_fast_eV," << fP.Ea_fast_eV << "\n";
    out << "seedCapturedElectrons," << fSeedCapturedElectrons << "\n";
    out << "totalCreated," << fTotalCreated << "\n";
    out << "nPrimaries," << nPrimaries << "\n";
    out << "createdPerPrimary," << (nPrimaries>0 ? (double)fTotalCreated/(double)nPrimaries : 0.0) << "\n";
}
