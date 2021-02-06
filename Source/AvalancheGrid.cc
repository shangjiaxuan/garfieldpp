#include "Garfield/AvalancheGrid.hh"

#include <TF1.h>

#include <algorithm>
#include <cmath>
#include <iostream>

#include "Garfield/Medium.hh"
#include "Garfield/Random.hh"

namespace Garfield {
void AvalancheGrid::SetZGrid(Grid& av, const double ztop, const double zbottom,
                             const int zsteps) {
  // Creating the z-coordinate grid.
  av.zsteps = zsteps;
  av.zStepSize = (ztop - zbottom) / zsteps;
  // Loop bassed grid creation.
  for (int i = 0; i < zsteps; i++) {
    av.zgrid.push_back(zbottom + i * av.zStepSize);
  }
}

void AvalancheGrid::SetYGrid(Grid& av, const double ytop, const double ybottom,
                             const int ysteps) {
  // Idem to SetZGrid for the x-coordinate grid.
  av.ysteps = ysteps;
  av.yStepSize = (ytop - ybottom) / ysteps;

  if (av.yStepSize == 0) av.yStepSize = 1;

  for (int i = 0; i < ysteps; i++) {
    av.ygrid.push_back(ybottom + i * av.yStepSize);
  }
}

void AvalancheGrid::SetXGrid(Grid& av, const double xtop, const double xbottom,
                             const int xsteps) {
  // Idem to SetZGrid for the x-coordinate grid.
  av.xsteps = xsteps;
  av.xStepSize = (xtop - xbottom) / xsteps;

  if (av.xStepSize == 0) av.xStepSize = 1;

  for (int i = 0; i < xsteps; i++) {
    av.xgrid.push_back(xbottom + i * av.xStepSize);
  }

  std::vector<int> nhx(xsteps, 0);
  std::vector<std::vector<int>> nhy(av.ysteps, nhx);
  std::vector<std::vector<std::vector<int>>> nhz(av.zsteps, nhy);
  av.n = nhz;

  // Get the diffusion factors for neighboring points on the grid.
  DiffusionFactors(av);
}

void AvalancheGrid::SetGrid(const double xmin, const double xmax,
                            const int xsteps, const double ymin,
                            const double ymax, const int ysteps,
                            const double zmin, const double zmax,
                            const int zsteps) {
  m_avgrid.gridset = true;

  if (zmin >= zmax || zsteps <= 0 || xmin > xmax || xsteps <= 0 ||
      ymin > ymax || ysteps <= 0) {
    std::cerr << m_className
              << "::SetGrid:Error: Grid is not properly defined.\n";
    return;
  }

  // Setting grid

  SetZGrid(m_avgrid, zmax, zmin, zsteps);
  SetYGrid(m_avgrid, ymax, ymin, ysteps);
  SetXGrid(m_avgrid, xmax, xmin, xsteps);

  if (m_sensor) GetParametersFromSensor();
}

int AvalancheGrid::GetAvalancheSize(double dx, const int nsize,
                                    const double alpha, const double eta) {
  // Algorithm to get the size of the avalanche after it has propagated over a
  // distance dx.

  int newnsize = 0;  // Holder for final size.

  const double k = eta / alpha;
  const double ndx = exp((alpha - eta) *
                         dx);  // Scaling Townsend and Attachment coef. to 1/mm.
  // If the size is higher than 1e3 the central limit theorem will be used to
  // describe the growth of the Townsend avalanche.
  if (nsize < 1e3) {
    // Running over all electrons in the avalanche.
    for (int i = 0; i < nsize; i++) {
      // Draw a random number from the uniform distribution (0,1).
      double s = RndmUniformPos();
      // Condition to which the random number will be compared. If the number is
      // smaller than the condition, nothing happens. Otherwise, the single
      // electron will be attached or retrieve additional electrons from the
      // gas.
      double condition = k * (ndx - 1) / (ndx - k);

      if (s >= condition)
        newnsize += (int)(1 + log((ndx - k) * (1 - s) / (ndx * (1 - k))) /
                                  log(1 - (1 - k) / (ndx - k)));
    }

  } else {
    // Central limit theorem.
    const double sigma = sqrt((1 + k) * ndx * (ndx - 1) / (1 - k));
    newnsize = RndmGaussian(nsize * ndx, sqrt(nsize) * sigma);
  }

  return newnsize;
}

bool AvalancheGrid::SnapToGrid(Grid& av, const double x, const double y,
                               const double z, const double v, const int n) {
  // Snap electron from AvalancheMicroscopic to the predefined grid.
  if (!av.gridset) {
    std::cerr << m_className << "::SnapToGrid:Error: grid is not defined.\n";
    return false;
  }
  // Finding the z position on the grid.

  int indexX, indexY, indexZ = 0;

  if (m_velNormal[0] != 0) {
    indexX = m_velNormal[0] < 0 ? floor((x - av.xgrid.front()) / av.xStepSize)
                                : ceil((x - av.xgrid.front()) / av.xStepSize);
    indexY = round((y - av.ygrid.front()) / av.yStepSize);
    indexZ = round((z - av.zgrid.front()) / av.zStepSize);
  } else if (m_velNormal[1] != 0) {
    indexX = round((x - av.xgrid.front()) / av.xStepSize);
    indexY = m_velNormal[1] < 0 ? floor((y - av.ygrid.front()) / av.yStepSize)
                                : ceil((y - av.ygrid.front()) / av.yStepSize);
    indexZ = round((z - av.zgrid.front()) / av.zStepSize);
  } else {
    indexX = round((x - av.xgrid.front()) / av.xStepSize);
    indexY = round((y - av.ygrid.front()) / av.yStepSize);
    indexZ = m_velNormal[2] < 0 ? floor((z - av.zgrid.front()) / av.zStepSize)
                                : ceil((z - av.zgrid.front()) / av.zStepSize);
  }

  if (indexX < 0 || indexX >= av.xsteps || indexY < 0 || indexY >= av.ysteps ||
      indexZ < 0 || indexZ >= av.zsteps)
    return false;

  av.gridPosition[2].push_back(indexX);
  av.gridPosition[1].push_back(indexY);
  av.gridPosition[0].push_back(indexZ);

  // When snapping the electron to the grid the distance traveled can yield
  // additional electrons or attachment.

  double step = z - av.zgrid[indexZ];

  if (m_velNormal[0] != 0) {
    step = x - av.xgrid[indexX];
  } else if (m_velNormal[1] != 0) {
    step = y - av.ygrid[indexY];
  }

  const int nholder = GetAvalancheSize(step, n, m_Townsend, m_Attachment);

  // av.N += nholder;
  av.N += n;

  av.n[indexZ][indexY][indexX] += nholder;
  if (m_debug)
    std::cerr << m_className << "::SnapToGrid: n from 1 to " << nholder
              << ".\n";

  if (m_debug)
    std::cerr << m_className << "::SnapToGrid: Snapped to (x,y,z) = (" << x
              << " -> " << av.xgrid[indexX] << ", " << y << " -> "
              << av.ygrid[indexY] << ", " << z << " -> " << av.zgrid[indexZ]
              << ").\n";
  return true;
}

void AvalancheGrid::NextAvalancheGridPoint(Grid& av) {
  // This main function propagates the electrons and applies the avalanche
  // statistics.
  int Nholder = 0;  // Holds the avalanche size before propagating it to the
  // next point in the grid.

  av.run = false;

  for (int& iz : av.gridPosition[0]) {  // For every avalanche position on the
    // z-coorindate grid.

    // Check if it has reached the bottom resistive plate.

    if ((m_velNormal[2] < 0 && iz == 0) ||
        (m_velNormal[2] > 0 && iz == av.zsteps - 1)) {
      continue;
    } else {
      av.run = true;
    }

    for (int& iy : av.gridPosition[1]) {
      if ((m_velNormal[1] < 0 && iy == 0) ||
          (m_velNormal[1] > 0 && iy == av.ysteps - 1)) {
        av.run = false;
        continue;
      } else {
        av.run = true;
      }

      for (int& ix : av.gridPosition[2]) {
        if ((m_velNormal[0] < 0 && ix == 0) ||
            (m_velNormal[0] > 0 && ix == av.xsteps - 1)) {
          av.run = false;
          continue;
        } else {
          av.run = true;
        }

        // Get avalanche size at z-index= iz and x-index=ix.
        Nholder = av.n[iz][iy][ix];

        if (Nholder == 0) continue;  // If empty go to next point.
        if (m_diffusion) {
          // Idem to the else part of this function, but with the additional
          // step that after the new avalanche size is calculated the charges
          // will be spread over the neighboring x-points following the normal
          // distribution given by the transverse diffusion coefficient.
          double holdnsize = 0.;
          if (av.N < m_MaxSize) {
            double step = av.zStepSize;

            if (m_velNormal[0] != 0) {
              step = av.xStepSize;
            } else if (m_velNormal[1] != 0) {
              step = av.yStepSize;
            }

            holdnsize = GetAvalancheSize(step, av.n[iz][iy][ix], m_Townsend,
                                         m_Attachment);

            if (m_MaxSize - av.N < holdnsize - av.n[iz][iy][ix])
              holdnsize = m_MaxSize - av.N + av.n[iz][iy][ix];

          } else {
            holdnsize = av.n[iz][iy][ix];
            m_Saturated = true;

            if (m_SaturationTime == -1)
              m_SaturationTime = av.time + std::abs(av.zStepSize / av.velocity);
          }

          int chargeRemaining =
              holdnsize;  // Keeps charge of the amount of charge that is left
          // over after it has spread to neighboring points. This
          // will be the amount that will stay at the same
          // x-coordinate grid point. This way no charge is lost
          // during the diffusion step.

          for (int j = av.transverseDiffusion.size() - 1; j >= 0; j--) {
            int jx = 0;
            int jy = 0;
            int jz = 0;

            if (m_velNormal[2] != 0) {
              jy = j;
              jz = j;
            } else if (m_velNormal[1] != 0) {
              jx = j;
              jz = j;
            } else {
              jx = j;
              jy = j;
            }

            if (ix - jx < 0 || ix + jx > av.xsteps - 1)
              continue;  // If in grid.
            if (iy - jy < 0 || iy + jy > av.ysteps - 1) continue;
            if (iz - jz < 0 || iz + jz > av.xsteps - 1) continue;

            if (j > 0) {  // For all neighboring points

              int nxd = (int)(av.transverseDiffusion[j] * holdnsize);

              av.n[iz - 1][iy][ix - j] += nxd;

              av.n[iz - 1][iy][ix + j] += nxd;

              m_sensor->AddSignal(-(nxd + Nholder) / 2, av.time,
                                  av.time + av.zStepSize / av.velocity,
                                  av.xgrid[ix], av.ygrid[iy], av.zgrid[iz],
                                  av.xgrid[ix - j], av.ygrid[iy],
                                  av.zgrid[iz - 1], false, true);

              m_sensor->AddSignal(-(nxd + Nholder) / 2, av.time,
                                  av.time + av.zStepSize / av.velocity,
                                  av.xgrid[ix], av.ygrid[iy], av.zgrid[iz],
                                  av.xgrid[ix + j], av.ygrid[iy],
                                  av.zgrid[iz - 1], false, true);
              chargeRemaining -= 2 * nxd;

            } else {  // For the initial x-position.
              av.n[iz - 1][iy][ix] += chargeRemaining;

              m_sensor->AddSignal(-(chargeRemaining + Nholder) / 2, av.time,
                                  av.time + av.zStepSize / av.velocity,
                                  av.xgrid[ix], av.ygrid[iy], av.zgrid[iz],
                                  av.xgrid[ix], av.ygrid[iy], av.zgrid[iz - 1],
                                  false, true);
            }
          }
          av.N += holdnsize - Nholder;

        } else {
          // If the total avalanche size is smaller than the set saturation
          // limit the GetAvalancheSize function is utilized to obtain the size
          // after its propagation to the next z-coordinate grid point. Else,
          // the size will be kept constant under the propagation.
          if (av.N < m_MaxSize) {
            int holdnsize = GetAvalancheSize(av.zStepSize, av.n[iz][iy][ix],
                                             m_Townsend, m_Attachment);

            if (m_MaxSize - av.N < holdnsize - av.n[iz][iy][ix])
              holdnsize = m_MaxSize - av.N + av.n[iz][iy][ix];

            av.n[iz + m_velNormal[2]][iy + m_velNormal[1]]
                [ix + m_velNormal[0]] = holdnsize;

          } else {
            av.n[iz + m_velNormal[2]][iy + m_velNormal[1]]
                [ix + m_velNormal[0]] = av.n[iz][iy][ix];
            m_Saturated = true;

            if (m_SaturationTime == -1)
              m_SaturationTime = av.time + std::abs(av.zStepSize / av.velocity);
          }
          // Produce induced signal on readout electrodes.

          m_sensor->AddSignal(
              -(av.n[iz][iy][ix] +
                av.n[iz + m_velNormal[2]][iy + m_velNormal[1]]
                    [ix + m_velNormal[0]]) /
                  2,
              av.time, av.time + av.zStepSize / av.velocity, av.xgrid[ix],
              av.ygrid[iy], av.zgrid[iz], av.xgrid[ix + m_velNormal[0]],
              av.ygrid[iy + m_velNormal[1]], av.zgrid[iz - 1], false, true);

          av.N += av.n[iz + m_velNormal[2]][iy + m_velNormal[1]]
                      [ix + m_velNormal[0]] -
                  Nholder;  // Update total number of electrons.
        }

        av.n[iz][iy][ix] = 0;  // Clear previous z-coordinate grid point.
      }
    }
  }

  // Update position index.
  if (m_velNormal[2] != 0) {
    for (int& iz : av.gridPosition[0])
      if ((m_velNormal[2] < 0 && iz != 0) ||
          (m_velNormal[2] > 0 && iz != av.zsteps - 1))
        iz += m_velNormal[2];
  } else if (m_velNormal[1] != 0) {
    for (int& iy : av.gridPosition[1])
      if ((m_velNormal[1] < 0 && iy != 0) ||
          (m_velNormal[1] > 0 && iy != av.ysteps - 1))
        iy += m_velNormal[1];
  } else {
    for (int& ix : av.gridPosition[2])
      if ((m_velNormal[0] < 0 && ix != 0) ||
          (m_velNormal[0] > 0 && ix != av.xsteps - 1))
        ix += m_velNormal[0];
  }
  // After all active grid points have propagated, update the time.

  av.time += std::abs(av.zStepSize / av.velocity);
}

void AvalancheGrid::StartGridAvalanche() {
  // Start the AvalancheGrid algorithm.
  if ((!m_avmc && !m_driftAvalanche) || !m_sensor) return;

  if (!m_importAvalanche && m_avmc) GetElectronsFromAvalancheMicroscopic();

  std::cerr << m_className
            << "::StartGridAvalanche::Starting grid based simulation with "
            << m_avgrid.N << " initial electrons.\n";
  if (m_avgrid.N <= 0) {
    std::cerr << m_className << "::StartGridAvalanche::Cancelled.\n";
    return;
  }

  m_nestart = m_avgrid.N;

  // The vector containing the indexes of the z-coordinate grid points of the
  // initial electrons needs to be ordered from small to large values. All
  // duplicate values need to be removed.

  GetParametersFromSensor();

  SortPositionVector();

  if (m_debug) {
    std::cerr << m_className
              << "::StartGridAvalanche::m_avgrid.gridPosition at iz = ";
    for (int i = 0; i < m_avgrid.gridPosition[0].size(); i++) {
      std::cerr << m_avgrid.gridPosition[0][i] << ",";
    }
    std::cerr << ".\n";
  }

  // Set velocity if given.
  m_avgrid.velocity = m_Velocity;
  // Main loop.
  while (m_avgrid.run == true) {
    NextAvalancheGridPoint(m_avgrid);
  }

  if (m_Saturated)
    std::cerr << m_className
              << "::StartGridAvalanche::Avalanche maximum size of " << m_MaxSize
              << " electrons reached at " << m_SaturationTime << " ns.\n";

  std::cerr << m_className
            << "::StartGridAvalanche::Final avalanche size = " << m_avgrid.N
            << " at t = " << m_avgrid.time << " ns.\n";

  return;
}

void AvalancheGrid::DiffusionFactors(Grid& av) {
  if (!m_diffusion) {
    av.transverseDiffusion.push_back(1);
    return;
  }

  // Get transverse diffusion factors, yielding in the spreading of charge.
  if (!av.gridset || av.xStepSize <= 0) return;

  auto cdfunctop = TF1("cdftop", "ROOT::Math::normal_cdf(x, [0],[1])", -5, 5);

  cdfunctop.SetParameters(m_DiffSigma, 0.0);

  double factor = 1;
  int index = 0;

  double x1, x2;
  while (factor > 1e-3) {  // 1e-3 is the precision cutoff.

    x1 = -av.xStepSize / 2 + index * av.xStepSize;
    x2 = av.xStepSize / 2 + index * av.xStepSize;

    factor = (std::erf(x2 / (m_DiffSigma * std::sqrt(2))) -
              std::erf(x1 / (m_DiffSigma * std::sqrt(2)))) /
             2;

    if (m_debug)
      std::cerr << m_className
                << "::DiffusionFactors::Transvers diffusion factor: " << factor
                << ", comparison: "
                << cdfunctop.Eval(0 - av.xStepSize / 2 + index * av.xStepSize) -
                       cdfunctop.Eval(0 + av.xStepSize / 2 +
                                      index * av.xStepSize)
                << ", top: "
                << cdfunctop.Eval(0 - av.xStepSize / 2 + index * av.xStepSize)
                << ", bottom: "
                << cdfunctop.Eval(0 + av.xStepSize / 2 + index * av.xStepSize)
                << ".\n";
    av.transverseDiffusion.push_back(factor);

    index++;
  }

  std::cerr << m_className
            << "::DiffusionFactors::Transvers diffusion spreads to "
            << av.transverseDiffusion.size() << " points.\n";

  return;
}

void AvalancheGrid::CreateAvalanche(const double x, const double y,
                                    const double z, const double t,
                                    const int n) {
  m_driftAvalanche = true;

  if (m_avgrid.time == 0 && m_avgrid.time != t && m_debug)
    std::cerr << m_className
              << "::CreateAvalanche::Overwriting start time of avalanche for t "
                 "= 0 to "
              << t << ".\n";
  m_avgrid.time = t;

  if (SnapToGrid(m_avgrid, x, y, z, 0, n) && m_debug)
    std::cerr << m_className
              << "::CreateAvalanche::Electron added at (t,x,y,z) =  (" << t
              << "," << x << "," << y << "," << z << ").\n";

  // std::cerr<< m_className<< "::CreateAvalanche::expected contribution is "<<
  // exp((m_Townsend-m_Attachment)*z) << ".\n";
}
void AvalancheGrid::GetElectronsFromAvalancheMicroscopic() {
  // Get the information of the electrons from the AvalancheMicroscopic class.
  if (!m_avmc) return;

  if (!m_importAvalanche) m_importAvalanche = true;

  int np = m_avmc->GetNumberOfElectronEndpoints();

  if (np == 0) return;

  // Get initial positions of electrons
  double x1, y1, z1, t1;
  double x2, y2, z2, t2;
  double e1, e2;
  int status;

  double vel = 0.;

  for (int i = 0; i < np; ++i) {
    m_avmc->GetElectronEndpoint(i, x1, y1, z1, t1, e1, x2, y2, z2, t2, e2,
                                status);

    vel = (z2 - z1) / (t2 - t1);

    m_avgrid.time = t2;

    if (SnapToGrid(m_avgrid, x2, y2, z2, vel) && m_debug)
      std::cerr << m_className
                << "::GetElectronsFromAvalancheMicroscopic::Electron added at "
                   "(x,y,z) =  ("
                << x2 << "," << y2 << "," << z2 << ").\n";
  }
}

void AvalancheGrid::GetParametersFromSensor() {
  if (!m_sensor || m_SensorParameters) return;

  double e[3], v;
  int status;
  Medium* m = nullptr;

  m_sensor->ElectricField(
      m_avgrid.xgrid[m_avgrid.xsteps / 2], m_avgrid.ygrid[m_avgrid.ysteps / 2],
      m_avgrid.zgrid[m_avgrid.zsteps / 2], e[0], e[1], e[2], v, m, status);

  if (m_Townsend == -1)
    m->ElectronTownsend(e[0], e[1], e[2], 0., 0., 0., m_Townsend);

  if (m_Attachment == -1)
    m->ElectronAttachment(e[0], e[1], e[2], 0., 0., 0., m_Attachment);

  if (m_Velocity == 0) {
    double vx, vy, vz;
    m->ElectronVelocity(e[0], e[1], e[2], 0., 0., 0., vx, vy, vz);

    double vel = sqrt(vx * vx + vy * vy + vz * vz);
    if (vel != std::abs(vx) && vel != std::abs(vy) && vel != std::abs(vz)) return;
    int nx = (int)round(vx / vel);
    int ny = (int)round(vy / vel);
    int nz = (int)round(vz / vel);
    m_velNormal = {nx, ny, nz};
    m_Velocity = -std::abs(vel);
  }

  std::cerr << m_className << "::GetParametersFromSensor::Electric field = ("
            << e[0] / 1000 << ", " << e[1] / 1000 << ", " << e[2] / 1000
            << ") [kV/cm].\n";

  std::cerr << m_className
            << "::GetParametersFromSensor::Townsend = " << m_Townsend
            << " [1/cm], Attachment = " << m_Attachment
            << " [1/cm], Velocity = " << m_Velocity << " [cm/ns].\n";

  m_SensorParameters = true;
}

void AvalancheGrid::SortPositionVector() {
  for (int i = 0; i < 3; i++) {
    sort(m_avgrid.gridPosition[i].begin(), m_avgrid.gridPosition[i].end());
    m_avgrid.gridPosition[i].erase(unique(m_avgrid.gridPosition[i].begin(),
                                          m_avgrid.gridPosition[i].end()),
                                   m_avgrid.gridPosition[i].end());
  }
}

void AvalancheGrid::Reset() {
  std::cerr << m_className << "::Reset::Resetting AvalancheGrid.\n";

  m_avgrid.n.clear();
  m_avgrid.transverseDiffusion.clear();
  m_avgrid.time = 0;
  m_avgrid.N = 0;
  m_avgrid.run = true;

  m_Saturated = false;
  m_SaturationTime = -1;

  m_driftAvalanche = false;

  std::vector<int> nhx(m_avgrid.xsteps, 0);
  std::vector<std::vector<int>> nhy(m_avgrid.ysteps, nhx);
  std::vector<std::vector<std::vector<int>>> nhz(m_avgrid.zsteps, nhy);
  m_avgrid.n = nhz;

  for (int i = 0; i < 3; i++) {
    m_avgrid.gridPosition[i].clear();
  }

  // Get the diffusion factors for neighboring points on the grid.
  DiffusionFactors(m_avgrid);
}
}  // namespace Garfield
