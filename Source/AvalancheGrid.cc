#include "Garfield/AvalancheGrid.hh"

#include <TF1.h>

#include <algorithm>
#include <cmath>
#include <iostream>

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
  // If diffusion is not enabled the z-coordinate grid is the only grid needed.
  // So all the additional parameters will be set accordingly.
  if (!m_diffusion) {
    // x-coordininate grid is only one point, x=0 in this case.
    av.xgrid.push_back(0.);
    std::vector<std::vector<int>> nh(zsteps, {0});
    av.n = nh;
    av.xsteps = 1;
    // The electrons will not defuse and thus stay on the same x-coordinate.
    av.transverseDiffusion = {1.};
  }

  av.gridset = true;
}

void AvalancheGrid::SetXGrid(Grid& av, const double xtop, const double xbottom,
                             const int xsteps) {
  // Idem to SetZGrid for the x-coordinate grid.
  av.xsteps = xsteps;
  av.xStepSize = (xtop - xbottom) / xsteps;

  for (int i = 0; i < xsteps; i++) {
    av.xgrid.push_back(xbottom + i * av.xStepSize);
  }

  std::vector<int> nh0(xsteps, 0);
  std::vector<std::vector<int>> nh(av.zsteps, nh0);
  av.n = nh;
  // Get the diffusion factors for neighboring points on the grid.
  DiffusionFactors(av);
}

int AvalancheGrid::GetAvalancheSize(double dx, const int nsize,
                                    const double alpha, const double eta) {
  // Algorithm to get the size of the avalanche after it has propagated over a
  // distance dx.
  dx = 10 * dx;  // Scaling to get the dimensions right [mm].

  int newnsize = 0;  // Holder for final size.

  const double k = eta / alpha;
  const double ndx = exp((alpha - eta) * dx);
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
      if (s >= condition) {
        newnsize += (int)(1 + log((ndx - k) * (1 - s) / (ndx * (1 - k))) /
                                log(1 - (1 - k) / (ndx - k)));
      }
    }

  } else {
    // Central limit theorem.
    const double sigma = sqrt((1 + k) * ndx * (ndx - 1) / (1 - k));
    newnsize = RndmGaussian(nsize * ndx, sqrt(nsize) * sigma);
  }

  return newnsize;
}

void AvalancheGrid::SnapToGrid(Grid& av, const double x, const double z,
                               const double v) {
  // Snap electron from AvalancheMicroscopic to the predefined grid.
  if (av.gridset == false) {
    std::cerr << m_className << "::SnapToGrid:Error: grid is not defined.\n";
    return;
  }
  // Get electron velocity from initial electrons, but this can be overwritten
  // using the SetElectronVelocity function.
  av.velocity = (av.velocity * av.N + v) / (av.N + 1);
  // Finding the z position on the grid.
  int index = av.zgrid.size() - 1;

  for (int i = 0; i < av.zsteps; i++) {
    if (z >= av.zgrid[index]) {
      break;

    } else {
      index -= 1;
    }

    if (index < 0) {
      index = 0;
      break;
    }
  }

  av.gridPosition.push_back(index);
  // When snapping the electron to the grid the distance traveled can yield
  // additional electrons or attachment.
  const double nholder = GetAvalancheSize(z - av.zgrid[index], 1, 13, 3.5);
  av.N += nholder;
  // std::cerr << m_className << "::SnapToGrid::size after snap "<<nholder<<",
  // for z = "<<z<<", and z[i] = "<<av.zgrid[index]<<".\n";
  // We need to do the same for the x position if diffusion is on. The
  // convention is the following: av.n[<index of z position>][<index of x
  // position>].
  if (!m_diffusion) {
    av.n[index][0] += nholder;

  } else {
    int index2 = 0;

    for (int i = 0; i < av.xsteps; i++) {
      if (x <= av.xgrid[i]) {
        index2 = i;

        if (i != 0 && (x - av.xgrid[i] < -x + av.xgrid[i - 1])) index2 -= 1;

        break;
      }

      if (i == av.xsteps - 1) {
        index2 = i;
      }
    }

    av.n[index][index2] += nholder;
  }
}

void AvalancheGrid::NextAvalancheGridPoint(Grid& av) {
  // This main function propagates the electrons and applies the avalanche
  // statistics.
  int Nholder = 0;  // Holds the avalanche size before propagating it to the
                    // next point in the grid.

  av.run = false;

  for (int& pos : av.gridPosition) {  // For every avalanche position on the
                                      // z-coorindate grid.

    if (pos == 0)
      continue;  // Check if it hqs reached the bottom resistive plate.

    av.run = true;

    for (int i = 0; i < av.xsteps; i++) {
      // Get avalanche size at z-index= pos and x-index=i.
      Nholder = av.n[pos][i];

      if (Nholder == 0) continue;  // If empty go to next point.

      if (m_diffusion) {
        // Idem to the else part of this function, but with the additional step
        // that after the new avalanche size is calculated the charges will be
        // spread over the neighboring x-points following the normal
        // distribution given by the transverse diffusion coefficient.
        double holdnsize = 0.;

        if (av.N < m_MaxSize) {
          holdnsize = GetAvalancheSize(av.zStepSize, av.n[pos][i], m_Townsend,
                                       m_Attachment);

        } else {
          holdnsize = av.n[pos][i];
          m_Saturated = true;

          if (m_SaturationTime == -1)
            m_SaturationTime = av.time + abs(av.zStepSize / av.velocity);
        }

        int chargeRemaining =
            holdnsize;  // Keeps charge of the amount of charge that is left
                        // over after it has spread to neighboring points. This
                        // will be the amount that will stay at the same
                        // x-coordinate grid point. This way no charge is lost
                        // during the diffusion step.

        for (int j = av.transverseDiffusion.size() - 1; j >= 0; j--) {
          if (i - j < 0 || i + j > av.xsteps) continue;  // If in grid.

          if (j > 0) {  // For all neighboring points

            int nxd = (int)(av.transverseDiffusion[j] * holdnsize);

            av.n[pos - 1][i - j] += nxd;

            av.n[pos - 1][i + j] += nxd;

            m_sensor->AddSignal(-(nxd + Nholder) / 2, av.time,
                                av.time + av.zStepSize / av.velocity,
                                av.xgrid[i], 0, av.zgrid[pos], av.xgrid[i - j],
                                0, av.zgrid[pos - 1], false, true);

            m_sensor->AddSignal(-(nxd + Nholder) / 2, av.time,
                                av.time + av.zStepSize / av.velocity,
                                av.xgrid[i], 0, av.zgrid[pos], av.xgrid[i + j],
                                0, av.zgrid[pos - 1], false, true);

            chargeRemaining -= 2 * nxd;

          } else {  // For the initial x-position.

            av.n[pos - 1][i] += chargeRemaining;

            m_sensor->AddSignal(-(chargeRemaining + Nholder) / 2, av.time,
                                av.time + av.zStepSize / av.velocity,
                                av.xgrid[i], 0, av.zgrid[pos], av.xgrid[i], 0,
                                av.zgrid[pos - 1], false, true);
          }
        }

        av.N += holdnsize - Nholder;

      } else {
        // If the total avalanche size is smaller than the set saturation limit
        // the GetAvalancheSize function is utilized to obtain the size after
        // its propagation to the next z-coordinate grid point. Else, the size
        // will be kept constant under the propagation.
        if (av.N < m_MaxSize) {
          av.n[pos - 1][i] = GetAvalancheSize(av.zStepSize, av.n[pos][i],
                                              m_Townsend, m_Attachment);
        } else {
          av.n[pos - 1][i] = av.n[pos][i];
          m_Saturated = true;
          if (m_SaturationTime == -1)
            m_SaturationTime = av.time + abs(av.zStepSize / av.velocity);
        }
        // Produce induced signal on readout electrodes.
        m_sensor->AddSignal(-(av.n[pos][i] + Nholder) / 2, av.time,
                            av.time + av.zStepSize / av.velocity, av.xgrid[i],
                            0, av.zgrid[pos], av.xgrid[i], 0, av.zgrid[pos - 1],
                            false, true);

        av.N +=
            av.n[pos - 1][i] - Nholder;  // Update total number of electrons.
      }

      av.n[pos][i] = 0;  // Clear previous z-coordinate grid point.
    }
    pos -= 1;  // Update position index.
  }
  // After all active grid points have propagated, update the time.
  av.time += abs(av.zStepSize / av.velocity);
}

void AvalancheGrid::StartGridAvalanche() {
  // Start the AvalancheGrid algorithm.
  if ((!m_avmc && !m_driftAvalanche) || !m_sensor) return;

  if (!m_importAvalanche && m_avmc) ImportElectronData();

  std::cerr << m_className
            << "::StartGridAvalanche::Starting grid based simulation with "
            << m_avgrid.N << " initial electrons.\n";

  // The vector containing the indexes of the z-coordinate grid points of the
  // initial electrons needs to be ordered from small to large values. All
  // duplicate values need to be removed.
  sort(m_avgrid.gridPosition.begin(), m_avgrid.gridPosition.end());
  m_avgrid.gridPosition.erase(
      unique(m_avgrid.gridPosition.begin(), m_avgrid.gridPosition.end()),
      m_avgrid.gridPosition.end());
  // Set velocity if given.
  if (m_Velocity != 0) m_avgrid.velocity = m_Velocity;
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

void AvalancheGrid::SetGrid(const double zmin, const double zmax,
                            const int zsteps, const double xmin,
                            const double xmax, const int xsteps) {
  if (zmin >= zmax || zsteps <= 0) return;
  // Set grid
  if (!m_diffusion) {
    SetZGrid(m_avgrid, zmax, zmin, zsteps);

  } else {
    SetZGrid(m_avgrid, zmax, zmin, zsteps);
    SetXGrid(m_avgrid, xmax, xmin, xsteps);
  }
}

void AvalancheGrid::DiffusionFactors(Grid& av) {
  // Get transverse diffusion factors, yielding in the spreading of charge.
  if (!av.gridset || av.xStepSize <= 0) return;

  auto cdfunctop = TF1("cdftop", "ROOT::Math::normal_cdf(x, [0],[1])", -5, 5);

  cdfunctop.SetParameters(m_DiffSigma, 0.0);

  double factor = 1;
  int index = 0;

  while (factor > 1e-3) {  // 1e-3 is the precision cutoff.

    factor = cdfunctop.Eval(0 + av.xStepSize / 2 + index * av.xStepSize) -
             cdfunctop.Eval(0 - av.xStepSize / 2 + index * av.xStepSize);
    std::cerr << m_className
              << "::DiffusionFactors::Transvers diffusion factor: " << factor
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

void AvalancheGrid::DriftAvalanche(const double x, const double z,
                                   const double vz, const double t) {
  std::cerr << m_className << "::DriftElectron::Start avalanche at (t,x,z) =  ("
            << t << "," << x << "," << z << "), and speed v = " << vz << ".\n";

  m_driftAvalanche = true;

  if (m_avgrid.time == 0)
    std::cerr
        << m_className
        << "::DriftElectron::Overzriting start time of avalanche for t = 0 to "
        << t << ".\n";
  m_avgrid.time = t;

  SnapToGrid(m_avgrid, x, z, vz);
}
void AvalancheGrid::ImportElectronData() {
  // Get the information of the electrons from the AvalancheMicroscopic class.
  if (!m_avmc) return;

  m_importAvalanche = true;

  int np = m_avmc->GetNumberOfElectronEndpoints();

  std::cerr << m_className
            << "::ImportElectronData::Number of initial electrons = " << np
            << ".\n";

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
    SnapToGrid(m_avgrid, x2, z2, vel);  // Snap electrons to grid
  }
}

}  // namespace Garfield
