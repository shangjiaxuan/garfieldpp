#include <fstream>
#include <iostream>
#include <algorithm>

#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLatex.h>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Numerics.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewBase.hh"
#include "Garfield/TrackSrim.hh"

namespace {

double StepVavilov(const double rkappa) {
  double xlmin = -3.7;
  if (rkappa < 0.1) {
    xlmin = -2.7;
  } else if (rkappa < 1) {
    xlmin = -2.9;
  } else if (rkappa < 2) {
    xlmin = -3.0;
  } else if (rkappa < 3) {
    xlmin = -3.1;
  } else if (rkappa < 4) {
    xlmin = -3.2;
  } else if (rkappa < 5) {
    xlmin = -3.3;
  } else if (rkappa < 6) {
    xlmin = -3.4;
  } else if (rkappa < 7) {
    xlmin = -3.5;
  } else if (rkappa < 8) {
    xlmin = -3.6;
  }
  return xlmin;
}

double Interpolate(const double x, const std::vector<double>& xtab,
                   const std::vector<double>& ytab) {
  if (x < xtab[0]) {
    return ytab[0];
  } else if (x > xtab.back()) {
    return ytab.back();
  }
  return Garfield::Numerics::Divdif(ytab, xtab, xtab.size(), x, 2);
}

void PrintSettings(const std::string& hdr, const double de, const double step,
                   const double ekin, const double beta2, const double gamma,
                   const double agas, const double zgas, const double density,
                   const double qp, const double mp, const double emax,
                   const double xi, const double kappa) {
  std::cout << hdr << "Settings:\n"
            << "    dE = " << de << " MeV,\n"
            << "    step = " << step << " cm.\n"
            << "    Ekin = " << ekin << " MeV,\n"
            << "    beta2 = " << beta2 << ",\n"
            << "    gamma = " << gamma << ".\n"
            << "    Agas = " << agas << ", Zgas = " << zgas << ",\n"
            << "    density = " << density << " g/cm3.\n"
            << "    Qpart = " << qp << ", mpart = " << 1.e-6 * mp << " MeV.\n"
            << "    Emax = " << emax << " MeV,\n"
            << "    xi = " << xi << " MeV,\n"
            << "    kappa = " << kappa << ".\n";
}
}

namespace Garfield {

TrackSrim::TrackSrim() : Track() { m_className = "TrackSrim"; }

bool TrackSrim::ReadFile(const std::string& file) {
  // SRMREA

  const std::string hdr = m_className + "::ReadFile:\n    ";
  // Open the material list.
  std::ifstream fsrim;
  fsrim.open(file.c_str(), std::ios::in);
  if (fsrim.fail()) {
    std::cerr << hdr << "Could not open SRIM file " << file
              << " for reading.\n    The file perhaps does not exist.\n";
    return false;
  }
  unsigned int nread = 0;

  // Read the header
  if (m_debug) {
    std::cout << hdr << "SRIM header records from file " << file << "\n";
  }
  constexpr size_t size = 100;
  char line[size];
  while (fsrim.getline(line, size, '\n')) {
    nread++;
    if (strstr(line, "SRIM version") != NULL) {
      if (m_debug) std::cout << "\t" << line << "\n";
    } else if (strstr(line, "Calc. date") != NULL) {
      if (m_debug) std::cout << "\t" << line << "\n";
    } else if (strstr(line, "Ion =") != NULL) {
      break;
    }
  }

  // Identify the ion
  char* token = NULL;
  token = strtok(line, " []=");
  token = strtok(NULL, " []=");
  token = strtok(NULL, " []=");
  // Set the ion charge.
  m_qion = std::atof(token);
  m_chargeset = true;
  token = strtok(NULL, " []=");
  token = strtok(NULL, " []=");
  token = strtok(NULL, " []=");
  // Set the ion mass (convert amu to eV).
  m_mion = std::atof(token) * AtomicMassUnitElectronVolt;

  // Find the target density
  if (!fsrim.getline(line, size, '\n')) {
    std::cerr << hdr << "Premature EOF looking for target density (line "
              << nread << ").\n";
    return false;
  }
  nread++;
  if (!fsrim.getline(line, size, '\n')) {
    std::cerr << hdr << "Premature EOF looking for target density (line "
              << nread << ").\n";
    return false;
  }
  nread++;
  const bool pre2013 = (strstr(line, "Target Density") != NULL); 
  token = strtok(line, " ");
  token = strtok(NULL, " ");
  token = strtok(NULL, " ");
  if (pre2013) token = strtok(NULL, " ");
  SetDensity(std::atof(token));

  // Check the stopping units
  while (fsrim.getline(line, size, '\n')) {
    nread++;
    if (strstr(line, "Stopping Units") == NULL) continue;
    if (strstr(line, "Stopping Units =  MeV / (mg/cm2)") != NULL ||
        strstr(line, "Stopping Units =  MeV/(mg/cm2)") != NULL) {
      if (m_debug) {
        std::cout << hdr << "Stopping units: MeV / (mg/cm2) as expected.\n";
      }
      break;
    }
    std::cerr << hdr << "Unknown stopping units. Aborting (line " << nread
              << ").\n";
    return false;
  }

  // Skip to the table
  while (fsrim.getline(line, size, '\n')) {
    nread++;
    if (strstr(line, "-----------") != NULL) break;
  }

  // Read the table line by line
  m_ekin.clear();
  m_emloss.clear();
  m_hdloss.clear();
  m_range.clear();
  m_transstraggle.clear();
  m_longstraggle.clear();
  unsigned int ntable = 0;
  while (fsrim.getline(line, size, '\n')) {
    nread++;
    if (strstr(line, "-----------") != NULL) break;
    // Energy
    token = strtok(line, " ");
    m_ekin.push_back(atof(token));
    token = strtok(NULL, " ");
    if (strcmp(token, "eV") == 0) {
      m_ekin[ntable] *= 1.0e-6;
    } else if (strcmp(token, "keV") == 0) {
      m_ekin[ntable] *= 1.0e-3;
    } else if (strcmp(token, "GeV") == 0) {
      m_ekin[ntable] *= 1.0e3;
    } else if (strcmp(token, "MeV") != 0) {
      std::cerr << hdr << "Unknown energy unit " << token << "; aborting\n";
      return false;
    }
    // EM loss
    token = strtok(NULL, " ");
    m_emloss.push_back(atof(token));
    // HD loss
    token = strtok(NULL, " ");
    m_hdloss.push_back(atof(token));
    // Projected range
    token = strtok(NULL, " ");
    m_range.push_back(atof(token));
    token = strtok(NULL, " ");
    if (strcmp(token, "A") == 0) {
      m_range[ntable] *= 1.0e-8;
    } else if (strcmp(token, "um") == 0) {
      m_range[ntable] *= 1.0e-4;
    } else if (strcmp(token, "mm") == 0) {
      m_range[ntable] *= 1.0e-1;
    } else if (strcmp(token, "m") == 0) {
      m_range[ntable] *= 1.0e2;
    } else if (strcmp(token, "km") == 0) {
      m_range[ntable] *= 1.0e5;
    } else if (strcmp(token, "cm") != 0) {
      std::cerr << hdr << "Unknown distance unit " << token << "; aborting\n";
      return false;
    }
    // Longitudinal straggling
    token = strtok(NULL, " ");
    m_longstraggle.push_back(atof(token));
    token = strtok(NULL, " ");
    if (strcmp(token, "A") == 0) {
      m_longstraggle[ntable] *= 1.0e-8;
    } else if (strcmp(token, "um") == 0) {
      m_longstraggle[ntable] *= 1.0e-4;
    } else if (strcmp(token, "mm") == 0) {
      m_longstraggle[ntable] *= 1.0e-1;
    } else if (strcmp(token, "m") == 0) {
      m_longstraggle[ntable] *= 1.0e2;
    } else if (strcmp(token, "km") == 0) {
      m_longstraggle[ntable] *= 1.0e5;
    } else if (strcmp(token, "cm") != 0) {
      std::cerr << hdr << "Unknown distance unit " << token << "; aborting\n";
      return false;
    }
    // Transverse straggling
    token = strtok(NULL, " ");
    m_transstraggle.push_back(atof(token));
    token = strtok(NULL, " ");
    if (strcmp(token, "A") == 0) {
      m_transstraggle[ntable] *= 1.0e-8;
    } else if (strcmp(token, "um") == 0) {
      m_transstraggle[ntable] *= 1.0e-4;
    } else if (strcmp(token, "mm") == 0) {
      m_transstraggle[ntable] *= 1.0e-1;
    } else if (strcmp(token, "m") == 0) {
      m_transstraggle[ntable] *= 1.0e2;
    } else if (strcmp(token, "km") == 0) {
      m_transstraggle[ntable] *= 1.0e5;
    } else if (strcmp(token, "cm") != 0) {
      std::cerr << hdr << "Unknown distance unit " << token << "; aborting\n";
      return false;
    }

    // Increment table line counter
    ++ntable;
  }

  // Find the scaling factor and convert to MeV/cm
  double scale = -1.;
  while (fsrim.getline(line, size, '\n')) {
    nread++;
    if (strstr(line, "=============") != NULL) {
      break;
    } else if (strstr(line, "MeV / (mg/cm2)") != NULL ||
               strstr(line, "MeV/(mg/cm2)") != NULL) {
      token = strtok(line, " ");
      scale = std::atof(token);
    }
  }
  if (scale < 0) {
    std::cerr << hdr << "Did not find stopping unit scaling; aborting.\n";
    return false;
  }
  scale *= 1.e3;
  for (unsigned int i = 0; i < ntable; ++i) {
    m_emloss[i] *= scale;
    m_hdloss[i] *= scale;
  }

  // Seems to have worked
  if (m_debug) {
    std::cout << hdr << "Successfully read " << file << "(" << nread
              << " lines).\n";
  }
  return true;
}

void TrackSrim::Print() {
  std::cout << "TrackSrim::Print:\n    SRIM energy loss table\n\n"
            << "    Energy     EM Loss     HD loss       Range  "
            << "l straggle  t straggle\n"
            << "     [MeV]    [MeV/cm]    [MeV/cm]        [cm] "
            << "      [cm]        [cm]\n\n";
  const unsigned int nPoints = m_emloss.size();
  for (unsigned int i = 0; i < nPoints; ++i) {
    printf("%10g  %10g  %10g  %10g  %10g  %10g\n", m_ekin[i],
           m_emloss[i] * m_density, m_hdloss[i] * m_density, m_range[i],
           m_longstraggle[i], m_transstraggle[i]);
  }
  std::cout << "\n";
  printf("    Work function:  %g eV\n", m_work);
  printf("    Fano factor:    %g\n", m_fano);
  printf("    Ion charge:     %g\n", m_qion);
  printf("    Mass:           %g MeV\n", 1.e-6 * m_mion);
  printf("    Density:        %g g/cm3\n", m_density);
  printf("    A, Z:           %g, %g\n", m_a, m_z);
}

void TrackSrim::PlotEnergyLoss() {

  const unsigned int nPoints = m_ekin.size();
  std::vector<double> yE;
  std::vector<double> yH;
  std::vector<double> yT;
  for (unsigned int i = 0; i < nPoints; ++i) {
    const double em = m_emloss[i] * m_density;
    const double hd = m_hdloss[i] * m_density;
    yE.push_back(em);
    yH.push_back(hd);
    yT.push_back(em + hd);
  }
  const double xmin = *std::min_element(std::begin(m_ekin), std::end(m_ekin));
  const double xmax = *std::max_element(std::begin(m_ekin), std::end(m_ekin));
  const double ymax = *std::max_element(std::begin(yT), std::end(yT));  
  // Prepare a plot frame.
  const std::string name = ViewBase::FindUnusedCanvasName("cSRIM"); 
  TCanvas* celoss = new TCanvas(name.c_str(), "Energy loss");
  celoss->SetLogx();
  celoss->SetGridx();
  celoss->SetGridy();
  celoss->DrawFrame(xmin, 0., xmax, 1.05 * ymax, ";Ion energy [MeV];Energy loss [MeV/cm]");

  // Make a graph for the 3 curves to plot.
  TGraph gr;
  gr.SetLineStyle(kSolid);
  gr.SetLineWidth(2);
  gr.SetMarkerStyle(21);
  gr.SetLineColor(kBlue + 1);
  gr.SetMarkerColor(kBlue + 1);
  gr.DrawGraph(nPoints, m_ekin.data(), yE.data(), "plsame");

  gr.SetLineColor(kGreen + 2);
  gr.SetMarkerColor(kGreen + 2);
  gr.DrawGraph(nPoints, m_ekin.data(), yH.data(), "plsame");

  gr.SetLineColor(kOrange - 3);
  gr.SetMarkerColor(kOrange - 3);
  gr.DrawGraph(nPoints, m_ekin.data(), yT.data(), "plsame");

  TLatex label;
  double xLabel = 0.4 * xmax;
  double yLabel = 0.9 * ymax;
  label.SetTextColor(kBlue + 1);
  label.SetText(xLabel, yLabel, "EM energy loss");
  label.DrawLatex(xLabel, yLabel, "EM energy loss");
  yLabel -= 1.5 * label.GetYsize();
  label.SetTextColor(kGreen + 2);
  label.DrawLatex(xLabel, yLabel, "HD energy loss");
  yLabel -= 1.5 * label.GetYsize();
  label.SetTextColor(kOrange - 3);
  label.DrawLatex(xLabel, yLabel, "Total energy loss");
  celoss->Update();
}

void TrackSrim::PlotRange() {

  const double xmin = *std::min_element(std::begin(m_ekin), std::end(m_ekin));
  const double xmax = *std::max_element(std::begin(m_ekin), std::end(m_ekin));
  const double ymax = *std::max_element(std::begin(m_range), std::end(m_range));

  // Prepare a plot frame.
  const std::string name = ViewBase::FindUnusedCanvasName("cSRIM");
  TCanvas* crange = new TCanvas(name.c_str(), "Range");
  crange->SetLogx();
  crange->SetGridx();
  crange->SetGridy();
  crange->DrawFrame(xmin, 0., xmax, 1.05 * ymax, ";Ion energy [MeV];Projected range [cm]");
  // Make a graph.
  TGraph gr;
  gr.SetLineColor(kOrange - 3);
  gr.SetMarkerColor(kOrange - 3);
  gr.SetLineStyle(kSolid);
  gr.SetLineWidth(2);
  gr.SetMarkerStyle(21);
  gr.DrawGraph(m_ekin.size(), m_ekin.data(), m_range.data(), "plsame");
  crange->Update();
}

void TrackSrim::PlotStraggling() {

  const double xmin = *std::min_element(std::begin(m_ekin), std::end(m_ekin));
  const double xmax = *std::max_element(std::begin(m_ekin), std::end(m_ekin));
  const double ymax = std::max(*std::max_element(std::begin(m_longstraggle),
                                                 std::end(m_longstraggle)),
                                *std::max_element(std::begin(m_transstraggle),
                                                  std::end(m_transstraggle)));
  // Prepare a plot frame.
  const std::string name = ViewBase::FindUnusedCanvasName("cSRIM");
  TCanvas* cstraggle = new TCanvas(name.c_str(), "Straggling");
  cstraggle->SetLogx();
  cstraggle->SetGridx();
  cstraggle->SetGridy();
  cstraggle->DrawFrame(xmin, 0., xmax, 1.05 * ymax, ";Ion energy [MeV];Straggling [cm]");

  // Make a graph for the 2 curves to plot.
  const unsigned int nPoints = m_ekin.size();
  TGraph gr;
  gr.SetLineStyle(kSolid);
  gr.SetLineWidth(2);
  gr.SetMarkerStyle(21);

  gr.SetLineColor(kOrange - 3);
  gr.SetMarkerColor(kOrange - 3);
  gr.DrawGraph(nPoints, m_ekin.data(), m_longstraggle.data(), "plsame");

  gr.SetLineColor(kGreen + 2);
  gr.SetMarkerColor(kGreen + 2);
  gr.DrawGraph(nPoints, m_ekin.data(), m_transstraggle.data(), "plsame");

  TLatex label;
  double xLabel = 1.2 * xmin;
  double yLabel = 0.9 * ymax;
  label.SetTextColor(kOrange - 3);
  label.SetText(xLabel, yLabel, "Longitudinal");
  label.DrawLatex(xLabel, yLabel, "Longitudinal");
  yLabel -= 1.5 * label.GetYsize();
  label.SetTextColor(kGreen + 2);
  label.DrawLatex(xLabel, yLabel, "Transverse");
  cstraggle->Update();
}

double TrackSrim::DedxEM(const double e) const {
  return Interpolate(e, m_ekin, m_emloss);
}

double TrackSrim::DedxHD(const double e) const {
  return Interpolate(e, m_ekin, m_hdloss);
}

double TrackSrim::Xi(const double x, const double beta2) const {

  constexpr double fconst = 1.e-6 * TwoPi * (
    FineStructureConstant * FineStructureConstant * HbarC * HbarC) / 
    (ElectronMass * AtomicMassUnit);
  return fconst * m_qion * m_qion * m_z * m_density * x / (m_a * beta2);
}

bool TrackSrim::PreciseLoss(const double step, const double estart,
                            double& deem, double& dehd) const {
  // SRMRKS

  const std::string hdr = m_className + "::PreciseLoss: ";
  // Debugging
  if (m_debug) {
    std::cout << hdr << "\n"
              << "    Initial energy: " << estart << " MeV\n"
              << "    Step: " << step << " cm\n";
  }
  // Precision aimed for.
  const double eps = 1.0e-2;
  // Number of intervals.
  unsigned int ndiv = 1;
  // Loop until precision achieved
  const unsigned int nMaxIter = 10;
  bool converged = false;
  for (unsigned int iter = 0; iter < nMaxIter; ++iter) {
    double e4 = estart;
    double e2 = estart;
    deem = 0.;
    dehd = 0.;
    // Compute rk2 and rk4 over the number of sub-divisions
    const double s = m_density * step / ndiv;
    for (unsigned int i = 0; i < ndiv; i++) {
      // rk2: initial point
      const double de21 = s * (DedxEM(e2) + DedxHD(e2));
      // Mid-way point
      const double em22 = s * DedxEM(e2 - 0.5 * de21);
      const double hd22 = s * DedxHD(e2 - 0.5 * de21);
      // Trace the rk2 energy
      e2 -= em22 + hd22;
      // rk4: initial point
      const double em41 = s * DedxEM(e4);
      const double hd41 = s * DedxHD(e4);
      const double de41 = em41 + hd41;
      // Mid-way point
      const double em42 = s * DedxEM(e4 - 0.5 * de41);
      const double hd42 = s * DedxHD(e4 - 0.5 * de41);
      const double de42 = em42 + hd42;
      // Second mid-point estimate
      const double em43 = s * DedxEM(e4 - 0.5 * de42);
      const double hd43 = s * DedxHD(e4 - 0.5 * de42);
      const double de43 = em43 + hd43;
      // End point estimate
      const double em44 = s * DedxEM(e4 - de43);
      const double hd44 = s * DedxHD(e4 - de43);
      const double de44 = em44 + hd44;
      // Store the energy loss terms (according to rk4)
      deem += (em41 + em44) / 6. + (em42 + em43) / 3.;
      dehd += (hd41 + hd44) / 6. + (hd42 + hd43) / 3.;
      // Store the new energy computed with rk4
      e4 -= (de41 + de44) / 6. + (de42 + de43) / 3.;
    }
    if (m_debug) {
      std::cout << hdr << "\n    Iteration " << iter << " has " << ndiv
                << " division(s). Losses:\n";
      printf("\tde4 = %12g, de2 = %12g MeV\n", estart - e2, estart - e4);
      printf("\tem4 = %12g, hd4 = %12g MeV\n", deem, dehd);
    }
    // Compare the two estimates
    if (fabs(e2 - e4) > eps * (fabs(e2) + fabs(e4) + fabs(estart))) {
      // Repeat with twice the number of steps.
      ndiv *= 2;
    } else {
      converged = true;
      break;
    }
  }

  if (!converged) {
    std::cerr << hdr << "No convergence achieved integrating energy loss.\n";
  } else if (m_debug) {
    std::cout << hdr << "Convergence at eps = " << eps << "\n";
  }
  return converged;
}

bool TrackSrim::EstimateRange(const double ekin, const double step,
                              double& stpmax) {
  // Find distance over which the ion just does not lose all its energy
  // ekin       : Kinetic energy [MeV]
  // step       : Step length as guessed [cm]
  // stpmax     : Maximum step
  // SRMDEZ

  const std::string hdr = m_className + "::EstimateRange: ";
  // Initial estimate
  stpmax = step;

  // Find the energy loss expected for the present step length.
  double st1 = step;
  double deem = 0., dehd = 0.;
  PreciseLoss(st1, ekin, deem, dehd);
  double de1 = deem + dehd;
  // Do nothing if this is ok
  if (de1 < ekin) {
    if (m_debug) std::cout << hdr << "Initial step OK.\n";
    return true;
  }
  // Find a smaller step for which the energy loss is less than EKIN.
  double st2 = 0.5 * step;
  double de2 = de1;
  const unsigned int nMaxIter = 20;
  for (unsigned int iter = 0; iter < nMaxIter; ++iter) {
    // See where we stand
    PreciseLoss(st2, ekin, deem, dehd);
    de2 = deem + dehd;
    // Below the kinetic energy: done
    if (de2 < ekin) break;
    // Not yet below the kinetic energy: new iteration.
    st1 = st2;
    de1 = de2;
    st2 *= 0.5;
  }
  if (de2 >= ekin) {
    std::cerr << hdr << "\n    Did not find a smaller step in " << nMaxIter
              << " iterations. Abandoned.\n";
    stpmax = 0.5 * (st1 + st2);
    return false;
  }
  if (m_debug)
    printf("\tstep 1 = %g cm, de 1 = %g MeV\n\tstep 2 = %g cm, de 2 = %g MeV\n",
           st1, de1 - ekin, st2, de2 - ekin);

  // Now perform a bisection
  for (unsigned int iter = 0; iter < nMaxIter; ++iter) {
    // Avoid division by zero.
    if (de2 == de1) {
      if (m_debug) {
        std::cerr << hdr << "Bisection failed due to equal energy loss for "
                  << "two step sizes. Abandoned.\n";
      }
      stpmax = 0.5 * (st1 + st2);
      return false;
    }
    // Estimate step to give total energy loss.
    double st3;
    if ((fabs(de1 - ekin) < 0.01 * fabs(de2 - de1)) ||
        (fabs(de1 - ekin) > 0.99 * fabs(de2 - de1))) {
      st3 = 0.5 * (st1 + st2);
    } else {
      st3 = st1 - (st2 - st1) * (de1 - ekin) / (de2 - de1);
    }
    // See how well we are doing.
    PreciseLoss(st3, ekin, deem, dehd);
    const double de3 = deem + dehd;
    if (m_debug) {
      std::printf("\tStep 1 = %g cm, dE 1 = %g MeV\n", st1, de1 - ekin);
      std::printf("\tStep 2 = %g cm, dE 2 = %g MeV\n", st2, de2 - ekin);
      std::printf("\tStep 3 = %g cm, dE 3 = %g MeV\n", st3, de3 - ekin);
    }
    //  Update the estimates above and below.
    if (de3 > ekin) {
      st1 = st3;
      de1 = de3;
    } else {
      st2 = st3;
      de2 = de3;
    }
    // See whether we've converged.
    if (fabs(de3 - ekin) < 1e-3 * (fabs(de3) + fabs(ekin)) ||
        fabs(st1 - st2) < 1e-3 * (fabs(st1) + fabs(st2))) {
      stpmax = st1 - (st2 - st1) * (de1 - ekin) / (de2 - de1);
      return true;
    }
  }
  if (m_debug) {
    std::cout << hdr << "Bisection did not converge in " << nMaxIter
              << " steps. Abandoned.\n";
  }
  stpmax = st1 - (st2 - st1) * (de1 - ekin) / (de2 - de1);
  return false;
}

bool TrackSrim::NewTrack(const double x0, const double y0, const double z0,
                         const double t0, const double dx0, const double dy0,
                         const double dz0) {
  // Generates electrons for a SRIM track
  // SRMGEN
  const std::string hdr = m_className + "::NewTrack: ";

  // Verify that a sensor has been set.
  if (!m_sensor) {
    std::cerr << hdr << "Sensor is not defined.\n";
    return false;
  }

  // Get the bounding box.
  double xmin = 0., ymin = 0., zmin = 0.;
  double xmax = 0., ymax = 0., zmax = 0.;
  if (!m_sensor->GetArea(xmin, ymin, zmin, xmax, ymax, zmax)) {
    std::cerr << hdr << "Drift area is not set.\n";
    return false;
  } else if (x0 < xmin || x0 > xmax || y0 < ymin || y0 > ymax || z0 < zmin ||
             z0 > zmax) {
    std::cerr << hdr << "Initial position outside bounding box.\n";
    return false;
  }

  // Make sure the initial position is inside an ionisable medium.
  Medium* medium = nullptr;
  if (!m_sensor->GetMedium(x0, y0, z0, medium)) {
    std::cerr << hdr << "No medium at initial position.\n";
    return false;
  } else if (!medium->IsIonisable()) {
    std::cerr << hdr << "Medium at initial position is not ionisable.\n";
    return false;
  }

  // Normalise and store the direction.
  const double normdir = sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0);
  double xdir = dx0;
  double ydir = dy0;
  double zdir = dz0;
  if (normdir < Small) {
    if (m_debug) {
      std::cout << hdr << "Direction vector has zero norm.\n"
                << "    Initial direction is randomized.\n";
    }
    // Null vector. Sample the direction isotropically.
    RndmDirection(xdir, ydir, zdir);
  } else {
    // Normalise the direction vector.
    xdir /= normdir;
    ydir /= normdir;
    zdir /= normdir;
  }

  // Make sure all necessary parameters have been set.
  if (m_mion < Small) {
    std::cerr << hdr << "Particle mass not set.\n";
    return false;
  } else if (!m_chargeset) {
    std::cerr << hdr << "Particle charge not set.\n";
    return false;
  } else if (m_energy < Small) {
    std::cerr << hdr << "Initial particle energy not set.\n";
    return false;
  } else if (m_work < Small) {
    std::cerr << hdr << "Work function not set.\n";
    return false;
  } else if (m_a < Small || m_z < Small) {
    std::cerr << hdr << "A and/or Z not set.\n";
    return false;
  }
  // Check the initial energy (in MeV).
  const double ekin0 = 1.e-6 * GetKineticEnergy();
  if (ekin0 < 1.e-14 * m_mion || ekin0 < 1.e-6 * m_work) {
    std::cerr << hdr << "Initial kinetic energy E = " << ekin0
              << " MeV such that beta2 = 0 or E << W; particle stopped.\n";
    return true;
  }

  // Get an upper limit for the track length.
  const double tracklength = 10 * Interpolate(ekin0, m_ekin, m_range);

  // Header of debugging output.
  if (m_debug) {
    std::cout << hdr << "Track generation with the following parameters:\n";
    const unsigned int nTable = m_ekin.size();
    printf("      Table size           %u\n", nTable);
    printf("      Particle kin. energy %g MeV\n", ekin0);
    printf("      Particle mass        %g MeV\n", 1.e-6 * m_mion);
    printf("      Particle charge      %g\n", m_qion);
    printf("      Work function        %g eV\n", m_work);
    if (m_fano > 0.) {
      printf("      Fano factor          %g\n", m_fano);
    } else {
      std::cout << "      Fano factor          Not set\n";
    }
    printf("      Long. straggling:    %d\n", m_useLongStraggle);
    printf("      Trans. straggling:   %d\n", m_useTransStraggle);
    printf("      Cluster size         %d\n", m_nsize);
  }

  // Plot.
  if (m_viewer) PlotNewTrack(x0, y0, z0);
 
  // Reset the cluster count.
  m_currcluster = 0;
  m_clusters.clear();

  // Initial situation: starting position
  double x = x0;
  double y = y0;
  double z = z0;
  double t = t0;
  // Store the energy [MeV].
  double e = ekin0;
  // Total distance covered
  double dsum = 0.0;
  // Pool of unused energy
  double epool = 0.0;

  // Loop generating clusters
  int iter = 0;
  while (iter < m_maxclusters || m_maxclusters < 0) {
    // Work out what the energy loss per cm, straggling and projected range are
    // at the start of the step.
    const double dedxem = DedxEM(e) * m_density;
    const double dedxhd = DedxHD(e) * m_density;
    const double prange = Interpolate(e, m_ekin, m_range);
    double strlon = Interpolate(e, m_ekin, m_longstraggle);
    double strlat = Interpolate(e, m_ekin, m_transstraggle);

    if (!m_useLongStraggle) strlon = 0;
    if (!m_useTransStraggle) strlat = 0;

    if (m_debug) {
      std::cout << hdr << "\n    Energy = " << e
                << " MeV,\n    dEdx em, hd = " << dedxem << ", " << dedxhd
                << " MeV/cm,\n    e-/cm = " << 1.e6 * dedxem / m_work
                << ".\n    Straggling long/lat: " << strlon << ", " << strlat
                << " cm\n";
    }
    // Find the step size for which we get approximately the target # clusters.
    double step;
    if (m_nsize > 0) {
      step = m_nsize * 1.e-6 * m_work / dedxem;
    } else {
      const double ncls = m_maxclusters > 0 ? 0.5 * m_maxclusters : 100;
      step = ekin0 / (ncls * (dedxem + dedxhd));
    }
    // Truncate if this step exceeds the length.
    bool finish = false;
    // Make an accurate integration of the energy loss over the step.
    double deem = 0., dehd = 0.;
    PreciseLoss(step, e, deem, dehd);
    // If the energy loss exceeds the particle energy, truncate step.
    double stpmax;
    if (deem + dehd > e) {
      EstimateRange(e, step, stpmax);
      step = stpmax;
      PreciseLoss(step, e, deem, dehd);
      deem = e * deem / (dehd + deem);
      dehd = e - deem;
      finish = true;
      if (m_debug) std::cout << hdr << "Finish raised. Track length reached.\n";
    } else {
      stpmax = tracklength - dsum;
    }
    if (m_debug) {
      std::cout << hdr << "Maximum step size set to " << stpmax << " cm.\n";
    }
    // Ensure that this is larger than the minimum modelable step size.
    double stpmin;
    if (!SmallestStep(e, deem, step, stpmin)) {
      std::cerr << hdr << "Failure computing the minimum step size."
                << "\n    Clustering abandoned.\n";
      return false;
    }

    double eloss;
    if (stpmin > stpmax) {
      // No way to find a suitable step size: use fixed energy loss.
      if (m_debug) std::cout << hdr << "stpmin > stpmax. Deposit all energy.\n";
      eloss = deem;
      if (e - eloss - dehd < 0) eloss = e - dehd;
      finish = true;
      if (m_debug) std::cout << hdr << "Finish raised. Single deposit.\n";
    } else if (step < stpmin) {
      // If needed enlarge the step size
      if (m_debug) std::cout << hdr << "Enlarging step size.\n";
      step = stpmin;
      PreciseLoss(step, e, deem, dehd);
      if (deem + dehd > e) {
        if (m_debug) std::cout << hdr << "Excess loss. Recomputing stpmax.\n";
        EstimateRange(e, step, stpmax);
        step = stpmax;
        PreciseLoss(step, e, deem, dehd);
        deem = e * deem / (dehd + deem);
        dehd = e - deem;
        eloss = deem;
      } else {
        eloss = RndmEnergyLoss(e, deem, step);
      }
    } else {
      // Draw an actual energy loss for such a step.
      if (m_debug) std::cout << hdr << "Using existing step size.\n";
      eloss = RndmEnergyLoss(e, deem, step);
    }
    // Ensure we are neither below 0 nor above the total energy.
    if (eloss < 0) {
      if (m_debug) std::cout << hdr << "Truncating negative energy loss.\n";
      eloss = 0;
    } else if (eloss > e - dehd) {
      if (m_debug) std::cout << hdr << "Excess energy loss, using mean.\n";
      eloss = deem;
      if (e - eloss - dehd < 0) {
        eloss = e - dehd;
        finish = true;
        if (m_debug) std::cout << hdr << "Finish raised. Using mean energy.\n";
      }
    }
    if (m_debug) {
      std::cout << hdr << "Step length = " << step << " cm.\n    "
                << "Mean loss =   " << deem << " MeV.\n    "
                << "Actual loss = " << eloss << " MeV.\n";
    }

    // Check that the cluster is in an ionisable medium and within bounding box
    if (!m_sensor->GetMedium(x, y, z, medium)) {
      if (m_debug) {
        std::cout << hdr << "No medium at position (" << x << "," << y << ","
                  << z << ").\n";
      }
      break;
    } else if (!medium->IsIonisable()) {
      if (m_debug) {
        std::cout << hdr << "Medium at (" << x << "," << y << "," << z
                  << ") is not ionisable.\n";
      }
      break;
    } else if (!m_sensor->IsInArea(x, y, z)) {
      if (m_debug) {
        std::cout << hdr << "Cluster at (" << x << "," << y << "," << z
                  << ") outside bounding box.\n";
      }
      break;
    }
    // Add a cluster.
    Cluster cluster;
    cluster.x = x;
    cluster.y = y;
    cluster.z = z;
    cluster.t = t;
    if (m_fano < Small) {
      // No fluctuations.
      cluster.electrons = int((eloss + epool) / (1.e-6 * m_work));
      cluster.ec = m_work * cluster.electrons;
    } else {
      double ecl = 1.e6 * (eloss + epool);
      cluster.electrons = 0.0;
      cluster.ec = 0.0;
      while (true) {
        // if (cluster.ec < 100) printf("ec = %g\n", cluster.ec);
        const double ernd1 = RndmHeedWF(m_work, m_fano);
        if (ernd1 > ecl) break;
        cluster.electrons++;
        cluster.ec += ernd1;
        ecl -= ernd1;
      }
      if (m_debug)
        std::cout << hdr << "EM + pool: " << 1.e6 * (eloss + epool)
                  << " eV, W: " << m_work
                  << " eV, E/w: " << (eloss + epool) / (1.0e-6 * m_work)
                  << ", n: " << cluster.electrons << ".\n";
    }
    cluster.kinetic = e;
    epool += eloss - 1.e-6 * cluster.ec;
    if (m_debug) {
      std::cout << hdr << "Cluster " << m_clusters.size() << "\n    at ("
                << cluster.x << ", " << cluster.y << ", " << cluster.z
                << "),\n    e = " << cluster.ec
                << ",\n    n = " << cluster.electrons
                << ",\n    pool = " << epool << " MeV.\n";
    }
    m_clusters.push_back(std::move(cluster));
    if (m_viewer) PlotCluster(x, y, z);

    // Keep track of the length and energy
    dsum += step;
    e -= eloss + dehd;
    if (finish) {
      // Stop if the flag is raised
      if (m_debug) std::cout << hdr << "Finishing flag raised.\n";
      break;
    } else if (e < ekin0 * 1.e-9) {
      // No energy left
      if (m_debug) std::cout << hdr << "Energy exhausted.\n";
      break;
    }
    // Draw scattering distances
    const double scale = sqrt(step / prange);
    const double sigt1 = RndmGaussian(0., scale * strlat);
    const double sigt2 = RndmGaussian(0., scale * strlat);
    const double sigl = RndmGaussian(0., scale * strlon);
    if (m_debug)
      std::cout << hdr << "sigma l, t1, t2: " << sigl << ", " << sigt1 << ", "
                << sigt2 << "\n";
    // Rotation angles to bring z-axis in line
    double theta, phi;
    if (xdir * xdir + zdir * zdir <= 0) {
      if (ydir < 0) {
        theta = -HalfPi;
      } else if (ydir > 0) {
        theta = +HalfPi;
      } else {
        std::cerr << hdr << "Zero step length; clustering abandoned.\n";
        return false;
      }
      phi = 0;
    } else {
      phi = atan2(xdir, zdir);
      theta = atan2(ydir, sqrt(xdir * xdir + zdir * zdir));
    }

    // Update the position.
    const double cp = cos(phi);
    const double ct = cos(theta);
    const double sp = sin(phi);
    const double st = sin(theta);
    x += step * xdir + cp * sigt1 - sp * st * sigt2 + sp * ct * sigl;
    y += step * ydir + ct * sigt2 + st * sigl;
    z += step * zdir - sp * sigt1 - cp * st * sigt2 + cp * ct * sigl;
    // Update the time.
    const double rk = 1.e6 * e / m_mion;
    const double gamma = 1. + rk;
    const double beta2 = rk > 1.e-5 ? 1. - 1. / (gamma * gamma) : 2. * rk;
    const double vmag = sqrt(beta2) * SpeedOfLight; 
    if (vmag > 0.) t += step / vmag;
    // (Do not) update direction
    if (false) {
      xdir = step * xdir + cp * sigt1 - sp * st * sigt2 + sp * ct * sigl;
      ydir = step * ydir + ct * sigt2 + st * sigl;
      zdir = step * zdir - sp * sigt1 - cp * st * sigt2 + cp * ct * sigl;
      double dnorm = sqrt(xdir * xdir + ydir * ydir + zdir * zdir);
      if (dnorm <= 0) {
        std::cerr << hdr << "Zero step length; clustering abandoned.\n";
        return false;
      }
      xdir = xdir / dnorm;
      ydir = ydir / dnorm;
      zdir = zdir / dnorm;
    }
    // Next cluster
    iter++;
  }
  if (iter == m_maxclusters) {
    std::cerr << hdr << "Exceeded maximum number of clusters.\n";
  }
  return true;
  // finished generating
}

bool TrackSrim::SmallestStep(const double ekin, double de, double step,
                             double& stpmin) {
  // Determines the smallest step size for which there is little
  // or no risk of finding negative energy fluctuations.
  // SRMMST

  const std::string hdr = m_className + "::SmallestStep: ";
  constexpr double expmax = 30;

  // By default, assume the step is right.
  stpmin = step;
  // Check correctness.
  if (ekin <= 0 || de <= 0 || step <= 0) {
    std::cerr << hdr << "Input parameters not valid.\n    Ekin = " << ekin
              << " MeV, dE = " << de << " MeV, step length = " << step
              << " cm.\n";
    return false;
  } else if (m_mion <= 0 || fabs(m_qion) <= 0) {
    std::cerr << hdr
              << "Track parameters not valid.\n    Mass = " << 1.e-6 * m_mion
              << " MeV, charge = " << m_qion << ".\n";
    return false;
  } else if (m_a <= 0 || m_z <= 0 || m_density <= 0) {
    std::cerr << hdr << "Gas parameters not valid.\n    A = " << m_a
              << ", Z = " << m_z << " density = " << m_density << " g/cm3.\n";
    return false;
  }

  // Basic kinematic parameters
  const double rkin = 1.e6 * ekin / m_mion;
  const double gamma = 1. + rkin;
  const double beta2 = rkin > 1.e-5 ? 1. - 1. / (gamma * gamma) : 2. * rkin;

  // Compute the maximum energy transfer [MeV]
  const double rm = ElectronMass / m_mion;
  const double emax = 2 * ElectronMass * 1.e-6 * beta2 * gamma * gamma /
                      (1. + 2 * gamma * rm + rm * rm);
  // Compute the Rutherford term
  double xi = Xi(step, beta2);
  // Compute the scaling parameter
  double rkappa = xi / emax;
  // Step size and energy loss
  double denow = de;
  double stpnow = step;
  constexpr unsigned int nMaxIter = 10;
  for (unsigned int iter = 0; iter < nMaxIter; ++iter) {
    bool retry = false;
    // Debugging output.
    if (m_debug) {
      PrintSettings(hdr, denow, stpnow, ekin, beta2, gamma, m_a, m_z, m_density,
                    m_qion, m_mion, emax, xi, rkappa);
    }
    double xinew = xi;
    double rknew = rkappa;
    if (m_model <= 0 || m_model > 4) {
      // No fluctuations: any step is permitted
      stpmin = stpnow;
    } else if (m_model == 1) {
      // Landau distribution
      constexpr double xlmin = -3.;
      const double exponent = -xlmin - 1. + Gamma - beta2 - denow / xi;
      const double rklim = exponent < -expmax ? 0. : exp(exponent);
      stpmin = stpnow * (rklim / rkappa);
      if (m_debug) {
        std::cout << hdr << "Landau distribution is imposed.\n    kappa_min = "
                  << rklim << ", d_min = " << stpmin << " cm.\n";
      }
    } else if (m_model == 2) {
      // Vavilov distribution, ensure we're in range.
      const double xlmin = StepVavilov(rkappa);
      const double exponent = -xlmin - 1. + Gamma - beta2 - denow / xi;
      const double rklim = exponent < -expmax ? 0. : exp(exponent);
      stpmin = stpnow * (rklim / rkappa);
      xinew = Xi(stpmin, beta2);
      rknew = xinew / emax;
      if (m_debug) {
        std::cout << hdr << "Vavilov distribution is imposed.\n    kappa_min = "
                  << rklim << ", d_min = " << stpmin
                  << " cm\n    kappa_new = " << rknew << ", xi_new = " << xinew
                  << " MeV.\n";
      }
      if (stpmin > stpnow * 1.1) {
        if (m_debug) std::cout << hdr << "Step size increase. New pass.\n";
        retry = true;
      }
    } else if (m_model == 3) {
      // Gaussian model
      const double sigma2 = xi * emax * (1 - 0.5 * beta2);
      stpmin = stpnow * 16 * sigma2 / (denow * denow);
      if (m_debug) {
        const double sigmaMin2 = Xi(stpmin, beta2) * emax * (1 - 0.5 * beta2);
        std::cout << hdr << "Gaussian distribution is imposed.\n    "
                  << "d_min = " << stpmin << " cm.\n    sigma/mu (old) = "
                  << sqrt(sigma2) / de << ",\n    sigma/mu (min) = " 
                  << sqrt(sigmaMin2) / (stpmin * denow / stpnow) << "\n";
      }
    } else if (rkappa < 0.05) {
      // Combined model: for low kappa, use the Landau distribution.
      constexpr double xlmin = -3.;
      const double exponent = -xlmin - 1. + Gamma - beta2 - denow / xi;
      const double rklim = exponent < -expmax ? 0. : exp(exponent);
      stpmin = stpnow * (rklim / rkappa);
      xinew = Xi(stpmin, beta2);
      rknew = xinew / emax;
      if (m_debug) {
        std::cout << hdr << "Landau distribution automatic.\n    kappa_min = " 
                  << rklim << ", d_min = " << stpmin << " cm.\n";
      }
      if (rknew > 0.05 || stpmin > stpnow * 1.1) {
        retry = true;
        if (m_debug) {
          std::cout << hdr << "Model change or step increase. New pass.\n";
        }
      }
    } else if (rkappa < 5) {
      // For medium kappa, use the Vavilov distribution
      const double xlmin = StepVavilov(rkappa);
      const double exponent = -xlmin - 1. + Gamma - beta2 - denow / xi;
      const double rklim = exponent < -expmax ? 0. : exp(exponent);
      stpmin = stpnow * (rklim / rkappa);
      xinew = Xi(stpmin, beta2);
      rknew = xinew / emax;
      if (m_debug) {
        std::cout << hdr << "Vavilov distribution automatic.\n    kappa_min = "
                  << rklim << ", d_min = " << stpmin << " cm\n    kappa_new = " 
                  << rknew << ", xi_new = " << xinew << " MeV.\n";
      }
      if (rknew > 5 || stpmin > stpnow * 1.1) {
        retry = true;
        if (m_debug) {
          std::cout << hdr << "Model change or step increase. New pass.\n";
        }
      }
    } else {
      // And for large kappa, use the Gaussian values.
      const double sigma2 = xi * emax * (1 - 0.5 * beta2);
      stpmin = stpnow * 16 * sigma2 / (denow * denow);
      if (m_debug) {
        const double sigmaMin2 = Xi(stpmin, beta2) * emax * (1 - 0.5 * beta2);
        std::cout << hdr << "Gaussian distribution automatic.\n    "
                  << "d_min = " << stpmin << " cm.\n    sigma/mu (old) = "
                  << sqrt(sigma2) / de << ",\n    sigma/mu (min) = " 
                  << sqrt(sigmaMin2) / (stpmin * denow / stpnow) << "\n";
      }
    }
    // See whether we should do another pass.
    if (stpnow > stpmin) {
      if (m_debug) {
        std::cout << hdr << "Step size ok, minimum: " << stpmin << " cm\n";
      }
      break;
    }
    if (!retry) {
      if (m_debug) {
        std::cerr << hdr << "Step size must be increased to " << stpmin
                  << "cm.\n";
      }
      break;
    }
    // New iteration
    rkappa = rknew;
    xi = xinew;
    denow *= stpmin / stpnow;
    stpnow = stpmin;
    if (m_debug) std::cout << hdr << "Iteration " << iter << "\n";
    if (iter == nMaxIter - 1) {
      // Need interation, but ran out of tries
      std::cerr << hdr << "No convergence reached on step size.\n";
    }
  }
  return true;
}

double TrackSrim::RndmEnergyLoss(const double ekin, const double de,
                                 const double step) const {
  //   RNDDE  - Generates a random energy loss.
  //   VARIABLES : EKIN       : Kinetic energy [MeV]
  //            DE         : Mean energy loss over the step [MeV]
  //            STEP       : Step length [cm]
  //            BETA2      : Velocity-squared
  //            GAMMA      : Projectile gamma
  //            EMAX       : Maximum energy transfer per collision [MeV]
  //            XI         : Rutherford term [MeV]
  //            FCONST     : Proportionality constant
  //            EMASS      : Electron mass [MeV]

  const std::string hdr = "TrackSrim::RndmEnergyLoss: ";
  // Check correctness.
  if (ekin <= 0 || de <= 0 || step <= 0) {
    std::cerr << hdr << "Input parameters not valid.\n    Ekin = " << ekin
              << " MeV, dE = " << de << " MeV, step length = " << step
              << " cm.\n";
    return 0.;
  } else if (m_mion <= 0 || fabs(m_qion) <= 0) {
    std::cerr << hdr << "Track parameters not valid.\n    Mass = " 
              << m_mion << " MeV, charge = " << m_qion << ".\n";
    return 0.;
  } else if (m_a <= 0 || m_z <= 0 || m_density <= 0) {
    std::cerr << hdr << "Material parameters not valid.\n    A = " << m_a
              << ", Z = " << m_z << ", density = " << m_density << " g/cm3.\n";
    return 0.;
  }
  // Basic kinematic parameters
  const double rkin = 1.e6 * ekin / m_mion;
  const double gamma = 1. + rkin;
  const double beta2 = rkin > 1.e-5 ? 1. - 1. / (gamma * gamma) : 2. * rkin;

  // Compute maximum energy transfer
  const double rm = ElectronMass / m_mion;
  const double emax = 2 * ElectronMass * 1.e-6 * beta2 * gamma * gamma /
                      (1. + 2 * gamma * rm + rm * rm);
  // Compute the Rutherford term
  const double xi = Xi(step, beta2);
  // Compute the scaling parameter
  const double rkappa = xi / emax;
  // Debugging output.
  if (m_debug) {
    PrintSettings(hdr, de, step, ekin, beta2, gamma, m_a, m_z, m_density, 
                  m_qion, m_mion, emax, xi, rkappa);
  }
  double rndde = de;
  if (m_model <= 0 || m_model > 4) {
    // No fluctuations.
    if (m_debug) std::cout << hdr << "Fixed energy loss.\n";
  } else if (m_model == 1) {
    // Landau distribution
    if (m_debug) std::cout << hdr << "Landau imposed.\n";
    const double xlmean = -(log(rkappa) + beta2 + 1. - Gamma);
    rndde += xi * (RndmLandau() - xlmean);
  } else if (m_model == 2) {
    // Vavilov distribution, ensure we are in range.
    if (m_debug) std::cout << hdr << "Vavilov imposed.\n";
    if (rkappa > 0.01 && rkappa < 12) {
      const double xvav = RndmVavilov(rkappa, beta2);
      rndde += xi * (xvav + log(rkappa) + beta2 + (1 - Gamma));
    }
  } else if (m_model == 3) {
    // Gaussian model
    if (m_debug) std::cout << hdr << "Gaussian imposed.\n";
    rndde += RndmGaussian(0., sqrt(xi * emax * (1 - 0.5 * beta2)));
  } else if (rkappa < 0.05) {
    // Combined model: for low kappa, use the landau distribution.
    if (m_debug) std::cout << hdr << "Landau automatic.\n";
    const double xlmean = -(log(rkappa) + beta2 + (1 - Gamma));
    const double par[] = {0.50884,    1.26116, 0.0346688,  1.46314,
                          0.15088e-2, 1.00324, -0.13049e-3};
    const double xlmax = par[0] + par[1] * xlmean + par[2] * xlmean * xlmean +
                         par[6] * xlmean * xlmean * xlmean +
                         (par[3] + xlmean * par[4]) * exp(par[5] * xlmean);
    double xlan = RndmLandau();
    for (unsigned int iter = 0; iter < 100; ++iter) {
      if (xlan < xlmax) break;
      xlan = RndmLandau();
    }
    rndde += xi * (xlan - xlmean);
  } else if (rkappa < 5) {
    // For medium kappa, use the Vavilov distribution.
    if (m_debug) std::cout << hdr << "Vavilov fast automatic.\n";
    const double xvav = RndmVavilov(rkappa, beta2);
    rndde += xi * (xvav + log(rkappa) + beta2 + (1 - Gamma));
  } else {
    // And for large kappa, use the Gaussian values.
    if (m_debug) std::cout << hdr << "Gaussian automatic.\n";
    rndde = RndmGaussian(de, sqrt(xi * emax * (1 - 0.5 * beta2)));
  }
  // Debugging output
  if (m_debug) {
    std::cout << hdr << "Energy loss generated = " << rndde << " MeV.\n";
  }
  return rndde;
}

bool TrackSrim::GetCluster(double& xcls, double& ycls, double& zcls,
                           double& tcls, int& n, double& e, double& extra) {
  if (m_debug) {
    std::cout << m_className << "::GetCluster: Cluster " << m_currcluster
              << " of " << m_clusters.size() << "\n";
  }
  // Stop if we have exhausted the list of clusters.
  if (m_currcluster >= m_clusters.size()) return false;

  const auto& cluster = m_clusters[m_currcluster];
  xcls = cluster.x;
  ycls = cluster.y;
  zcls = cluster.z;
  tcls = cluster.t;

  n = cluster.electrons;
  e = cluster.ec;
  extra = cluster.kinetic;
  // Move to next cluster
  ++m_currcluster;
  return true;
}
}
