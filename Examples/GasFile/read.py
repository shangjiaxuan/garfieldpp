import ROOT
import os, sys
import ctypes

path = os.getenv('GARFIELD_INSTALL')
if sys.platform == 'darwin':
  ROOT.gSystem.Load(path + '/lib/libmagboltz.dylib')
  ROOT.gSystem.Load(path + '/lib/libGarfield.dylib')
else:
  ROOT.gSystem.Load(path + '/lib/libmagboltz.so')
  ROOT.gSystem.Load(path + '/lib/libGarfield.so')

# Set up the gas.
gas = ROOT.Garfield.MediumMagboltz()
gas.LoadGasFile('ar_80_co2_20_2T.gas')
gas.LoadIonMobility(path + '/share/Garfield/Data/IonMobility_Ar+_Ar.txt')
gas.PrintGas()

view = ROOT.Garfield.ViewMedium()
view.SetMedium(gas)
view.SetMagneticField(2.)
  
cV = ROOT.TCanvas('cV', '', 600, 600)
view.SetCanvas(cV)
view.PlotElectronVelocity('e')

cD = ROOT.TCanvas('cD', '', 600, 600)
view.SetCanvas(cD)
view.PlotElectronDiffusion('e')

cT = ROOT.TCanvas('cT', '', 600, 600)
view.SetCanvas(cT)
view.PlotElectronTownsend('e')

cA = ROOT.TCanvas('cA', '', 600, 600)
view.SetCanvas(cA)
view.PlotElectronAttachment('e')

