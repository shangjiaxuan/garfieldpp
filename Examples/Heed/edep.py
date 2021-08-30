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

ROOT.Garfield.SetDefaultStyle()

# Histograms
ROOT.TH1.StatOverflows(True)
hElectrons = ROOT.TH1F("hElectrons", "Number of electrons", 200, 0, 200)
hEdep = ROOT.TH1F("hEdep", "Energy Loss", 100, 0., 10.)

gas = ROOT.Garfield.MediumMagboltz()
gas.SetComposition("ar", 90., "co2", 10.)
gas.SetTemperature(293.15)
gas.SetPressure(760.)

# Width of the gas gap [cm].
width = 1.

# Make a component.
cmp = ROOT.Garfield.ComponentConstant()
cmp.SetArea(0., -10., -10., width, 10., 10.)
cmp.SetMedium(gas)
cmp.SetElectricField(100., 0., 0.)

# Make a sensor.
sensor = ROOT.Garfield.Sensor()
sensor.AddComponent(cmp)

# Set up HEED.
track = ROOT.Garfield.TrackHeed()
track.SetSensor(sensor)
track.SetParticle("pi")
track.SetMomentum(120.e9)

nEvents = 10000
for i in range(nEvents):
  if i % 1000 == 0: print i, "/", nEvents
  # Initial position and direction 
  x0 = 0.
  y0 = 0.
  z0 = 0.
  t0 = 0.
  dx0 = 1.
  dy0 = 0.
  dz0 = 0.
  track.NewTrack(x0, y0, z0, t0, dx0, dy0, dz0)
  # Cluster coordinates
  xc = ctypes.c_double(0.)
  yc = ctypes.c_double(0.)
  zc = ctypes.c_double(0.)
  tc = ctypes.c_double(0.)
  # Energy loss in a collision
  ec = ctypes.c_double(0.)
  # Dummy variable (not used at present)
  extra = ctypes.c_double(0.)
  # Number of electrons produced in a collision
  nc = ctypes.c_int(0)
  # Total energy loss along the track
  esum = 0.
  # Total number of electrons produced along the track
  nsum = 0
  # Loop over the clusters.
  while track.GetCluster(xc, yc, zc, tc, nc, ec, extra):
    esum += ec.value
    nsum += nc.value
  hElectrons.Fill(nsum)
  hEdep.Fill(esum * 1.e-3);
 
c1 = ROOT.TCanvas()
hElectrons.GetXaxis().SetTitle("number of electrons") 
hElectrons.Draw()
c1.SaveAs("ne.pdf")

c2 = ROOT.TCanvas()
hEdep.GetXaxis().SetTitle("energy loss [keV]")
hEdep.Draw()
c2.SaveAs("edep.pdf")

