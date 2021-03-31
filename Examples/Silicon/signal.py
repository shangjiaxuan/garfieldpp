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

si = ROOT.Garfield.MediumSilicon()
si.SetTemperature(293.)

mediumView = ROOT.Garfield.ViewMedium()
cM = ROOT.TCanvas('cM', '', 600, 600)
mediumView.SetCanvas(cM)
mediumView.SetMedium(si)
mediumView.PlotElectronVelocity('e')
mediumView.PlotHoleVelocity('e', True)

# Thickness of the silicon [cm]
d = 100.e-4
box = ROOT.Garfield.SolidBox(0, 0.5 * d, 0, 2 * d, 0.5 * d, 2 * d)
geo = ROOT.Garfield.GeometrySimple()
geo.AddSolid(box, si)

# Make a component with constant drift field and weighting field.
# Bias voltage [V]
vbias = -50.;
uniformField = ROOT.Garfield.ComponentConstant()
uniformField.SetGeometry(geo);
uniformField.SetElectricField(0, vbias / d, 0)
uniformField.SetWeightingField(0, -1. / d, 0, 'pad')

# Make a component with analytic weighting field for a strip or pixel.
pitch = 55.e-4
halfpitch = 0.5 * pitch
wField = ROOT.Garfield.ComponentAnalyticField()
wField.SetGeometry(geo)
wField.AddPlaneY(0, vbias, 'back')
wField.AddPlaneY(d, 0, 'front')
wField.AddStripOnPlaneY('z', d, -halfpitch, halfpitch, 'strip')
wField.AddPixelOnPlaneY(d, -halfpitch, halfpitch, 
                           -halfpitch, halfpitch, 'pixel')
wField.AddReadout('strip');
wField.AddReadout('pixel');
wField.AddReadout('front');

# Create a sensor. 
sensor = ROOT.Garfield.Sensor()
sensor.AddComponent(uniformField)
label = 'strip'
sensor.AddElectrode(wField, label)

# Plot the weighting potential.
fieldView = ROOT.Garfield.ViewField()
cF = ROOT.TCanvas('cF', '', 600, 600)
fieldView.SetCanvas(cF)
fieldView.SetComponent(wField)
fieldView.SetArea(-0.5 * d, 0, 0.5 * d, d)
fieldView.PlotContourWeightingField('strip', 'v')

# Set the time bins.
nTimeBins = 1000
tmin =  0.
tmax = 10.
tstep = (tmax - tmin) / nTimeBins
sensor.SetTimeWindow(tmin, tstep, nTimeBins)

# Set up Heed.
track = ROOT.Garfield.TrackHeed()
track.SetSensor(sensor)
# Set the particle type and momentum [eV/c].
track.SetParticle('pion')
track.SetMomentum(180.e9)

# Simulate electron/hole drift lines using MC integration.
drift = ROOT.Garfield.AvalancheMC()
drift.SetSensor(sensor)
# Use steps of 1 micron.
drift.SetDistanceSteps(1.e-4)
drift.EnableSignalCalculation()

# Plot the signal if requested.
signalView = ROOT.Garfield.ViewSignal()
signalView.SetSensor(sensor)
cS = ROOT.TCanvas('cS', '', 600, 600)
signalView.SetCanvas(cS)
signalView.SetRangeX(0., 6.)
plotSignal = True

driftView = ROOT.Garfield.ViewDrift()
driftView.SetArea(-0.5 * d, 0, -0.5 * d, 0.5 * d, d, 0.5 * d)
cD = ROOT.TCanvas('cD', '', 600, 600)
driftView.SetCanvas(cD)
plotDrift = True
if plotDrift:
  track.EnablePlotting(driftView)
  drift.EnablePlotting(driftView)

# Flag to randomise the position of the track.  
smearx = True
nEvents = 10
for i in range(nEvents):
  print i, '/', nEvents
  if plotDrift: driftView.Clear()
  # Simulate a charged-particle track.
  xt = 0.;
  if smearx: xt = -0.5 * pitch + ROOT.Garfield.RndmUniform() * pitch
  track.NewTrack(xt, 0, 0, 0, 0, 1, 0)
  xc = ctypes.c_double(0.)
  yc = ctypes.c_double(0.)
  zc = ctypes.c_double(0.)
  tc = ctypes.c_double(0.)
  ec = ctypes.c_double(0.)
  extra = ctypes.c_double(0.)
  ne = ctypes.c_int(0)
  # Retrieve the clusters along the track.
  while track.GetCluster(xc, yc, zc, tc, ne, ec, extra):
    # Loop over the electrons in the cluster.
    for j in range(ne.value):
      xe = ctypes.c_double(0.)
      ye = ctypes.c_double(0.)
      ze = ctypes.c_double(0.)
      te = ctypes.c_double(0.)
      ee = ctypes.c_double(0.)
      dx = ctypes.c_double(0.)
      dy = ctypes.c_double(0.)
      dz = ctypes.c_double(0.)
      track.GetElectron(j, xe, ye, ze, te, ee, dx, dy, dz)
      # Simulate the electron and hole drift lines.
      if plotDrift:
        drift.DisablePlotting()
        if ROOT.Garfield.RndmUniform() < 0.01:
          drift.EnablePlotting(driftView)
      drift.DriftElectron(xe.value, ye.value, ze.value, te.value)
      drift.DriftHole(xe.value, ye.value, ze.value, te.value)
  if plotSignal:
    signalView.PlotSignal(label)
    signalView.GetCanvas().Update()
    ROOT.gSystem.ProcessEvents()
  if plotDrift:
    driftView.Plot(True)
    driftView.GetCanvas().Update()
    ROOT.gSystem.ProcessEvents()
