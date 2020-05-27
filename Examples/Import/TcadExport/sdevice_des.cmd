File {
  	Grid    = "../../Meshes/2D/2D_3pixel_msh.tdr"
	Parameter= "sdevice.par"
	Plot=   "@tdrdat@"
	Current="@plot@"
	Output= "@log@"
}

Electrode {
	{Name="pix1_electrode"	Voltage=0.0 Material = "Aluminum"}
	{Name="pix2_electrode" Voltage= 0.0 Material = "Aluminum"}
	{Name="pix3_electrode"	Voltage=0.0 Material = "Aluminum"}
	{Name="bot_electrode"	Voltage=0.0 Material = "Aluminum"}
}

Physics	{
	Fermi
	Temperature = @Temperature@
	Mobility(
	 	eHighFieldSaturation
	 	hHighFieldSaturation
	 	PhuMob( Phosphorus Klaassen )
	 )
	Recombination(
			SRH(
				DopingDependence 		
				TempDependence		
				ElectricField(Lifetime=Hurkx DensityCorrection=none )
			)
	 		eAvalanche (vanOverstraeten Eparallel)
			hAvalanche (vanOverstraeten Eparallel)  
		)	
	EffectiveIntrinsicDensity(BandGapNarrowing(Slotboom))
}

Physics(MaterialInterface="Oxide/Silicon") {
		Charge(Conc = 1e12 ) 					
}
#if @fluence@ != 0
Physics (material = "Silicon"){
	Traps(
		(Donor Level
		
		fromValBand
		Add2TotalDoping
		Conc = @<4*fluence>@
		EnergyMid = 0.48
		eXsection = 2.0e-14
		hXsection = 1e-14)
	
		(Acceptor Level
		
		fromCondBand
		Add2TotalDoping
		Conc = @<0.75*fluence>@
		EnergyMid = 0.525
		eXsection = 5e-15
		hXsection = 1e-14)

		(Acceptor Level
		
		fromValBand
		Add2TotalDoping
		Conc = @<36*fluence>@
		EnergyMid = 0.90
		eXsection = 1e-16
		hXsection = 1e-16)
	)
}
#endif
Plot {
	Potential	ElectricField/Vector eMobility hMobility 
	eLifetime hLifetime hDriftVelocity/Vector hDriftVelocity/Vector
	

}

Math{
	Method = blocked
	SubMethod = pardiso
	NumberOfThreads = maximum
	Digits=5
	Extrapolate
	RelErrControl
	Derivatives
	Notdamped = 50
	Iterations = 15
}

Solve {
	Coupled(iterations = 600){Poisson}
	Coupled(iterations = 600){Poisson Electron Hole}
	Quasistationary(
		InitialStep=1e-3
		Maxstep=0.04
		MinStep=1e-5
		Increment=1.4
		Decrement=2
		Goal{ name="bot_electrode" voltage= @voltage@}
	){ 
		Coupled {Poisson Electron Hole}
		Plot ( FilePrefix = "n@node@_0V" Time = (1.0) NoOverwrite ) 
	}
	Quasistationary(
		InitialStep=1e-3
		Maxstep=0.04
		MinStep=1e-7
		Increment=1.4
		Decrement=2
		Goal{ name="pix2_electrode" voltage= 1.0 }
	){ 
		Coupled {Poisson Electron Hole}
		Plot ( FilePrefix = "n@node@_1V" Time = (1.0) NoOverwrite ) 
	}
	
}
