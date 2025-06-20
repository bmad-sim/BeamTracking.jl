using Beamlines

# quad1 =  Quadrupole( L = 1.8, K1 = -0.2291420342)
drift1 = Drift(L = 0.5)
# quad2 =  Quadrupole( L = 1.4, K1 = 0.2267785688)
# drift2 =  Drift( L = 0.5)
# sbend1 =  SBend( L = 5.50007539103, g = -3.2977170394029E-3, e1 = -9.0688461675E-3, e2 = -9.0688461675E-3)

ring = Beamline([drift1], Brho_ref=59.52872449027632)