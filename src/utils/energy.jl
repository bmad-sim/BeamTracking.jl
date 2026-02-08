# Particle energy conversions =============================================================
E_to_R(species::Species, E) = @FastGTPSA massof(species)*sinh(acosh(E/massof(species)))/C_LIGHT/chargeof(species) 
E_to_c_beta(species::Species, E) = @FastGTPSA sqrt(E^2 - massof(species)^2) * C_LIGHT / E

R_to_E(species::Species, R) = @FastGTPSA sqrt((R*C_LIGHT*chargeof(species))^2 + massof(species)^2)
R_to_gamma(species::Species, R) = @FastGTPSA sqrt((R*C_LIGHT/massof(species))^2+1)
R_to_pc(species::Species, R) = @FastGTPSA R*chargeof(species)*C_LIGHT
R_to_beta_gamma(species::Species, R) = @FastGTPSA R*chargeof(species)*C_LIGHT/massof(species)
R_to_v(species::Species, R) = @FastGTPSA chargeof(species)*C_LIGHT / sqrt(1+(massof(species)/(R*C_LIGHT))^2)

beta_gamma_to_v(beta_gamma) = @FastGTPSA C_LIGHT*beta_gamma/sqrt(1+beta_gamma^2)
pc_to_R(species::Species, pc) = @FastGTPSA pc/C_LIGHT/chargeof(species)
