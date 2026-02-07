#---------------------------------------------------------------------------------------------------

#@inline 
function RFcavity(tm, bunch, bmultipoleparams, rfparams, beamlineparams, L)
    rf_omega = rf_omega_calc(rfparams, beamlineparams.beamline.line[end].s_downstream, bunch.species, bunch.p_over_q_ref)
    t_phi0 = rf_phi0_calc_old(rfparams, beamlineparams.beamline.species_ref) / rf_omega

    if !isactive(bmultipoleparams)
      return pure_rf(tm, bunch, rfparams, rf_omega, t_phi0, L)
    else
      return bmultipole_rf(tm, bunch, bmultipoleparams, rfparams, rf_omega, t_phi0, L)
    end
end