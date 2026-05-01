using BeamTracking: C_LIGHT, Species, chargeof, massof, gyromagnetic_anomaly, 
E_to_R, E_to_v, implicit_integrator!, order_eight_integrator! 

@testset "Implicit" begin
  sp = Species("electron")
  q = chargeof(sp)
  mc2 = massof(sp)
  a = gyromagnetic_anomaly(sp)
  E_ref = 18e9
  beta_0 = E_to_v(sp, E_ref)/C_LIGHT
  tilde_m = sqrt(1 - beta_0^2)/beta_0
  p_over_q_ref = E_to_R(sp, E_ref)
  radiation_params = (q, mc2, E_ref)
  w_id = (1.0, 0.0, 0.0, 0.0)

  function sol(x, y, s, t, Ksol)
    potential = (0.0, -Ksol*y/2, Ksol*x/2, 0.0)
    jac = (0.0,    0.0,     0.0, 0.0,
           0.0,    -Ksol/2, 0.0, 0.0,
           Ksol/2, 0.0,     0.0, 0.0,
           0.0,    0.0,     0.0, 0.0)
    return potential, jac
  end
  b0 = Bunch([0.01 0.02 0.03 0.04 0.05 0.06], [1.0 0.0 0.0 0.0])
  params = (radiation_params, beta_0, tilde_m, a, 0.0, w_id, w_id, sol, 0.3, p_over_q_ref, Val{true}())
  order_eight_integrator!(1, b0.coords, implicit_integrator!, params, nothing, 1.2/3, 3, nothing, Val{false}(), Val{false}(), 1.2)
  v_expected = [0.04457384910680009 0.02571695892367658 0.06811660659555369 0.034813336046770935 0.04888640130696242 0.05998151571786634]
  q_expected = [0.9342295999456364 0.2091256485383741 0.23055437398128034 -0.17414418838112483]
  @test b0.coords.v ≈ v_expected
  @test b0.coords.q ≈ q_expected (rtol = 5e-8)

  function dipole(x, y, s, t, p)
    g, Kn0 = p
    potential = (0.0, 0.0, 0.0, -Kn0*(x + g*x^2/2))
    jac = (0.0,            0.0, 0.0, 0.0,
           0.0,            0.0, 0.0, 0.0,
           0.0,            0.0, 0.0, 0.0,
           -Kn0*(1 + g*x), 0.0, 0.0, 0.0)
    return potential, jac
  end
  b0 = Bunch([0.01 0.02 0.03 0.04 0.05 0.06], [1.0 0.0 0.0 0.0])
  params = (radiation_params, beta_0, tilde_m, a, -0.2, w_id, w_id, dipole, (-0.2, 0.1), p_over_q_ref, Val{true}())
  order_eight_integrator!(1, b0.coords, implicit_integrator!, params, nothing, 1.2/3, 3, nothing, Val{false}(), Val{false}(), 1.2)
  v_expected = [-0.18731064753292975 -0.350614104089283 0.07662549801271587 0.03995709101566334 0.014424302720948396 0.05886291189268685]
  q_expected = [-0.9040459577632979 -0.0034797733367084084 -0.4268728902595191 0.021641002538076357]
  @test b0.coords.v ≈ v_expected
  @test b0.coords.q ≈ q_expected (rtol = 3e-7)

  function oct(x, y, s, t, Kn3)
    potential = (0.0, 0.0, 0.0, -Kn3/24*(x^4 - 6*x^2*y^2 + y^4))
    jac = (0.0,            0.0, 0.0, 0.0,
           0.0,            0.0, 0.0, 0.0,
           0.0,            0.0, 0.0, 0.0,
           -Kn3/6*(x^3 - 3*x*y^2), -Kn3/6*y*(y^2 - 3*x^2), 0.0, 0.0)
    return potential, jac
  end
  b0 = Bunch([0.01 0.02 0.03 0.04 0.05 0.06], [1.0 0.0 0.0 0.0])
  params = (radiation_params, beta_0, tilde_m, a, 0.0, w_id, w_id, oct, 10.1, p_over_q_ref, Val{true}())
  order_eight_integrator!(1, b0.coords, implicit_integrator!, params, nothing, 1.2/3, 3, nothing, Val{false}(), Val{false}(), 1.2)
  v_expected = [0.032802114477371816 0.020410479075338005 0.07525971038509004 0.03983249578484203 0.048930310980042475 0.05999997774858176]
  q_expected = [0.9999570879314574 0.003495937725397606 0.008570532232252856 -0.00038300428514522394]
  @test b0.coords.v ≈ v_expected
  @test b0.coords.q ≈ q_expected 

  function crazy(x, y, s, t, p)
    a0, c0 = p
    potential = a0 .* (C_LIGHT*sin(x)*cos(y)*sinh(s)*cosh(t), cosh(x)*sin(y)*cos(s)*sinh(t), sinh(x)*cosh(y)*sin(s)*cos(t), c0*cos(x)*sinh(y)*cosh(s)*sin(t))
    jac =  a0 .* (C_LIGHT*cos(x)*cos(y)*cosh(t)*sinh(s), -C_LIGHT*cosh(t)*sin(x)*sin(y)*sinh(s), C_LIGHT*cos(y)*cosh(s)*cosh(t)*sin(x), C_LIGHT*cos(y)*sin(x)*sinh(s)*sinh(t),
           cos(s)*sin(y)*sinh(t)*sinh(x),  cos(s)*cos(y)*cosh(x)*sinh(t), -cosh(x)*sin(s)*sin(y)*sinh(t),  cos(s)*cosh(t)*cosh(x)*sin(y),
           cos(t)*cosh(x)*cosh(y)*sin(s),  cos(t)*sin(s)*sinh(x)*sinh(y),  cos(s)*cos(t)*cosh(y)*sinh(x), -cosh(y)*sin(s)*sin(t)*sinh(x),
          -c0*cosh(s)*sin(t)*sin(x)*sinh(y), c0*cos(x)*cosh(s)*cosh(y)*sin(t), c0*cos(x)*sin(t)*sinh(s)*sinh(y), c0*cos(t)*cos(x)*cosh(s)*sinh(y))
    return potential, jac
  end
  b0 = Bunch([0.01 0.02 0.03 0.04 0.05 0.06], [1.0 0.0 0.0 0.0])
  params = (radiation_params, beta_0, tilde_m, a, 0.0, w_id, w_id, crazy, (1.0, 1.0), p_over_q_ref, Val{true}())
  order_eight_integrator!(1, b0.coords, implicit_integrator!, params, nothing, 1.2/60, 60, nothing, Val{false}(), Val{false}(), 1.2)
  v_expected = [-0.24794117122676387 -0.7621581346830383 0.12143614880778619 0.024481117969862107 -0.006481844718705498 -0.0003806157078307523]
  q_expected = [-0.8671095776608468 -0.16135398678020024 -0.405824923995721 0.23956627964481797]
  @test b0.coords.v ≈ v_expected
  @test b0.coords.q ≈ q_expected 

  b0 = Bunch([0.01 0.02 0.03 0.04 0.05 0.06] .+ collect(transpose(@vars(D1))), TPS64{D1}[1 0 0 0])
  args = (implicit_integrator!, params, nothing, 1.2/60, 60, nothing, Val{false}(), Val{false}(), 1.2)
  order_eight_integrator!(1, b0.coords, args...)
  @test scalar.(b0.coords.v) ≈ v_expected
  @test scalar.(b0.coords.q) ≈ q_expected 
  M = [0.6236345934584919    0.9495961035774303   0.017035902012427193 0.2187104923599301  -2.73297873027427e-10   0.21831631202222596; 
      -0.4843268398176658    0.7981101244239478   0.0624521841956804   0.5724693392197097  -1.5970299450143453e-9  0.02736760564700434; 
      -0.436011338811749    -0.3594099028586768   0.9890748851735769   1.0691744171696649  -2.3407456473471284e-9 -0.12579170052828884; 
       0.11654780028265424   0.10595575615675462 -0.08783741159302688  0.8816704143893558  -4.194902351203336e-9   0.01412623663158183; 
      -0.11945699132057884   0.17613239866393063  0.008157852668804522 0.005162919229188737 1.0000000000830762     0.09432913318067639; 
      -0.010855520486550325 -0.03388639005364529  0.013570539648297493 0.09507272679052414 -1.3953727419977462e-10 0.87897338710916]
  Q = [-4.003668537643902  -4.323651212047558   0.5348998994545897   4.717968619077585  -2.49853749526394e-9  -0.9964924657421544; 
       -0.34214759237863596 1.7511043319960213 -0.08354112709849622 -1.0771679496754156 -2.8648881796123093e-9 0.4101322809068979; 
        11.59702678163501   9.644617457227584  -1.2422798535283495  -10.309730971563797  1.4248220979198173e-8 2.0405723195357095; 
        4.923632647634007   1.8679360247524797 -0.22461842867509274 -1.1134931264839765  1.3163439971223364e-8 0.1261589669099078]
  @test GTPSA.jacobian(b0.coords.v) ≈ M
  @test GTPSA.jacobian(b0.coords.q) ≈ Q

  function track(v)
    b0 = Bunch(v')
    order_eight_integrator!(1, b0.coords, args...)
    return b0.coords.v'
  end
  @test ForwardDiff.jacobian(track, [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]) ≈ M

  M = [0.5642361394656592    0.9695143326521845   0.0008233255950771141 0.23459732508210268 -3.98096261096381e-10  0.26202087937470775; 
      -0.543530422071873     0.762355986589084    0.005384939535613539  0.5568726180598322  -1.5243574558899086e-9 0.036470556982992404; 
      -0.4593283172457024   -0.3984638458095007   0.9811769291798218    1.100107521964648   -2.298597550516041e-9 -0.12126875320594324; 
       0.08747063284126881   0.07920527444449403 -0.1220421218444098    0.8468859184498521  -4.040295096962048e-9  0.011352492304963603; 
      -0.1600021639255181    0.19633471917846448  0.003002641020079756  0.03951167031586203 1.0000000000195635     0.1128782075288401; 
      -0.012282105647725347 -0.0399005425614098   0.008678224403495206  0.09108613577838762 -9.043193961707872e-11 0.8815006978818158]
  #test_matrix(M, KernelCall(order_eight_integrator!, args); type_stable=true, no_scalar_allocs=true)

  a0 = 1.0 + (@params D1_1)[1]
  b0 = Bunch([0.01 0.02 0.03 0.04 0.05 0.06] .+ collect(transpose(@vars(D1_1))), TPS64{D1_1}[1 0 0 0])
  params = (radiation_params, beta_0, tilde_m, a, 0.0, w_id, w_id, crazy, (a0, 1.0), p_over_q_ref, Val{true}())
  order_eight_integrator!(1, b0.coords, implicit_integrator!, params, nothing, 1.2/60, 60, nothing, Val{false}(), Val{false}(), 1.2)
  @test scalar.(b0.coords.v) ≈ v_expected
  @test scalar.(b0.coords.q) ≈ q_expected 
  M = [0.6236345934584699    0.9495961035774574   0.017035902012431346 0.21871049235992143 -2.732978730272376e-10  0.21831631202219653  -0.26766551627369584; 
      -0.484326839817701     0.7981101244239058   0.062452184195696855 0.5724693392197072  -1.5970299450142425e-9  0.027367605646998156 -0.7690749351321912; 
      -0.4360113388116895   -0.35940990285867236  0.989074885173537    1.0691744171696107  -2.3407456473471214e-9 -0.12579170052835006   0.10154148911372539; 
       0.11654780028268388   0.10595575615679668 -0.0878374115930651   0.8816704143893476  -4.194902351203437e-9   0.014126236631594152 -0.04760209051623857; 
      -0.11945699132061917   0.17613239866389865  0.008157852668812367 0.005162919229213165 1.0000000000830866     0.09432913318073993  -0.10803633514243256; 
      -0.010855520486557295 -0.03388639005364263  0.013570539648297658 0.0950727267905285  -1.395372741998208e-10  0.8789733871091513   -0.10664797131666366]
  Q = [-4.003668911274518   -4.323651522672887   0.5348999386629655   4.7179689442268415 -2.498537952817113e-9 -0.9964925152787328 -5.153645361761044; 
       -0.3421475312171915   1.7511043071864707 -0.08354111705747419 -1.0771678773979212 -2.8648877574033663e-9 0.41013226757794063 1.2892905015542115; 
        11.597026589966077   9.644617287278072  -1.2422798341661323  -10.309730791491138  1.4248220656589593e-8 2.0405722532749664  11.406794337190332; 
        4.923632757147277    1.8679360966994514 -0.22461843686744848 -1.1134931787038658  1.3163440144406119e-8 0.12615898461801006 1.5378576321517365]
  @test GTPSA.jacobian(b0.coords.v, include_params=true) ≈ M
  @test GTPSA.jacobian(b0.coords.q, include_params=true) ≈ Q (rtol = 3e-8)

  a0 = 1.0 + (@params D1_2)[1]
  c0 = 1.0 + (@params D1_2)[2]
  b0 = Bunch([0.01 0.02 0.03 0.04 0.05 0.06] .+ collect(transpose(@vars(D1_2))), TPS64{D1_2}[1 0 0 0])
  params = (radiation_params, beta_0, tilde_m, a, 0.0, w_id, w_id, crazy, (a0, c0), p_over_q_ref, Val{true}())
  order_eight_integrator!(1, b0.coords, implicit_integrator!, params, nothing, 1.2/60, 60, nothing, Val{false}(), Val{false}(), 1.2)
  @test scalar.(b0.coords.v) ≈ v_expected
  @test scalar.(b0.coords.q) ≈ q_expected 
  M = [0.6236345934584648    0.9495961035774615   0.017035902012432144 0.21871049235991938 -2.732978730271998e-10 0.21831631202219096   -0.2676655162736905   7.257029615247699e-11; 
      -0.48432683981770763   0.7981101244238971   0.062452184195700186 0.5724693392197057  -1.5970299450142197e-9  0.02736760564699731  -0.7690749351321868   7.539602174052556e-10; 
      -0.43601133881167686  -0.3594099028586719   0.9890748851735325   1.0691744171695987  -2.3407456473471268e-9 -0.12579170052836247   0.10154148911373952  9.861736852807539e-10; 
       0.11654780028268955   0.10595575615680505 -0.08783741159307329  0.8816704143893457  -4.1949023512034535e-9  0.014126236631596725 -0.04760209051624329  2.9693961461517193e-9; 
      -0.11945699132062708   0.17613239866389174  0.00815785266881394  0.005162919229218142 1.0000000000830964     0.09432913318075303  -0.10803633514244579 -8.416261746589163e-11; 
      -0.010855520486558544 -0.03388639005364217  0.013570539648297774 0.0950727267905295  -1.3953727419983035e-10 0.8789733871091538   -0.106647971316665   -3.5488236986610555e-10]
  Q = [-4.003668825703249  -4.323651467126826   0.5348999292958801   4.717968868793816  -2.4985378311294204e-9 -0.9964924968020226 -5.153645274147697  -2.468861713639696e-9; 
       -0.34214754507748046 1.751104315348463  -0.08354111910494988 -1.0771678911737181 -2.864887853554742e-9   0.41013226812915626 1.2892905128065009 -1.1056102995505371e-9; 
        11.597026640077837  9.644617354640323  -1.2422798369951702  -10.30973082488686   1.4248220727836011e-8  2.0405722513247886  11.40679436251822   6.647715697066145e-9; 
        4.923632751344744   1.8679360850878741 -0.22461843487926048 -1.1134931676013486  1.3163440162595715e-8  0.12615897844622773 1.5378576115996576  1.5805257053628144e-9]
  @test GTPSA.jacobian(b0.coords.v, include_params=true) ≈ M
  @test GTPSA.jacobian(b0.coords.q, include_params=true) ≈ Q (rtol = 3e-8)

  function crazy_unnormalized(x, y, s, t, p)
    potential = (C_LIGHT*sin(x)*cos(y)*sinh(s)*cosh(t), cosh(x)*sin(y)*cos(s)*sinh(t), sinh(x)*cosh(y)*sin(s)*cos(t), cos(x)*sinh(y)*cosh(s)*sin(t))
    jac = (C_LIGHT*cos(x)*cos(y)*cosh(t)*sinh(s), -C_LIGHT*cosh(t)*sin(x)*sin(y)*sinh(s), C_LIGHT*cos(y)*cosh(s)*cosh(t)*sin(x), C_LIGHT*cos(y)*sin(x)*sinh(s)*sinh(t),
           cos(s)*sin(y)*sinh(t)*sinh(x),  cos(s)*cos(y)*cosh(x)*sinh(t), -cosh(x)*sin(s)*sin(y)*sinh(t),  cos(s)*cosh(t)*cosh(x)*sin(y),
           cos(t)*cosh(x)*cosh(y)*sin(s),  cos(t)*sin(s)*sinh(x)*sinh(y),  cos(s)*cos(t)*cosh(y)*sinh(x), -cosh(y)*sin(s)*sin(t)*sinh(x),
          -cosh(s)*sin(t)*sin(x)*sinh(y),  cos(x)*cosh(s)*cosh(y)*sin(t),  cos(x)*sin(t)*sinh(s)*sinh(y),  cos(t)*cos(x)*cosh(s)*sinh(y))
    return p_over_q_ref .* potential, p_over_q_ref .* jac
  end
  b0 = Bunch([0.01 0.02 0.03 0.04 0.05 0.06], [1.0 0.0 0.0 0.0])
  params = (radiation_params, beta_0, tilde_m, a, 0.0, w_id, w_id, crazy_unnormalized, nothing, p_over_q_ref, Val{false}())
  order_eight_integrator!(1, b0.coords, implicit_integrator!, params, nothing, 1.2/60, 60, nothing, Val{false}(), Val{false}(), 1.2)
  v_expected = [-0.24794117122676387 -0.7621581346830383 0.12143614880778619 0.024481117969862107 -0.006481844718705498 -0.0003806157078307523]
  q_expected = [-0.8671095776608468 -0.16135398678020024 -0.405824923995721 0.23956627964481797]
  @test b0.coords.v ≈ v_expected
  @test b0.coords.q ≈ q_expected 

end