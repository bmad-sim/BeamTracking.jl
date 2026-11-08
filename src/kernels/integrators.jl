#
# ===============  I N T E G R A T O R S  ===============
#

# Tracks n_steps integration steps, each made of the substeps of lengths ws (the scheme of the
# integrator), with the fringes, radiation and callbacks. `ker` has a single call site, so it is
# inlined once: written out per substep it was inlined at every substep (up to 70 copies of the
# kernel for order 10), which multiplied the compile time and the size of the compiled code.
# `ker::K` forces specialization on the kernel (it is not called directly, only splatted)
function integrate_steps!(i, coords::Coords, ker::K, params, photon_params, ds_step, n_steps, edge_params, ::Val{fringe_in}, ::Val{fringe_out}, ws) where {K, fringe_in, fringe_out}
  @inbounds begin
    if !isnothing(edge_params) && fringe_in
      fringe!(i, coords, edge_params..., Val{1}())
    end
    s = zero(first(ws))
    if !isnothing(photon_params)
      stochastic_radiation!(i, coords, s, photon_params..., ds_step / 2)
    end
    for step in 1:n_steps
      for w in ws
        ker(i, coords, s, params..., w)
        s += w
      end
      if step < n_steps
        if !isnothing(photon_params)
          stochastic_radiation!(i, coords, s, photon_params..., ds_step)
        end
        dt_ref = compute_dt_ref(s, ker, params)
        execute_callbacks(i, coords, s, dt_ref)
      end
    end
    if !isnothing(photon_params)
      stochastic_radiation!(i, coords, s, photon_params..., ds_step / 2)
    end
    if !isnothing(edge_params) && fringe_out
      fringe!(i, coords, edge_params..., Val{-1}())
    end
  end
  return nothing
end

@inline function order_two_integrator!(i, coords::Coords, ker, params, photon_params, ds_step, n_steps, edge_params, fringe_in, fringe_out, use_optimized_schemes, L)
  return integrate_steps!(i, coords, ker, params, photon_params, ds_step, n_steps, edge_params, fringe_in, fringe_out, (ds_step,))
end

@inline @generated function order_four_integrator!(i, coords::Coords{<:Any,V}, ker, params, photon_params, ds_step, n_steps, edge_params, ::Val{fringe_in}, ::Val{fringe_out}, ::Val{use_optimized_schemes}, L) where {V, fringe_in, fringe_out, use_optimized_schemes}
  if use_optimized_schemes
    c = (
      -0.6579630871775028477799196480191312730312347412109375,
      0.414490771794375711944979912004782818257808685302734375,
    )
  else
    c = (
      -1.7024143839193153215916254339390434324741363525390625,
      1.3512071919596577718181151794851757586002349853515625,
    )
  end

  T = ForwardDiff.valtype(eltype(V))
  if T == Float16 || T == Float32
    c = T.(c)
  end
  return quote
    if use_optimized_schemes
      w0 = $(c[1]) * ds_step
      w1 = $(c[2]) * ds_step
      ws = (w1, w1, w0, w1, w1)
    else
      w0 = $(c[1]) * ds_step
      w1 = $(c[2]) * ds_step
      ws = (w1, w0, w1)
    end
    integrate_steps!(i, coords, ker, params, photon_params, ds_step, n_steps, edge_params, Val{fringe_in}(), Val{fringe_out}(), ws)
  end
end


@inline @generated function order_six_integrator!(i, coords::Coords{<:Any,V}, ker, params, photon_params, ds_step, n_steps, edge_params, ::Val{fringe_in}, ::Val{fringe_out}, ::Val{use_optimized_schemes}, L) where {V, fringe_in, fringe_out, use_optimized_schemes}
  if use_optimized_schemes
    c = (
      0.39907751534871587459988795520665,
      0.11805530653002387170273438954049,
      -0.35000324893920896516170830911323,
      0.12961893756907034772505366537091,
      0.13070531011449225190542755785015,
      0.13346562851074760407046858832209,
      0.13861930854051695245808013042625,
    )
  else
    c = (
      1.315186320683911169737712043570355,
      -1.17767998417887100694641568096432,
      0.235573213359358133684793182978535,
      0.784513610477557263819497633866351,
      NaN,
      NaN,
      NaN,
    )
  end

  T = ForwardDiff.valtype(eltype(V))
  if T == Float16 || T == Float32
    c = T.(c)
  end
  return quote
    if use_optimized_schemes
      w0 = $(c[1]) * ds_step
      w1 = $(c[2]) * ds_step
      w2 = $(c[3]) * ds_step
      w3 = $(c[4]) * ds_step
      w4 = $(c[5]) * ds_step
      w5 = $(c[6]) * ds_step
      w6 = $(c[7]) * ds_step
      ws = (w6, w5, w4, w3, w2, w1, w0, w1, w2, w3, w4, w5, w6)
    else
      w0 = $(c[1]) * ds_step
      w1 = $(c[2]) * ds_step
      w2 = $(c[3]) * ds_step
      w3 = $(c[4]) * ds_step
      ws = (w3, w2, w1, w0, w1, w2, w3)
    end
    integrate_steps!(i, coords, ker, params, photon_params, ds_step, n_steps, edge_params, Val{fringe_in}(), Val{fringe_out}(), ws)
  end
end


@inline @generated function order_eight_integrator!(i, coords::Coords{<:Any,V}, ker, params, photon_params, ds_step, n_steps, edge_params, ::Val{fringe_in}, ::Val{fringe_out}, ::Val{use_optimized_schemes}, L) where {V, fringe_in, fringe_out, use_optimized_schemes}
  if use_optimized_schemes
    c = (
      -0.38263596012643665350944670744040,
      0.11699135019217642180722881433533,
      0.12581718736176041804392391641587,
      0.12603912321825988140305670268365,
      0.11892905625000350062692972283951,
      0.11317848435755633314700952515599,
      -0.24445266791528841269462171413216,
      -0.23341414023165082198780281128319,
      0.35337821052654342419534541324080,
      0.10837408645835726397433410591546,
      0.10647728984550031823931967854896,
    )
  else
    c = (
        1.7084530707869978,
        0.102799849391985,
        -1.96061023297549,
        1.93813913762276,
        -0.158240635368243,
        -1.44485223686048,
        0.253693336566229,
        0.914844246229740,
        NaN,
        NaN,
        NaN
    )
  end

  T = ForwardDiff.valtype(eltype(V))
  if T == Float16 || T == Float32
    c = T.(c)
  end
  return quote
    if use_optimized_schemes
      w0 = $(c[1]) * ds_step
      w1 = $(c[2]) * ds_step
      w2 = $(c[3]) * ds_step
      w3 = $(c[4]) * ds_step
      w4 = $(c[5]) * ds_step
      w5 = $(c[6]) * ds_step
      w6 = $(c[7]) * ds_step
      w7 = $(c[8]) * ds_step
      w8 = $(c[9]) * ds_step
      w9 = $(c[10]) * ds_step
      w10 = $(c[11]) * ds_step
      ws = (w10, w9, w8, w7, w6, w5, w4, w3, w2, w1, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10)
    else
      w0 = $(c[1]) * ds_step
      w1 = $(c[2]) * ds_step
      w2 = $(c[3]) * ds_step
      w3 = $(c[4]) * ds_step
      w4 = $(c[5]) * ds_step
      w5 = $(c[6]) * ds_step
      w6 = $(c[7]) * ds_step
      w7 = $(c[8]) * ds_step
      ws = (w7, w6, w5, w4, w3, w2, w1, w0, w1, w2, w3, w4, w5, w6, w7)
    end
    integrate_steps!(i, coords, ker, params, photon_params, ds_step, n_steps, edge_params, Val{fringe_in}(), Val{fringe_out}(), ws)
  end
end


@inline @generated function order_ten_integrator!(i, coords::Coords{<:Any,V}, ker, params, photon_params, ds_step, n_steps, edge_params, ::Val{fringe_in}, ::Val{fringe_out}, use_optimized_schemes, L) where {V, fringe_in, fringe_out}
  c = (
    0.049317735759594537917680008339338 ,
    0.049674370639729879054568800279461 ,
    0.050665090759924496335874344156866 ,
    0.051942502962449647037182904015976 ,
    -0.39203335370863990644808193642610  ,
    -0.0048663605831352617621956593099771,
    0.41143087395589023782070411897608  ,
    0.10308739852747107731580277001372  ,
    -0.39910563013603589787862981058340  ,
    0.36613344954622675119314812353150  ,
    0.11199342399981020488957508073640  ,
    0.074973343155891435666137105641410 ,
    -0.26973340565451071434460973222411  ,
    0.13096206107716486317465685927961  ,
    -0.22959284159390709415121339679655  ,
    0.027918383235078066109520273275299 ,
    0.31309610341510852776481247192647  ,
    0.078795722521686419263907679337684 ,
  )
  T = ForwardDiff.valtype(eltype(V))
  if T == Float16 || T == Float32
    c = T.(c)
  end
  return quote
    w0 = $(c[1]) * ds_step
    w1 = $(c[2]) * ds_step
    w2 = $(c[3]) * ds_step
    w3 = $(c[4]) * ds_step
    w4 = $(c[5]) * ds_step
    w5 = $(c[6]) * ds_step
    w6 = $(c[7]) * ds_step
    w7 = $(c[8]) * ds_step
    w8 = $(c[9]) * ds_step
    w9 = $(c[10]) * ds_step
    w10 = $(c[11]) * ds_step
    w11 = $(c[12]) * ds_step
    w12 = $(c[13]) * ds_step
    w13 = $(c[14]) * ds_step
    w14 = $(c[15]) * ds_step
    w15 = $(c[16]) * ds_step
    w16 = $(c[17]) * ds_step
    w17 = $(c[18]) * ds_step
    ws = (w17, w16, w15, w14, w13, w12, w11, w10, w9, w8, w7, w6, w5, w4, w3, w2, w1, w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15, w16, w17)
    integrate_steps!(i, coords, ker, params, photon_params, ds_step, n_steps, edge_params, Val{fringe_in}(), Val{fringe_out}(), ws)
  end
end
