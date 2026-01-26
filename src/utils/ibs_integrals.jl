"""
Computes the Carlson symmetric form R_F(x,y,z) with error bounded by r.
The algorithm is from arXiv:math/9409227v1.
"""
function Carlson_RF(x, y, z; r=1e-11)
  A0 = (x+y+z)/3
  Q = (3*r)^(-1/6)*max(abs(A0-x), abs(A0-y), abs(A0-z))

  Am = A0
  lambdam = 0
  four_m = 1
  xm, ym, zm = x, y, z
  while Q >= four_m*abs(Am)
    xm_s, ym_s, zm_s = sqrt(xm), sqrt(ym), sqrt(zm)
    lambdam = xm_s*ym_s + xm_s*zm_s + ym_s*zm_s
    four_m *= 4
    Am = (Am + lambdam)/4
    xm = (xm + lambdam)/4
    ym = (ym + lambdam)/4
    zm = (zm + lambdam)/4
  end
  X = (A0-x)/(four_m*Am)
  Y = (A0-y)/(four_m*Am)
  Z = -X-Y
  E2 = X*Y - Z^2
  E3 = X*Y*Z
  result = 1 - (E2/10) + (E3/14) + (E2^2/24) - (3*E2*E3/44)
  result /= sqrt(Am)
  return result
end


"""
Computes the Carlson symmetric form R_D(x,y,z) with error bounded by r.
The algorithm is from arXiv:math/9409227v1.
"""
function Carlson_RD(x, y, z; r=1e-11)
  A0 = (x+y+3*z)/5
  Q = (r/4)^(-1/6)*max(abs(A0-x), abs(A0-y), abs(A0-z))

  Am = A0
  four_m = 1
  xm, ym, zm = x, y, z
  lambdam = 0
  sum = 0
  while Q >= four_m*abs(Am)
    xm_s, ym_s, zm_s = sqrt(xm), sqrt(ym), sqrt(zm)
    lambdam = xm_s*ym_s + xm_s*zm_s + ym_s*zm_s
    sum += 3/(four_m*zm_s*(zm+lambdam))
    four_m *= 4
    Am = (Am + lambdam)/4
    xm = (xm + lambdam)/4
    ym = (ym + lambdam)/4
    zm = (zm + lambdam)/4
  end
  X = (A0-x)/(four_m*Am)
  Y = (A0-y)/(four_m*Am)
  Z = -(X+Y)/3
  E2 = X*Y - 6*Z^2
  E3 = (3*X*Y - 8Z^2)*Z
  E4 = 3*(X*Y - Z^2)*Z^2
  E5 = X*Y*Z^3
  result = 1 - (3/14)*E2 + E3/6 + (9/88)*E2^2 -(3/22)*E4 - (9/52)*E2*E3 + (3/26)*E5
  result /= four_m*Am^(3/2)
  result += sum
  return result
end


"""
Computes the three integrals used in the computation of IBS kicks.
"""
function IBS_integrals(λ1, λ2, λ3)
  if λ1 ≈ λ2 ≈ λ3
    I = 4*pi/(3*λ1)
    return I, I, I
  elseif λ1 ≈ λ2
    ratio = λ3/λ1
    value = ratio - 1
    if value >= 0
      sq = sqrt(value)
      arctan = sq*atan(sq)
    else
      sq = sqrt(-value)
      arctan = -sq*atanh(sq)
    end
    I1 = 2*pi/(λ1-λ3)*(1+λ3/(λ1-λ3)*arctan)
    R_D = Carlson_RD(ratio, ratio, 1)
    I3 = 4*pi/(3*λ1)*R_D
    return I1, I1, I3
  else
    args = (λ3/λ1, λ3/λ2, 1)
    R_F = Carlson_RF(args...)
    R_D = Carlson_RD(args...)
    I1 = 4*pi/(λ2-λ1)*(sqrt(λ2/λ1)*R_F-1/3*sqrt(λ2/λ1)*(1-λ3/λ2)*R_D-1)
    I2 = 4*pi/(λ1-λ2)*(sqrt(λ1/λ2)*R_F-1/3*sqrt(λ1/λ2)*(1-λ3/λ1)*R_D-1)
    I3 = 4*pi/(3*sqrt(λ1*λ2))*R_D
    return I1, I2, I3
  end
end