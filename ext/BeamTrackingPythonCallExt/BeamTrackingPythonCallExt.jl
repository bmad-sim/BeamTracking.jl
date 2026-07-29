module BeamTrackingPythonCallExt
using PythonCall
import BeamTracking: TimeDependentParam, BatchParam

# Construction of BeamTracking's deferred-parameter types from Python objects.
#
# `TimeDependentParam` and `BatchParam` are not `Number` subtypes, so juliacall
# neither overloads their operators nor converts Python callables/sequences into
# them automatically. These methods let a Python callable become a
# time-dependent parameter, and a Python sequence become a batch parameter,
# converting values to `T` (default `Float64`) at evaluation time.
#
# Operator overloading on the *Python* side (so `2 * Time() + 1` works from
# Python) is provided by the operators these Julia types already define; a Python
# wrapper only needs to forward `+ - * / ^` onto them.

_pycallable(f::Py) = pytruth(pybuiltins.callable(f))

# --- TimeDependentParam ------------------------------------------------------

function TimeDependentParam(f::Py, ::Type{T}) where {T}
  if _pycallable(f)
    return TimeDependentParam(t -> pyconvert(T, f(t)), false)
  else
    return TimeDependentParam(pyconvert(T, f))  # a constant
  end
end

TimeDependentParam(f::Py) = TimeDependentParam(f, Float64)

# --- BatchParam --------------------------------------------------------------

function BatchParam(batch::Py, ::Type{T}) where {T}
  if pyhasattr(batch, "__len__")           # a sequence (list/tuple/ndarray/...)
    return BatchParam(pyconvert(Vector{T}, batch))
  else
    return BatchParam(pyconvert(T, batch))  # a scalar (degenerate batch)
  end
end

BatchParam(batch::Py) = BatchParam(batch, Float64)

end
