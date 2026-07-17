# Exercises BeamTrackingPythonCallExt: constructing TimeDependentParam / BatchParam
# from Python objects. Requires PythonCall (which provisions Python via CondaPkg
# on first use).
@testset "BeamTrackingPythonCallExt" begin
  # TimeDependentParam from a Python callable f(t)
  tdp = TimeDependentParam(pyeval(Py, "lambda t: 2*t + 1", Main))
  @test tdp(3.0) == 7.0

  # constant TimeDependentParam from a Python value
  @test TimeDependentParam(pyeval(Py, "5.0", Main))(0.0) == 5.0

  # BatchParam from a Python sequence
  bp = BatchParam(pyeval(Py, "[0.1, 0.2, 0.3]", Main))
  @test collect(Float64, bp.batch) == [0.1, 0.2, 0.3]

  # BatchParam scalar (degenerate batch) from a Python value
  @test BatchParam(pyeval(Py, "0.5", Main)).batch == 0.5
end
