package gonumeth

import "math"

const hsolve float64 = 0.01

// TODO: Finish caching implementation

type CachedSingleVarFunction SingleVarFunction

func CacheFunction(f SingleVarFunction) (fres CachedSingleVarFunction) {
	m := make(map[float64]float64)
	return func(x float64) float64 {
		var (
			res float64
			ok  bool
		)
		if res, ok = m[x]; ok {
			return res
		}
		res = f(x)
		m[x] = res
		return res
	}
}

// A common type that all simple solvers implement
type NSimpleSolver func(f SingleVarFunction, x0 float64, maxIterations int,
	epsilon float64) float64

// The iteration function for the bisection method
func bisectionIteration(f SingleVarFunction, xi float64,
	xi_1 float64) (xi_new float64, xi1 float64) {
	xi1 = 0.5 * (xi + xi_1)
	fi1 := f(xi1)
	fi := f(xi)
	if fi1*fi <= 0 {
		xi_new = xi
	} else {
		xi_new = xi_1
	}
	return
}

// NSimpleSolveBisection attempts to find a root of the function f starting
// at x0 by using the Bisection method.
// A `root` value of NaN means the function failed.
func NSimpleSolveBisection(f SingleVarFunction, x0 float64, maxIterations int,
	epsilon float64) (result float64) {
	step := hsolve
	xi := x0 + step
	xi_1 := x0
	fi := f(xi_1)
	var i int = 0
	for ; i < int(float64(maxIterations)*0.25) || maxIterations == 0; i++ {
		if fi*f(xi) <= 0 {
			break
		}
		step += hsolve
		xi = x0 + step
	}
	step = hsolve
	for ; i < int(float64(maxIterations)*0.35) || maxIterations == 0; i++ {
		if fi*f(xi) <= 0 {
			break
		}
		step *= 2
		xi = x0 + step
	}
	for ; i < maxIterations || maxIterations == 0; i++ {
		fi = f(0.5 * (xi + xi_1))
		if math.Abs(fi) < epsilon {
			return 0.5 * (xi + xi_1)
		}
		xi_1, xi = bisectionIteration(f, xi, xi_1)
	}
	return math.NaN()
}

// The iteration function for the Newton method
func newtonIteration(f SingleVarFunction, xi float64) float64 {
	deriv := NDifferentiateCentral(f, xi, hsolve)
	if deriv == 0 {
		return math.NaN()
	}
	fi := f(xi)
	return xi - fi/deriv
}

// NSimpleSolveNewton attempts to find a root of the function f starting
// at x0 using the Newton method. This method generally requires a "good"
// behavior of f locally (it being differentiable and having non-zero deriv.).
// If this is not true, the function might fail to produce the proper result.
// A `root` value of NaN means the function failed.
func NSimpleSolveNewton(f SingleVarFunction, x0 float64, maxIterations int,
	epsilon float64) (result float64) {
	var (
		xi float64 = x0
		fi float64
	)
	for i := 0; i < maxIterations || maxIterations == 0; i++ {
		if math.IsNaN(xi) {
			return xi
		}
		fi = f(xi)
		if math.Abs(fi) < epsilon {
			return xi
		}
		xi = newtonIteration(f, xi)
	}
	return math.NaN()
}

// Iteration function for the Halley method
func halleyIteration(f SingleVarFunction, fprime SingleVarFunction,
	xi float64) float64 {
	fi := f(xi)
	f_di := fprime(xi)
	f_d2i := NDifferentiateCentral(fprime, xi, hsolve)
	factor := 2*f_di*f_di - fi*f_d2i
	if factor == 0 {
		return math.NaN()
	}
	return xi - 2*fi*f_di/factor
}

// NSimpleSolveHalley attempts to find a root of the function f starting
// at x0 using the Halley method. This method is the (potentially) fastest
// method provided in this package, however it requires "very good" behavior
// of f locally (in particular it being differentiable and a non-zero value
// of the function and its first derivative).
// A `root` value of NaN means the function failed.
func NSimpleSolveHalley(f SingleVarFunction, x0 float64, maxIterations int,
	epsilon float64) (root float64) {
	var (
		xi     float64 = x0
		fi     float64
		fderiv SingleVarFunction = NDerivative(f, hsolve)
	)
	for i := 0; i < maxIterations || maxIterations == 0; i++ {
		if math.IsNaN(xi) {
			return xi
		}
		fi = f(xi)
		if math.Abs(fi) < epsilon {
			return xi
		}
		xi = halleyIteration(f, fderiv, xi)
	}
	return math.NaN()
}

// Iteration function for the secant method
func secantIteration(f SingleVarFunction, xi float64, xi_1 float64) float64 {
	fi := f(xi)
	fi_1 := f(xi_1)
	if fi == fi_1 {
		return math.NaN()
	}
	return xi - fi*(xi-xi_1)/(fi-fi_1)
}

// NSimpleSolveSecant attempts to find a root of the function f starting
// at x0 using the Secant method. The function chooses a second starting point.
// A `root` value of NaN means the function failed.
func NSimpleSolveSecant(f SingleVarFunction, x0 float64, maxIterations int,
	epsilon float64) (root float64) {
	var (
		xi  float64 = x0
		xi1 float64 = x0 + hsolve
		fi  float64
	)
	for i := 0; i < maxIterations || maxIterations == 0; i++ {
		if math.IsNaN(xi) {
			return xi
		}
		fi = f(xi)
		if math.Abs(fi) < epsilon {
			return xi
		}
		xi = secantIteration(f, xi1, xi)
	}
	return math.NaN()
}

// TODO: Implement Generic solve
// Use the iteration functions for each solver
