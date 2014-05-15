package gonumeth

import (
	"math"
	//"reflect"
	//"strconv"
	//"fmt"
)

const hsolve float64 = 0.01

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
func bisectionIteration(f CachedSingleVarFunction, xi float64,
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

func bisectionFindInterval(f CachedSingleVarFunction, x0 float64,
	maxIterations int) (left float64, right float64, usedSteps int) {
	var (
		step                float64 = hsolve
		xi                  float64 = x0 + step
		fi                  float64 = f(x0)
		i                   int     = 0
		pseudoMaxIterations int
	)
	if maxIterations != 0 {
		pseudoMaxIterations = maxIterations
	} else {
		pseudoMaxIterations = 100
	}
	for ; i < int(float64(pseudoMaxIterations)*0.25); i++ {
		if fi*f(xi) <= 0 {
			return x0, xi, i
		}
		step += hsolve
		xi = x0 + step
	}
	step = hsolve
	for ; i < int(float64(pseudoMaxIterations)*0.35) || maxIterations == 0; i++ {
		if fi*f(xi) <= 0 {
			return x0, xi, i
		}
		step *= 2
		xi = x0 + step
	}
	return math.NaN(), math.NaN(), pseudoMaxIterations
}

// NSimpleSolveBisection attempts to find a root of the function f starting
// at x0 by using the Bisection method.
// A `root` value of NaN means the function failed.
func NSimpleSolveBisection(f0 SingleVarFunction, x0 float64, maxIterations int,
	epsilon float64) (result float64) {
	f := CacheFunction(f0)
	xi_1, xi, i := bisectionFindInterval(f, x0, maxIterations)
	if math.IsNaN(xi) || math.IsNaN(xi_1) {
		return math.NaN()
	}
	var fi float64
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
func newtonIteration(f CachedSingleVarFunction, xi float64) float64 {
	deriv := NDifferentiateCentral(SingleVarFunction(f), xi, hsolve)
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
func NSimpleSolveNewton(f0 SingleVarFunction, x0 float64, maxIterations int,
	epsilon float64) (result float64) {
	var (
		xi float64 = x0
		fi float64
		f  CachedSingleVarFunction = CacheFunction(f0)
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
func halleyIteration(f CachedSingleVarFunction, fprime CachedSingleVarFunction,
	xi float64) float64 {
	fi := f(xi)
	f_di := fprime(xi)
	f_d2i := NDifferentiateCentral(SingleVarFunction(fprime), xi, hsolve)
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
func NSimpleSolveHalley(f0 SingleVarFunction, x0 float64, maxIterations int,
	epsilon float64) (root float64) {
	var (
		f  CachedSingleVarFunction = CacheFunction(f0)
		xi float64                 = x0
		fi float64
	)
	fderiv := CacheFunction(NDerivative(SingleVarFunction(f), hsolve))
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
func secantIteration(f CachedSingleVarFunction, xi float64,
	xi_1 float64) float64 {
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
func NSimpleSolveSecant(f1 SingleVarFunction, x0 float64, maxIterations int,
	epsilon float64) (root float64) {
	var (
		f   CachedSingleVarFunction = CacheFunction(f1)
		xi  float64                 = x0
		xi1 float64                 = x0 + hsolve
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

// A dedicated type that indicates the how greedy the algorithm for iteration
// should be.
type solverGreed int16

const (
	greedLowest  solverGreed = iota
	greedMedium              = iota
	greedHigh                = iota
	greedHighest             = iota
)

// Attempts to "raise" the "greediness". Returns the same if already max.
func greedier(level solverGreed) solverGreed {
	switch level {
	case greedLowest:
		return greedMedium
	case greedMedium:
		return greedHigh
	default:
		return greedHighest
	}
}

// Inner function to be used for dispatching the correct iteration function.
// Used when the bisection method is "fair game", that is the function should
// be pretty "bad" in behavior.
func genericSolveInnerAliasWithBisection(f CachedSingleVarFunction,
	fprime CachedSingleVarFunction, x0 float64,
	x1 float64, level solverGreed) (root float64) {
	switch level {
	case greedMedium:
		xmid := 0.5 * (x0 + x1)
		return newtonIteration(f, xmid)
	case greedHigh:
		xmid := 0.5 * (x0 + x1)
		return secantIteration(f, xmid, xmid+hsolve)
	case greedHighest:
		xmid := 0.5 * (x0 + x1)
		return halleyIteration(f, fprime, xmid)
	case greedLowest:
		l, r := bisectionIteration(f, x0, x1)
		if l != x0 {
			return l
		} else {
			return r
		}
	default:
		panic("Wrong argument at genericSolveInnerAliasWithBisection")
	}
}

// Dispatcher for the correct iteration function, depending on the level of
// greed, exerted by the algorithm.
func genericSolveInnerAlias(f CachedSingleVarFunction,
	fprime CachedSingleVarFunction, xi float64,
	level solverGreed) (root float64) {
	switch level {
	case greedMedium:
		return newtonIteration(f, xi)
	case greedHigh:
		return secantIteration(f, xi, xi+hsolve)
	case greedHighest:
		return halleyIteration(f, fprime, xi)
	default:
		//panic(strconv.Itoa(int(level)) + strconv.Itoa(int(greedHighest)))
		panic("Wrong argument at genericSolveInnerAlias")
	}
}

// Work in progress! Don't rely on this
// NSimpleSolveGeneric attempts to find a root of the function f starting
// at x0 using a variety of methods. This is *NOT* the fastest function to
// do so. Use at your own discretion.
// A `root` value of NaN means the function failed.
func NSimpleSolveGeneric(f1 SingleVarFunction, x0 float64, maxIterations int,
	epsilon float64) (root float64) {
	var (
		f                   CachedSingleVarFunction = CacheFunction(f1)
		fprime              CachedSingleVarFunction = CacheFunction(NDerivative(SingleVarFunction(f), hsolve))
		xi                  float64                 = x0
		xi1                 float64
		fi                  float64
		i                   int
		currentGreed        solverGreed = greedHighest
		pseudoMaxIterations int
	)
	if maxIterations == 0 {
		pseudoMaxIterations = 200
	} else {
		pseudoMaxIterations = maxIterations
	}
	// We start with the greedier methods (no bisection)
	for i = 0; i < int(math.Floor(float64(pseudoMaxIterations)*0.6)); i++ {
		//fmt.Println(int(currentGreed))
		xi1 = genericSolveInnerAlias(f, fprime, xi, currentGreed)
		if math.IsNaN(xi1) || math.IsInf(xi1, 0) {
			if currentGreed > greedMedium {
				currentGreed = currentGreed - 1
				xi1 = xi
			} else {
				break
			}
		} else {
			//(?) Check whether to increase greed with only 1, or elevate to max
			currentGreed = greedier(currentGreed)
			fi = f(xi)
			if math.Abs(fi) < epsilon {
				return xi
			}
			xi = xi1
		}
	}
	// We have to use bisection
	var iterationsLeft int
	if maxIterations == 0 {
		iterationsLeft = 0
	} else {
		iterationsLeft = maxIterations - i
	}
	if iterationsLeft == 0 {
		return math.NaN()
	}
	left, right, spent := bisectionFindInterval(f, x0, iterationsLeft)
	i += spent
	if math.IsNaN(left) || math.IsNaN(right) {
		return math.NaN()
	}
	mid := (left + right) * 0.5
	fright := f(right)
	highGreedFailures := 0
	currentGreed = greedHighest
	for ; i < maxIterations || maxIterations == 0; i++ {
		xi1 = genericSolveInnerAliasWithBisection(f, fprime, left,
			right, currentGreed)
		if math.IsNaN(xi1) || math.IsInf(xi1, 0) {
			if currentGreed != greedLowest {
				currentGreed = currentGreed - 1
			}
			highGreedFailures++
		} else {
			//(?) Check whether to increase greed with only 1, or elevate to max
			if math.Abs(fi) < epsilon {
				return mid
			}
			currentGreed = greedier(currentGreed)
			if xi1 < right && xi1 < left {
				mid = xi1
			}
			fi = f(mid)
			if fi*fright > 0 {
				right = mid
				fright = fi
			} else {
				left = mid
			}
		}
	}
	return
}
