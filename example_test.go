package gonumeth

import (
	"fmt"
	"github.com/skelterjohn/go.matrix"
	"math"
)

// TODO: Check for results on other archs for other Example tests

const (
	maxIterations = 100
	defEpsilon    = 0.0001
	derivH        = 0.001
)

func sin(x float64) float64 {
	return math.Sin(x)
}

var sinFunc SingleVarFunction = sin

func ExampleNDifferentiateCentral() {
	var diff = NDifferentiateCentral(sinFunc, math.Pi/2, derivH)
	fmt.Println(diff)
	// Output: -2.7755575615628914e-14
}

func ExampleNDifferentiateForward() {
	var diff = NDifferentiateForward(sinFunc, math.Pi, derivH)
	fmt.Println(diff)
	// Output: -1.0000000000003382
}

func ExampleNDifferentiateBackward() {
	var diff = NDifferentiateBackward(sinFunc, 0, derivH)
	fmt.Println(diff)
	// Output: 1.0000000000003
}

func ExampleNIntegrateGaussKronrodNonAdaptive() {
	var res, _ = NIntegrateGaussKronrodNonAdaptive(sinFunc, 0,
		math.Pi, 0.001, 0.001)
	fmt.Println(res)
	// Output: 2.0000000000008953
}

func ExampleNIntegrateSimpsonAdaptive() {
	var res, _ = NIntegrateSimpsonAdaptive(sinFunc, 0,
		math.Pi/2, 0.001, 0.001)
	fmt.Println(res)
	// Output: 0.9999915654729925
}

func ExampleNSimpleSolveBisection() {
	var res = NSimpleSolveBisection(sinFunc, 3*math.Pi/4, maxIterations, defEpsilon)
	fmt.Println(res)
	// Output: 3.1415069901923447
}

func ExampleNSimpleSolveNewton() {
	var res = NSimpleSolveNewton(sinFunc, 5*math.Pi/4, maxIterations, defEpsilon)
	fmt.Println(res)
	// Output: 3.1415926409864197
}

func ExampleNSimpleSolveHalley() {
	var res = NSimpleSolveHalley(sinFunc, 7*math.Pi/4, maxIterations, defEpsilon)
	fmt.Println(res)
	// Output: 6.283185307175954
}

func ExampleNSimpleSolveSecant() {
	var res = NSimpleSolveSecant(sinFunc, 9*math.Pi/4, maxIterations, defEpsilon)
	fmt.Println(res)
	// Output: 6.283158738184373
}

func test2d(x matrix.Matrix) (result matrix.Matrix) {
	var (
		a = x.Get(0, 0)
		b = x.Get(0, 1)
	)
	result = matrix.Zeros(1, 2)
	result.Set(0, 0, a-math.Sinh(b))
	result.Set(0, 1, b-math.Cosh(a)/2)
	return
}

var testFunc MultiVarFunction = test2d

func ExampleNSolveSystemFixedPoint() {
	var x0 matrix.Matrix = matrix.Zeros(1, 2)
	x0.Set(0, 0, 0.6)
	x0.Set(0, 1, 0.6)
	var res = NSolveSystemFixedPoint(test2d, x0, maxIterations, defEpsilon)
	fmt.Printf("%v %v\n", res.Get(0, 0), res.Get(0, 1))
	//Output: 0.6462381444968512 0.6080337065319752
}
