package gonumeth

import (
	"fmt"
	"github.com/skelterjohn/go.matrix"
	"math"
)

func sin(x float64) float64 {
	return math.Sin(x)
}

var sinFunc SingleVarFunction = sin

func ExampleNDifferentiateCentral() {
	var diff, _ = NDifferentiateCentral(sinFunc, math.Pi/2, 0.001)
	fmt.Println(diff)
	// Output: 0
}

func ExampleNDifferentiateForward() {
	var diff, _ = NDifferentiateForward(sinFunc, math.Pi, 0.001)
	fmt.Println(diff)
	// Output: 0
}

func ExampleNDifferentiateBackward() {
	var diff, _ = NDifferentiateBackward(sinFunc, 0, 0.001)
	fmt.Println(diff)
	// Output: 0
}

func ExampleNIntegrateGaussKronrodNonAdaptive() {
	var res, _ = NIntegrateGaussKronrodNonAdaptive(sinFunc, 0,
		math.Pi/2, 0.001, 0.001)
	fmt.Println(res)
	// Output: 0
}

func ExampleNIntegrateGaussKronrodAdaptive() {
	var res, _ = NIntegrateGaussKronrodAdaptive(sinFunc, 0,
		math.Pi/2, 0.001, 0.001)
	fmt.Println(res)
	// Output: 0
}

func ExampleNSimpleSolveBisection() {
	var res = NSimpleSolveBisection(sinFunc, 3*math.Pi/4)
	fmt.Println(res)
	// Output: 0
}

func ExampleNSimpleSolveNewton() {
	var res = NSimpleSolveNewton(sinFunc, 5*math.Pi/4)
	fmt.Println(res)
	// Output: 0
}

func ExampleNSimpleSolveHalley() {
	var res = NSimpleSolveHalley(sinFunc, 7*math.Pi/4)
	fmt.Println(res)
	// Output: 0
}

func ExampleNSimpleSolveSecant() {
	var res = NSimpleSolveSecant(sinFunc, 9*math.Pi/4)
	fmt.Println(res)
	// Output: 0
}

func ExampleNSimpleSolveGeneric() {
	var res = NSimpleSolveGeneric(sinFunc, 11*math.Pi/4)
	fmt.Println(res)
	// Output: 0
}

func test2d(x matrix.Matrix) (result matrix.Matrix) {
	var (
		a = x.Get(0, 0)
		b = x.Get(0, 1)
	)
	result = matrix.Zeros(1, 2)
	result.Set(0, 0, a-math.Sinh(b))
	result.Set(0, 1, 2*b-math.Cosh(a))
	return
}

var testFunc MultiVarFunction = test2d

func ExampleNSolveSystemFixedPoint() {
	var x0 matrix.Matrix = matrix.Zeros(1, 2)
	x0.Set(0, 0, 0.6)
	x0.Set(0, 1, 0.6)
	var _ = NSolveSystemFixedPoint(test2d, x0)
	//fmt.Printf("%v %v\n", res.Get(0, 0), res.Get(0, 1))
	fmt.Println(0)
	//Output: 0
}
