// numsimplesolve_test.go
package gonumeth

import (
	//"fmt"
	"math"
	"reflect"
	"runtime"
	"testing"
)

// Useful to print which method failed on which function
func getFunctionName(i interface{}) string {
	return runtime.FuncForPC(reflect.ValueOf(i).Pointer()).Name()
}

// Sample values - useful for testing only
const (
	testepsilon    float64 = 0.001
	testiterations int     = 1000
)

func sqr(x float64) float64 {
	return x * x
}

func sqrp1(x float64) float64 {
	return x*x + 1
}

func sqsqrp1(x float64) float64 {
	return sqr(sqr(x)) + 1
}

// Functions and starting points to try and find a root for
var testFunctions = []struct {
	f  SingleVarFunction
	x0 float64
}{
	{math.Sin, 3.00},
	{math.Sin, 6.10},
	{math.Cos, 1.5},
	{math.Cos, 0},
	{math.Sin, -0.05},
}

// Functions that don't have a root - test graceful failure
var testFailFunctions = []struct {
	f  SingleVarFunction
	x0 float64
}{
	{sqrp1, 1.0},
	{sqsqrp1, 1.0},
}

// All solver methods provided
var solvers = []NSimpleSolver{
	NSimpleSolveBisection,
	NSimpleSolveHalley,
	NSimpleSolveNewton,
	NSimpleSolveSecant,
}

// Tests all solvers with the testing functions
func TestSolversTable(t *testing.T) {
	for _, tt := range testFunctions {
		for _, solver := range solvers {
			result := solver(tt.f, tt.x0, testiterations, testepsilon)
			if math.IsNaN(result) {
				t.Log("Method ", getFunctionName(solver),
					"failed to find a solution for function ",
					getFunctionName(tt.f))
			}
			if math.Abs(tt.f(result)) > testepsilon {
				t.Error("Method ", getFunctionName(solver),
					"produced wrong root ", result, " for function ",
					getFunctionName(tt.f))
			}
		}
	}
}

// Tests the correct behavior when the function has no root
func TestGracefulFailureTable(t *testing.T) {
	for _, tt := range testFailFunctions {
		for _, solver := range solvers {
			result := solver(tt.f, tt.x0, testiterations, testepsilon)
			if !math.IsNaN(result) {
				t.Error("Method ", getFunctionName(solver),
					"produced a result for a positive function ",
					getFunctionName(tt.f))
			}
		}
	}
}
