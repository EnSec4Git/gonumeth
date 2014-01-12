// Package gonumeth provides basic primitives for performing
// Numerical computations. In particular, routines for numerical
// differentiation, integration, equation and system solving are provided.
// Type is assumed to be float64.
package gonumeth

// SingleVarFunction is a type used to represent a function that takes a single
// parameter. It should take a float64 argument and return a float64 result.
// This primitive is used throughout the package.
type SingleVarFunction func(float64) float64

// NDifferentiateCentral takes a single variable function (f), a value (x) and
// a step (h) and returns the numerical value of the derivative of f at x
// along with its estimated error (result and errorEstimate respectively).
// Note that an appropriate value of h is essential. This function uses a
// 3-point rule with x at the center.
func NDifferentiateCentral(f SingleVarFunction, x float64, h float64) (result float64, errorEstimate float64) {
	return
}

// NDifferentiateForward works the same way as NDifferentiateCentral3, except
// that it uses points x, x+h, x+2*h and is suitable for function undefined
// for y<x.
func NDifferentiateForward(f SingleVarFunction, x float64, h float64) (result float64, errorEstimate float64) {
	return
}

// NDifferentiateForward works the same way as NDifferentiateCentral3, except
// that it uses points x, x-h, x-2*h and is suitable for function undefined
// for y>x.
func NDifferentiateBackward(f SingleVarFunction, x float64, h float64) (result float64, errorEstimate float64) {
	return
}
