// Package gonumeth provides basic primitives for performing
// Numerical computations. In particular, routines for numerical
// differentiation, integration, equation and system solving are provided.
// Type is assumed to be float64.
package gonumeth

const (
	d_c3_c_2 float64 = 1.0 / 12
	d_c3_c_1 float64 = -2.0 / 3
	d_c3_c1  float64 = 2.0 / 3
	d_c3_c2  float64 = -1.0 / 12
	d_f3_c0  float64 = -11.0 / 6
	d_f3_c1  float64 = 3.0
	d_f3_c2  float64 = -3.0 / 2
	d_f3_c3  float64 = 1.0 / 3
	d_b3_c0  float64 = 11.0 / 6
	d_b3_c_1 float64 = -3.0
	d_b3_c_2 float64 = 3.0 / 2
	d_b3_c_3 float64 = -1.0 / 3
)

// SingleVarFunction is a type used to represent a function that takes a single
// parameter. It should take a float64 argument and return a float64 result.
// This primitive is used throughout the package.
type SingleVarFunction func(float64) float64

// NDifferentiateCentral takes a single variable function (f), a value (x) and
// a step (h) and returns the numerical value of the derivative of f at x
// along with its estimated error (result and errorEstimate respectively).
// Note that an appropriate value of h is essential. This function uses a
// 4-point rule with x at the center.
func NDifferentiateCentral(f SingleVarFunction, x float64,
	h float64) (result float64) {
	x_2, x_1, x1, x2 := x-2*h, x-h, x+h, x+2*h
	f_2, f_1, f1, f2 := f(x_2), f(x_1), f(x1), f(x2)
	result = 1 / h * (d_c3_c_2*f_2 + d_c3_c_1*f_1 + d_c3_c1*f1 + d_c3_c2*f2)
	return
}

// NDifferentiateForward works the same way as NDifferentiateCentral, except
// that it uses points x, x+h, x+2*h, x+3*h and is suitable for function
// undefined for y<x.
func NDifferentiateForward(f SingleVarFunction, x float64,
	h float64) (result float64) {
	x0, x1, x2, x3 := x, x+h, x+2*h, x+3*h
	f0, f1, f2, f3 := f(x0), f(x1), f(x2), f(x3)
	result = 1 / h * (d_f3_c0*f0 + d_f3_c1*f1 + d_f3_c2*f2 + d_f3_c3*f3)
	return
}

// NDifferentiateForward works the same way as NDifferentiateCentral3, except
// that it uses points x, x-h, x-2*h, x-3*h and is suitable for function
// undefined for y>x.
func NDifferentiateBackward(f SingleVarFunction, x float64,
	h float64) (result float64) {
	x0, x_1, x_2, x_3 := x, x-h, x-2*h, x-3*h
	f0, f_1, f_2, f_3 := f(x0), f(x_1), f(x_2), f(x_3)
	result = 1 / h * (d_b3_c0*f0 + d_b3_c_1*f_1 + d_b3_c_2*f_2 + d_b3_c_3*f_3)
	return
}

// TODO: Implement 3-point rules
// TODO: Maybe implement n-point rules

// NDerivative returns a function that approximates the original function's
// derivative. The derivative is calculculated locally.
// Essentially this is a convenience wrapper for NDifferentiateCentral
func NDerivative(f SingleVarFunction, h float64) (f_prime SingleVarFunction) {
	f_prime = func(x float64) float64 {
		return NDifferentiateCentral(f, x, h)
	}
	return
}

// NDerivativeHigher returns the order-th derivative of f.
// It is used just like NDerivative
func NDerivativeHigher(f SingleVarFunction, order int,
	h float64) (f_ith SingleVarFunction) {
	if order == 0 {
		return f
	}
	f_ith = func(x float64) float64 {
		return NDifferentiateCentral(NDerivativeHigher(f, order-1, h), x, h)
	}
	return
}
