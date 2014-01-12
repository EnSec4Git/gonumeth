package gonumeth

import "github.com/skelterjohn/go.matrix"

// MultiVarFunction is a type used to indicate a function that takes multiple
// arguments and (in general) returns a vector. Vectors are represented as
// matrices, so the return type is Matrix.
type MultiVarFunction func(matrix.Matrix) matrix.Matrix

// NSolveSystemFixedPoint tries to solve a nonlinear system of equations
// by the fixed-point iteration method. The result is returned as `root`.
// A value of NaN indicates the function failed to find a root.
// The parameter x0 indicates the starting point of the iteration process
// as a vector.
func NSolveSystemFixedPoint(f MultiVarFunction, x0 matrix.Matrix) (root matrix.Matrix) {
	return
}

// NSolveSystemNewton tries to solve a nonlinear system of equations
// using the Newton's method for systems. The parameter x0 indicates the
// starting point of the iteration process as vector.
// The result is returned as `root`. A value of NaN indicates the function
// failed.
func NSolveSystemNewton(f MultiVarFunction, x0 matrix.Matrix) (root matrix.Matrix) {
	return
}

// NSolveSystemSeidel tries to solve a nonlinear system of equations
// using the Gauss-Seidel method (or rather its modification for nonlinear
// systems). The parameter x0 indicates the starting point of the iteration
// process as vector.
// The result is returned as `root`. A value of NaN indicates the function
// failed.
func NSolveSystemSeidel(f MultiVarFunction, x0 matrix.Matrix) (root matrix.Matrix) {
	return
}
