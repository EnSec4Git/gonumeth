package gonumeth

import (
	//"fmt"
	"github.com/skelterjohn/go.matrix"
	"math"
)

const (
	hf float64 = 0.001
)

// MultiVarFunction is a type used to indicate a function that takes multiple
// arguments and (in general) returns a vector. Vectors are represented as
// matrices, so the return type is Matrix.
type MultiVarFunction func(matrix.Matrix) matrix.Matrix

// Common type that encompasses all solvers for systems of equations
type NSolveSystemMethod func(f MultiVarFunction, x0 matrix.Matrix,
	maxIterations int, epsilon float64)

// Check if a matrix contains NaN or Inf
func matrixIsInvalid(m matrix.Matrix) bool {
	for i := 0; i < m.Rows(); i++ {
		for j := 0; j < m.Cols(); j++ {
			val := m.Get(i, j)
			if math.IsNaN(val) || math.IsInf(val, 0) {
				return true
			}
		}
	}
	return false
}

// Determines whether a matrix is zero with precision epsilon
func matrixIsZero(m matrix.Matrix, epsilon float64) bool {
	for i := 0; i < m.Rows(); i++ {
		for j := 0; j < m.Cols(); j++ {
			if math.Abs(m.Get(i, j)) > epsilon {
				return false
			}
		}
	}
	return true
}

// NSolveSystemFixedPoint tries to solve a nonlinear system of equations
// by the fixed-point iteration method. The result is returned as `root`.
// A value of nil indicates the function failed to find a root.
// The parameter x0 indicates the starting point of the iteration process
// as a vector.
func NSolveSystemFixedPoint(f MultiVarFunction, x0 matrix.Matrix,
	maxIterations int, epsilon float64) (root matrix.Matrix) {
	var (
		xi matrix.Matrix = x0
		fi matrix.Matrix
	)
	for i := 0; i < maxIterations || maxIterations == 0; i++ {
		if matrixIsInvalid(xi) {
			return nil
		}
		fi = f(xi)
		if matrixIsZero(fi, epsilon) {
			return xi
		}
		xi = matrix.Sum(xi, matrix.Scaled(fi, -1.0))
	}
	return nil
}

// Calculates the Jacobian matrix of f at x0.
func jacobianOfSystem(f MultiVarFunction,
	x0 matrix.Matrix) (result matrix.Matrix) {
	n := x0.Cols()
	result = matrix.Zeros(n, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			delta := matrix.Zeros(x0.Rows(), x0.Cols())
			delta.Set(0, i, 1)
			lambda := func(x float64) float64 {
				var interimResult matrix.Matrix = f(matrix.Sum(x0,
					matrix.Scaled(delta, x)))
				return interimResult.Get(0, j)
			}
			deriv := NDifferentiateCentral(lambda, 0, hf)
			result.Set(i, j, deriv)
		}
	}
	return
}

// NSolveSystemNewton tries to solve a nonlinear system of equations
// using the Newton's method for systems. The parameter x0 indicates the
// starting point of the iteration process as vector.
// The result is returned as `root`. A value of nil indicates the function
// failed.
func NSolveSystemNewton(f MultiVarFunction, x0 matrix.Matrix, maxIterations int,
	epsilon float64) (root matrix.Matrix) {
	var (
		xi    matrix.Matrix = x0
		fi    matrix.Matrix
		Ji    matrix.Matrix
		delta matrix.Matrix
	)
	for i := 0; i < maxIterations || maxIterations == 0; i++ {
		if matrixIsInvalid(xi) {
			return nil
		}
		fi = f(xi)
		if matrixIsZero(fi, epsilon) {
			return xi
		}
		Ji = jacobianOfSystem(f, xi)
		if Ji.Det() == 0 {
			return nil
		}
		inv := matrix.Inverse(Ji)
		delta = matrix.Scaled(matrix.Product(fi, inv), -1.0)
		xi = matrix.Sum(xi, delta)
	}
	return nil
}

// Calculates n derivatives (df_i / dx_i) of f at x0
func simpleDerivs(f MultiVarFunction, x0 matrix.Matrix) (result matrix.Matrix) {
	result = matrix.MakeDenseCopy(x0)
	for i := 0; i < result.Cols(); i++ {
		delta := matrix.Zeros(result.Rows(), result.Cols())
		delta.Set(0, i, 1)
		lambda := func(x float64) float64 {
			var interimResult matrix.Matrix = f(matrix.Sum(x0,
				matrix.Scaled(delta, x)))
			return interimResult.Get(0, i)
		}
		deriv := NDifferentiateCentral(lambda, 0, hf)
		result.Set(0, i, deriv)
	}
	return
}

// NSolveSystemDeriv tries to solve a nonlinear system of equations
// by calculating n derivatives at each step.
// The parameter x0 indicates the starting point of the iteration
// process as vector.
// The result is returned as `root`. A value of NaN indicates the function
// failed.
func NSolveSystemDeriv(f MultiVarFunction, x0 matrix.Matrix, maxIterations int,
	epsilon float64) (root matrix.Matrix) {
	var (
		xi     matrix.Matrix = x0
		fi     matrix.Matrix
		derivs matrix.Matrix
	)
	for i := 0; i < maxIterations || maxIterations == 0; i++ {
		if matrixIsInvalid(xi) {
			return nil
		}
		fi = f(xi)
		if matrixIsZero(fi, epsilon) {
			return xi
		}
		derivs = simpleDerivs(f, xi)
		for j := 0; j < xi.Rows(); j++ {
			for k := 0; k < xi.Cols(); k++ {
				newVal := xi.Get(j, k) - fi.Get(j, k)/derivs.Get(j, k)
				xi.Set(j, k, newVal)
			}
		}
	}
	return nil
}

// TODO: Implement NSolveSystemSeidel method
