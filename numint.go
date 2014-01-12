package gonumeth

// NIntegrateGaussKronrodNonAdaptive attempts to find the numeric value of the
// integral of f in the interval [a, b] using the Gauss-Kronrod rules for a
// (almost geometric) progression of the (algebraic) degree of the method.
// This method does not adapt to the "stiffness" of the function.
// If the error is greater than any of the targets (goalErrorAbs or
// goalErrorRel * result), the function has failed.
// The method will reuse all values of the function calculated.
// The function returns the result and error estimation as result and
// errorEstimate respectively.
func NIntegrateGaussKronrodNonAdaptive(f SingleVarFunction, a float64, b float64, goalErrorAbs float64, goalErrorRel float64) (result float64, errorEstimate float64) {
	return
}

// NIntegrateGaussKronrodAdaptive attempts to find the numeric value of the
// integral of f in the interval [a, b] using the Gauss-Kronrod rules by
// further "sectioning" subintervals that have highest absolute/relative error.
// The method adapts and will split intervals in which the function is stiff
// in more subintervals.
// If the error is greater than any of the targets (goalErrorAbs or
// goalErrorRel * result), the function has failed.
// The method reuses all values of the function calculated.
// The function returns the result and error estimation as result and
// errorEstimate respectively.
func NIntegrateGaussKronrodAdaptive(f SingleVarFunction, a float64, b float64, goalErroAbs float64, goalErroRel float64) (result float64, errorEstimate float64) {
	return
}
