package gonumeth

// NSimpleSolveBisection attempts to find a root of the function f starting
// at x0 by using the Bisection method.
// A `root` value of NaN means the function failed.
func NSimpleSolveBisection(f SingleVarFunction, x0 float64) (result float64) {
	return
}

// NSimpleSolveNewton attempts to find a root of the function f starting
// at x0 using the Newton method. This method generally requires a "good"
// behavior of f locally (it being differentiable and having non-zero deriv.).
// If this is not true, the function might fail to produce the proper result.
// A `root` value of NaN means the function failed.
func NSimpleSolveNewton(f SingleVarFunction, x0 float64) (result float64) {
	return
}

// NSimpleSolveHalley attempts to find a root of the function f starting
// at x0 using the Halley method. This method is the (potentially) fastest
// method provided in this package, however it requires "very good" behavior
// of f locally (in particular it being differentiable and a non-zero value
// of the function and its first derivative).
// A `root` value of NaN means the function failed.
func NSimpleSolveHalley(f SingleVarFunction, x0 float64) (root float64) {
	return
}

// NSimpleSolveSecant attempts to find a root of the function f starting
// at x0 using the Secant method.
// A `root` value of NaN means the function failed.
func NSimpleSolveSecant(f SingleVarFunction, x0 float64) (root float64) {
	return
}

// NSimpleSolveGeneric is the "generic" root-finding algorithm of the package.
// It will attempt to use the fastest available method that is usable for the
// particular function (that means ill-behaving functions might use bisector
// while "good" functions use the Halley method) and will try to find a root.
// Note that this function does not make any guarantees on the particular root
// it finds. A `root` value of NaN means the function failed.
func NSimpleSolveGeneric(f SingleVarFunction, x0 float64) (root float64) {
	return
}
