package gonumeth

import "math"

// Nodes for the non-adaptive Gauss-Kronrod method
const (
	gaussNodeCount   int = 7
	kronrodNodeCount int = 15
)

var gaussNodes = [gaussNodeCount][2]float64{
	{0.0000000000000000000000000, 0.4179591836734693877551020},
	{0.4058451513773971669066064, 0.3818300505051189449503698},
	{-0.4058451513773971669066064, 0.3818300505051189449503698},
	{0.7415311855993944398638648, 0.2797053914892766679014678},
	{-0.7415311855993944398638648, 0.2797053914892766679014678},
	{0.9491079123427585245261897, 0.1294849661688696932706114},
	{-0.9491079123427585245261897, 0.1294849661688696932706114}}

var kronrodNodes = [kronrodNodeCount][2]float64{
	{0.0000000000000000000000000, 0.2094821410847278280129992},
	{0.2077849550078984676006894, 0.2044329400752988924141620},
	{-0.2077849550078984676006894, 0.2044329400752988924141620},
	{0.4058451513773971669066064, 0.1903505780647854099132564},
	{-0.4058451513773971669066064, 0.1903505780647854099132564},
	{0.5860872354676911302941448, 0.1690047266392679028265834},
	{-0.5860872354676911302941448, 0.1690047266392679028265834},
	{0.7415311855993944398638648, 0.1406532597155259187451896},
	{-0.7415311855993944398638648, 0.1406532597155259187451896},
	{0.8648644233597690727897128, 0.1047900103222501838398763},
	{-0.8648644233597690727897128, 0.1047900103222501838398763},
	{0.9491079123427585245261897, 0.0630920926299785532907007},
	{-0.9491079123427585245261897, 0.0630920926299785532907007},
	{0.9914553711208126392068547, 0.0229353220105292249637320},
	{-0.9914553711208126392068547, 0.0229353220105292249637320}}

// NIntegrateGaussKronrodNonAdaptive attempts to find the numeric value of the
// integral of f in the interval [a, b] using the Gauss-Kronrod rules.
// This method does not adapt to the "stiffness" of the function.
// If the error is greater than any of the targets (goalErrorAbs or
// goalErrorRel * result), the function has failed.
// The method will reuse all values of the function calculated.
// The function returns the result and error estimation as result and
// errorEstimate respectively.
func NIntegrateGaussKronrodNonAdaptive(f SingleVarFunction, a float64,
	b float64, goalErrorAbs float64, goalErrorRel float64) (result float64,
	errorEstimate float64) {
	var (
		gaussApprox   float64 = 0
		kronrodApprox float64 = 0
		xpoint        float64
		fvalue        float64
	)
	for i := 0; i < gaussNodeCount; i++ {
		xpoint = gaussNodes[i][0]*(b-a)/2 + (a+b)/2
		fvalue = f(xpoint)
		gaussApprox += fvalue * gaussNodes[i][1]
	}
	for i := 0; i < kronrodNodeCount; i++ {
		xpoint = kronrodNodes[i][0]*(b-a)/2 + (a+b)/2
		fvalue = f(xpoint)
		kronrodApprox += fvalue * kronrodNodes[i][1]
	}
	errorEstimate = (b - a) * 100 * math.Abs(gaussApprox-kronrodApprox)
	result = (b - a) * 0.25 * (kronrodApprox + gaussApprox)
	if errorEstimate > goalErrorAbs || errorEstimate > result*goalErrorRel {
		result = math.NaN()
	}
	return
}

// Recursive function for calculating the integral by adaptive Simspon's method
func simpsonAdaptiveRec(f SingleVarFunction, a float64, b float64,
	goalErrorAbs float64, goalErrorRel float64, S float64, fa float64,
	fb float64, fc float64) (result float64, errorEstimate float64) {
	c := (a + b) / 2
	h := b - a
	d := (a + c) / 2
	e := (c + b) / 2
	fd := f(d)
	fe := f(e)
	S_left := (h / 12) * (fa + 4*fd + fc)
	S_right := (h / 12) * (fc + 4*fe + fb)
	S2 := S_left + S_right
	errorEstimate = (S2 - S) / 15
	err := math.Abs(errorEstimate)
	if err <= goalErrorAbs && err <= S*goalErrorRel && err <= S2*goalErrorRel {
		result = S2 + errorEstimate
		errorEstimate = err
		return
	}
	res1, err1 := simpsonAdaptiveRec(f, a, c, goalErrorAbs/2, goalErrorRel,
		S_left, fa, fc, fd)
	res2, err2 := simpsonAdaptiveRec(f, c, b, goalErrorAbs/2, goalErrorRel,
		S_right, fc, fb, fe)
	result = res1 + res2
	errorEstimate = err1 + err2
	return
}

// NIntegrateGaussKronrodAdaptive attempts to find the numeric value of the
// integral of f in the interval [a, b] using the Simspon rule by
// further "sectioning" subintervals until desired error goals are achieved.
// The function returns the result and error estimation as result and
// errorEstimate respectively.
func NIntegrateSimpsonAdaptive(f SingleVarFunction, a float64, b float64,
	goalErrorAbs float64, goalErrorRel float64) (result float64,
	errorEstimate float64) {
	c := (a + b) / 2
	h := b - a
	var (
		fa     float64 = f(a)
		fb     float64 = f(b)
		fc     float64 = f(c)
		S_init float64 = h / 6 * (fa + 4*fc + fb)
	)
	return simpsonAdaptiveRec(f, a, b, goalErrorAbs, goalErrorRel,
		S_init, fa, fb, fc)
}

// TODO: Implement adaptive Gauss-Kronrod
