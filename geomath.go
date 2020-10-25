package ggeodesic

import "math"

func cbrt(x float64) float64 {
	// Real cube root of a number
	y := math.Pow(math.Abs(x), 1.0/3.0)
	switch {
	case x > 0:
		return y
	case x < 0:
		return -y
	default:
		return x
	}
}

func polyval(N int, p []float64, s int, x float64) float64 {
	// Evaluate a polynomial.
	var y float64
	if N < 0 {
		y = 0
	} else {
		y = p[s]
	}
	for N > 0 {
		N -= 1
		s += 1
		y = y*x + p[s]
	}
	return y
}

func sum(u, v float64) (float64, float64) {
	// Error free transformation of a sum.
	// Note that t can be the same as one of the first two arguments.
	s := u + v
	up := s - v
	vpp := s - up
	up -= u
	vpp -= v
	t := -(up + vpp)
	// u + v =       s      + t
	//       = round(u + v) + t
	return s, t
}

func remainder(x, y float64) float64 {
	// remainder of x/y in the range [-y/2, y/2].
	z := math.NaN()
	if !math.IsInf(x, 0) {
		z = math.Mod(x, y)
	}
	switch {
	case z < -y/2:
		return z + y
	case z < y/2:
		return z
	default:
		return z - y
	}
}

func angNormalize(x float64) float64 {
	// reduce angle to (-180,180]
	y := remainder(x, 360)
	if y == -180 {
		return 180
	}
	return y
}

func angDiff(x, y float64) (float64, float64) {
	// compute y - x and reduce to [-180,180] accurately
	d, t := sum(angNormalize(-x), angNormalize(y))
	d = angNormalize(d)
	if d == 180 && t > 0 {
		return sum(-180, t)
	}
	return sum(d, t)
}

func angRound(x float64) float64 {
	// Round an angle so that small values underflow to zero."""
	// The makes the smallest gap in x = 1/16 - nextafter(1/16, 0) = 1/2^57
	// for reals = 0.7 pm on the earth if x is an angle in degrees.  (This
	// is about 1000 times more resolution than we get with angles around 90
	// degrees.)  We use this to avoid having to deal with near singular
	// cases when x is non-zero but tiny (e.g., 1.0e-200).
	z := 1 / 16.0
	y := math.Abs(x)
	// The compiler mustn't "simplify" z - (z - y) to y
	if y < z {
		y = z - (z - y)
	}
	switch {
	case x == 0:
		return 0.0
	case x < 0:
		return -y
	default:
		return y
	}
}

func radians(deg float64) float64 {
	return math.Pi * deg / 180
}

func degrees(rad float64) float64 {
	return rad * 180 / math.Pi
}

func atan2d(y, x float64) float64 {
	//"""compute atan2(y, x) with the result in degrees"""
	var q float64
	if math.Abs(y) > math.Abs(x) {
		q = 2
		x, y = y, x
	} else {
		q = 0
	}
	if x < 0 {
		q += 1
		x = -x
	}
	ang := degrees(math.Atan2(y, x))
	switch q {
	case 1:
		if y >= 0 {
			ang = 180 - ang
		} else {
			ang = -180 - ang
		}
	case 2:
		ang = 90 - ang
	case 3:
		ang = -90 + ang
	}
	return ang
}

func norm(x, y float64) (float64, float64) {
	//"""Private: Normalize a two-vector."""
	r := math.Hypot(x, y)
	return x / r, y / r
}

func latFix(x float64) float64 {
	// replace angles outside [-90,90] by NaN
	if math.Abs(x) > 90 {
		return math.NaN()
	}
	return x
}

func sincosd(x float64) (float64, float64) {
	//Compute sine and cosine of x in degrees.
	r := math.NaN()
	if !math.IsInf(x, 0) {
		r = math.Mod(x, 360)
	}
	q := 0
	if !math.IsNaN(r) {
		q = int(math.Round(r / 90))
	}
	r -= float64(90 * q)
	r = radians(r)
	s := math.Sin(r)
	c := math.Cos(r)
	q = q % 4
	switch q {
	case 1:
		s, c = c, -s
	case 2:
		s, c = -s, -c
	case 3:
		s, c = -c, s
	}
	if x == 0 {
		return x, c
	}
	return s, c
}

func makeRange(min, max float64) []float64 {
	a := make([]float64, int(max-min))
	for i := range a {
		a[i] = min + float64(i)
	}
	return a
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
