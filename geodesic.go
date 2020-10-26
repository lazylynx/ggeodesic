package ggeodesic

import (
	"math"
)

const geographicLibGeodesicOrder = 6
const nA1 = geographicLibGeodesicOrder
const nC1 = geographicLibGeodesicOrder
const nC1p = geographicLibGeodesicOrder
const nA2 = geographicLibGeodesicOrder
const nC2 = geographicLibGeodesicOrder
const nA3 = geographicLibGeodesicOrder
const nA3x = nA3
const nC3 = geographicLibGeodesicOrder
const nC3x = (nC3 * (nC3 - 1)) / 2
const nC4 = geographicLibGeodesicOrder
const nC4x = (nC4 * (nC4 + 1)) / 2
const digits = 53
const maxIt1 = 20
const maxIt2 = maxIt1 + digits + 10

type Geodesic struct {
	//Solve geodesic problems
	a, f, f1, e2, ep2, n, b, c2, etol2, tol0, tol1, tolb, xthresh, tiny float64
	a3x, c3x, c4x                                                       []float64
}

func WGS84() *Geodesic {
	// Instantiation for the WGS84 ellipsoid
	// the equatorial radius in meters of the WGS84 ellipsoid in meters
	a := 6378137.0
	// the flattening of the WGS84 ellipsoid, 1/298.257223563
	f := 1 / 298.257223563
	return New(a, f)
}

func New(a, f float64) *Geodesic {
	// Construct a Geodesic object
	// :param a: the equatorial radius of the ellipsoid in meters
	// :param f: the flattening of the ellipsoid
	//
	// An exception is thrown if *a* or the polar semi-axis *b* = *a* (1 -
	//    *f*) is not a finite positive quantity.
	epsilon := math.Pow(2.0, 1-digits)
	tol0_ := epsilon
	tol2_ := math.Sqrt(tol0_)
	e2 := f * (2 - f)
	b := a * (1 - f)

	if !(!math.IsInf(a, 0) && a > 0) {
		panic("Equatorial radius is not positive")
	}
	if !(!math.IsInf(b, 0) && b > 0) {
		panic("Polar semi-axis is not positive")
	}

	// authalic radius squared
	c2 := math.Pow(a, 2) + math.Pow(b, 2)
	switch {
	case e2 == 0:
		c2 /= 2
	case e2 > 0:
		c2 *= math.Atanh(math.Sqrt(e2))
		c2 /= math.Sqrt(math.Abs(e2))
		c2 /= 2
	default:
		c2 *= math.Atan(math.Sqrt(-e2))
		c2 /= math.Sqrt(math.Abs(e2))
		c2 /= 2
	}

	geod := Geodesic{
		a:   a,
		f:   f,
		f1:  1 - f,
		e2:  e2,
		ep2: e2 / math.Pow(1-f, 2),
		n:   f / (2 - f),
		b:   b,
		c2:  c2,
		// The sig12 threshold for "really short".  Using the auxiliary sphere
		// solution with dnm computed at (bet1 + bet2) / 2, the relative error in
		// the azimuth consistency check is sig12^2 * abs(f) * min(1, 1-f/2) / 2.
		// (Error measured for 1/100 < b/a < 100 and abs(f) >= 1/1000.  For a given
		// f and sig12, the max error occurs for lines near the pole.  If the old
		// rule for computing dnm = (dn1 + dn2)/2 is used, then the error increases
		// by a factor of 2.)  Setting this equal to epsilon gives sig12 = etol2.
		// Here 0.1 is a safety factor (error decreased by 100) and max(0.001,
		// abs(f)) stops etol2 getting too large in the nearly spherical case.
		etol2: 0.1 * tol2_ / math.Sqrt(math.Max(0.001, math.Abs(f))*
			math.Min(1.0, 1-f/2)/2),
		a3x:     makeRange(0, nA3x),
		c3x:     makeRange(0, nC3x),
		c4x:     makeRange(0, nC4x),
		tol0:    tol0_,
		tol1:    200 * tol0_,
		tolb:    tol0_ * tol2_,
		xthresh: 1000 * tol2_,
		tiny:    math.Sqrt(math.Pow(2.0, -1022)),
	}
	// init a3x c3x c4x
	geod.a3Coeff()
	geod.c3Coeff()
	geod.c4Coeff()
	return &geod
}

func (geod *Geodesic) a3Coeff() {
	// return coefficients for A3
	coeff := []float64{-3, 128, -2, -3, 64, -1, -3, -1, 16, 3, -1, -2, 8, 1, -1, 2, 1, 1}
	o := 0
	k := 0
	for j := nA3 - 1; j > -1; j-- { // coeff of eps^j
		m := min(nA3-j-1, j) // order of polynomial in n
		geod.a3x[k] = polyval(m, coeff, o, geod.n) / coeff[o+m+1]
		k += 1
		o += m + 2
	}
}

func (geod *Geodesic) c3Coeff() {
	// return coefficients for C3
	coeff := []float64{3, 128, 2, 5, 128, -1, 3, 3, 64, -1, 0, 1, 8, -1, 1, 4, 5, 256, 1, 3, 128, -3, -2, 3, 64, 1, -3, 2, 32, 7, 512, -10, 9, 384, 5, -9, 5, 192, 7, 512, -14, 7, 512, 21, 2560}
	o := 0
	k := 0
	for l := 1; l < nC3; l++ {
		for j := nC3 - 1; j > l-1; j-- {
			m := min(nC3-j-1, j) // order of polynomial in n
			geod.c3x[k] = polyval(m, coeff, o, geod.n) / coeff[o+m+1]
			k += 1
			o += m + 2
		}
	}
}

func (geod *Geodesic) c4Coeff() {
	// return coefficients for C4
	coeff := []float64{97, 15015, 1088, 156, 45045, -224, -4784, 1573, 45045, -10656, 14144, -4576, -858, 45045, 64, 624, -4576, 6864, -3003, 15015,
		100, 208, 572, 3432, -12012, 30030, 45045, 1, 9009, -2944, 468, 135135, 5792, 1040, -1287, 135135, 5952, -11648, 9152, -2574, 135135, -64, -624, 4576, -6864, 3003, 135135,
		8, 10725, 1856, -936, 225225, -8448, 4992, -1144, 225225, -1440, 4160, -4576, 1716, 225225,
		-136, 63063, 1024, -208, 105105, 3584, -3328, 1144, 315315, -128, 135135, -2560, 832, 405405, 128, 99099,
	}

	o := 0
	k := 0
	for l := 0; l < nC4; l++ { // l is index of C4[l]
		for j := nC4 - 1; j > l-1; j-- { // coeff of eps^j
			m := nC4 - j - 1 // order of polynomial in n
			geod.c4x[k] = polyval(m, coeff, o, geod.n) / coeff[o+m+1]
			k += 1
			o += m + 2
		}
	}
}

func (geod *Geodesic) Inverse(lat1, lon1, lat2, lon2 float64) map[string]float64 {
	return geod.inverse(lat1, lon1, lat2, lon2, STANDARD)
}

func (geod *Geodesic) inverse(lat1, lon1, lat2, lon2 float64, outmask int) map[string]float64 {
	//Solve the inverse geodesic problem
	//
	//:param lat1: latitude of the first point in degrees
	//:param lon1: longitude of the first point in degrees
	//:param lat2: latitude of the second point in degrees
	//:param lon2: longitude of the second point in degrees
	//:param outmask: the :ref:`output mask <outmask>`
	//:return: a :ref:`dict`
	//
	//Compute geodesic between (*lat1*, *lon1*) and (*lat2*, *lon2*).
	//The default value of *outmask* is STANDARD, i.e., the *lat1*,
	//*lon1*, *azi1*, *lat2*, *lon2*, *azi2*, *s12*, *a12* entries are
	//returned.
	a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12 := geod.genInverse(
		lat1, lon1, lat2, lon2, outmask)
	outmask &= outMask
	if outmask&longUnroll > 0 {
		lon12, e := angDiff(lon1, lon2)
		lon2 = (lon1 + lon12) + e
	} else {
		lon2 = angNormalize(lon2)
	}
	result := make(map[string]float64)
	result["lat1"] = latFix(lat1)
	if outmask&longUnroll > 0 {
		result["lon1"] = lon1
	} else {
		result["lon1"] = angNormalize(lon1)
	}
	result["lat2"] = latFix(lat2)
	result["lon2"] = lon2
	result["a12"] = a12
	if outmask&DISTANCE > 0 {
		result["s12"] = s12
	}
	if outmask&AZIMUTH > 0 {
		result["azi1"] = atan2d(salp1, calp1)
		result["azi2"] = atan2d(salp2, calp2)
	}
	if outmask&reducedLength > 0 {
		result["m12"] = m12
	}
	if outmask&geodesicScale > 0 {
		result["M12"] = M12
		result["M21"] = M21
	}
	if outmask&area > 0 {
		result["S12"] = S12
	}
	return result
}

func (geod *Geodesic) a3f(eps float64) float64 {
	// Evaluate A3
	return polyval(nA3-1, geod.a3x, 0, eps)
}

func (geod *Geodesic) c3f(eps float64, c []float64) {
	// Evaluate C3
	// Elements c[1] thru c[nC3 - 1] are set
	mult := 1.0
	o := 0
	for l := 1; l < nC3; l++ { // l is index of C3[l]
		m := nC3 - l - 1 // order of polynomial in eps
		mult *= eps
		c[l] = mult * polyval(m, geod.c3x, o, eps)
		o += m + 1
	}
}

func (geod *Geodesic) c4f(eps float64, c []float64) {
	// Evaluate C4 coeffs by Horner's method
	// Elements c[0] thru c[nC4 - 1] are set
	mult := 1.0
	o := 0
	for l := 0; l < nC4; l++ { // l is index of C4[l]
		m := nC4 - l - 1 // order of polynomial in eps
		c[l] = mult * polyval(m, geod.c4x, o, eps)
		o += m + 1
		mult *= eps
	}
}

func (geod *Geodesic) lengths(eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2 float64,
	outmask int, C1a, C2a []float64) (s12b, m12b, m0, M12, M21 float64) {
	// return a bunch of lengths
	// Return s12b, m12b, m0, M12, M21, where
	// m12b = (reduced lengths)/_b; s12b = distance/_b,
	// m0 = coefficient of secular term in expression for reduced lengths.
	outmask &= outMask
	// outmask & DISTANCE: set s12b
	// outmask & reducedLength: set m12b & m0
	// outmask & geodesicScale: set M12 & M21
	s12b = math.NaN()
	m12b = math.NaN()
	m0 = math.NaN()
	M12 = math.NaN()
	M21 = math.NaN()
	A1 := math.NaN()
	A2 := math.NaN()
	m0x := math.NaN()
	J12 := math.NaN()
	if outmask&(DISTANCE|reducedLength|geodesicScale) > 0 {
		A1 = a1m1f(eps)
		c1f(eps, C1a)
		if outmask&(reducedLength|geodesicScale) > 0 {
			A2 = a2m1f(eps)
			c2f(eps, C2a)
			m0x = A1 - A2
			A2 = 1 + A2
		}
		A1 = 1 + A1
	}
	if outmask&DISTANCE > 0 {
		B1 := sinCosSeries(true, ssig2, csig2, C1a) -
			sinCosSeries(true, ssig1, csig1, C1a)
		// Missing a factor of _b
		s12b = A1 * (sig12 + B1)
		if outmask&(reducedLength|geodesicScale) > 0 {
			B2 := sinCosSeries(true, ssig2, csig2, C2a) -
				sinCosSeries(true, ssig1, csig1, C2a)
			J12 = m0x*sig12 + (A1*B1 - A2*B2)
		}
	} else if outmask&(reducedLength|geodesicScale) > 0 {
		// Assume here that nC1 >= nC2
		for l := 1; l < nC2; l++ {
			C2a[l] = A1*C1a[l] - A2*C2a[l]
		}
		J12 = m0x*sig12 + (sinCosSeries(true, ssig2, csig2, C2a) -
			sinCosSeries(true, ssig1, csig1, C2a))
	}
	if outmask&reducedLength > 0 {
		m0 = m0x
		// Missing a factor of _b.
		// Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
		// accurate cancellation in the case of coincident points.
		m12b = dn2*(csig1*ssig2) - dn1*(ssig1*csig2) - csig1*csig2*J12
	}
	if outmask&geodesicScale > 0 {
		csig12 := csig1*csig2 + ssig1*ssig2
		t := geod.ep2 * (cbet1 - cbet2) * (cbet1 + cbet2) / (dn1 + dn2)
		M12 = csig12 + (t*ssig2-csig2*J12)*ssig1/dn1
		M21 = csig12 - (t*ssig1-csig1*J12)*ssig2/dn2
	}
	return
}

func (geod *Geodesic) inverseStart(sbet1, cbet1, dn1, sbet2, cbet2, dn2,
	lam12, slam12, clam12 float64, C1a, C2a []float64) (sig12, salp1, calp1, salp2, calp2, dnm float64) {
	// Find a starting value for Newton's method.
	// Return a starting point for Newton's method in salp1 and calp1 (function
	// value is -1).  If Newton's method doesn't need to be used, return also
	// salp2 and calp2 and function value is sig12.
	sig12 = -1
	salp2 = math.NaN()
	calp2 = math.NaN()
	dnm = math.NaN() // Return values
	// bet12 = bet2 - bet1 in [0, pi); bet12a = bet2 + bet1 in (-pi, 0]
	sbet12 := sbet2*cbet1 - cbet2*sbet1
	cbet12 := cbet2*cbet1 + sbet2*sbet1
	// Volatile declaration needed to fix inverse cases
	// 88.202499451857 0 -88.202499451857 179.981022032992859592
	// 89.262080389218 0 -89.262080389218 179.992207982775375662
	// 89.333123580033 0 -89.333123580032997687 179.99295812360148422
	// which otherwise fail with g++ 4.4.4 x86 -O3
	sbet12a := sbet2 * cbet1
	sbet12a += cbet2 * sbet1

	var somg12, comg12 float64
	shortline := cbet12 >= 0 && sbet12 < 0.5 && cbet2*lam12 < 0.5
	if shortline {
		sbetm2 := math.Pow(sbet1+sbet2, 2)
		// sin((bet1+bet2)/2)^2
		// =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2)
		sbetm2 /= sbetm2 + math.Pow(cbet1+cbet2, 2)
		dnm = math.Sqrt(1 + geod.ep2*sbetm2)
		omg12 := lam12 / (geod.f1 * dnm)
		somg12 = math.Sin(omg12)
		comg12 = math.Cos(omg12)
	} else {
		somg12 = slam12
		comg12 = clam12
	}

	salp1 = cbet2 * somg12
	if comg12 >= 0 {
		calp1 = sbet12 + cbet2*sbet1*math.Pow(somg12, 2)/(1+comg12)
	} else {
		calp1 = sbet12a - cbet2*sbet1*math.Pow(somg12, 2)/(1-comg12)
	}
	ssig12 := math.Hypot(salp1, calp1)
	csig12 := sbet1*sbet2 + cbet1*cbet2*comg12

	if shortline && ssig12 < geod.etol2 {
		// really short lines
		salp2 = cbet1 * somg12
		var mult float64
		if comg12 >= 0 {
			mult = math.Pow(somg12, 2) / (1 + comg12)
		} else {
			mult = 1 - comg12
		}
		calp2 = sbet12 - cbet1*sbet2*mult
		salp2, calp2 = norm(salp2, calp2)
		// Set return value
		sig12 = math.Atan2(ssig12, csig12)
	} else if math.Abs(geod.n) >= 0.1 || // Skip astroid calc if too eccentric
		csig12 >= 0 ||
		ssig12 >= 6*math.Abs(geod.n)*math.Pi*math.Pow(cbet1, 2) {
		// Nothing to do, zeroth order spherical approximation is OK
	} else {
		// Scale lam12 and bet2 to x, y coordinate system where antipodal point
		// is at origin and singular point is at y = 0, x = -1.
		// real y, lamscale, betscale
		// Volatile declaration needed to fix inverse case
		// 56.320923501171 0 -56.320923501171 179.664747671772880215
		// which otherwise fails with g++ 4.4.4 x86 -O3
		// volatile real x
		var x, y, lamscale float64
		lam12x := math.Atan2(-slam12, -clam12)
		if geod.f >= 0 {
			// In fact f == 0 does not get here
			// x = dlong, y = dlat
			k2 := math.Pow(sbet1, 2) * geod.ep2
			eps := k2 / (2*(1+math.Sqrt(1+k2)) + k2)
			lamscale = geod.f * cbet1 * geod.a3f(eps) * math.Pi
			betscale := lamscale * cbet1
			x = lam12x / lamscale
			y = sbet12a / betscale
		} else {
			// _f < 0
			// x = dlat, y = dlong
			cbet12a := cbet2*cbet1 - sbet2*sbet1
			bet12a := math.Atan2(sbet12a, cbet12a)
			// real m12b, m0, dummy
			// In the case of lon12 = 180, this repeats a calculation made in
			// Inverse.
			_, m12b, m0, _, _ := geod.lengths(
				geod.n, math.Pi+bet12a, sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
				cbet1, cbet2, reducedLength, C1a, C2a)
			x = -1 + m12b/(cbet1*cbet2*m0*math.Pi)
			var betscale float64
			if x < -0.01 {
				betscale = sbet12a / x
			} else {
				betscale = -geod.f * math.Pow(cbet1, 2) * math.Pi
			}
			lamscale = betscale / cbet1
			y = lam12x / lamscale
		}
		if y > -geod.tol1 && x > -1-geod.xthresh {
			// strip near cut
			if geod.f >= 0 {
				salp1 = math.Min(1.0, -x)
				calp1 = -math.Sqrt(1 - math.Pow(salp1, 2))
			} else {
				if x > -geod.tol1 {
					calp1 = math.Max(0.0, x)
				} else {
					calp1 = math.Max(-1.0, x)
				}
				salp1 = math.Sqrt(1 - math.Pow(calp1, 2))
			}
		} else {
			// Estimate alp1, by solving the astroid problem.
			//
			// Could estimate alpha1 = theta + pi/2, directly, i.e.,
			//   calp1 = y/k; salp1 = -x/(1+k);  for _f >= 0
			//   calp1 = x/(1+k); salp1 = -y/k;  for _f < 0 (need to check)
			//
			// However, it's better to estimate omg12 from astroid and use
			// spherical formula to compute alp1.  This reduces the mean number of
			// Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
			// (min 0 max 5).  The changes in the number of iterations are as
			// follows:
			//
			// change percent
			//    1       5
			//    0      78
			//   -1      16
			//   -2       0.6
			//   -3       0.04
			//   -4       0.002
			//
			// The histogram of iterations is (m = number of iterations estimating
			// alp1 directly, n = number of iterations estimating via omg12, total
			// number of trials = 148605):
			//
			//  iter    m      n
			//    0   148    186
			//    1 13046  13845
			//    2 93315 102225
			//    3 36189  32341
			//    4  5396      7
			//    5   455      1
			//    6    56      0
			//
			// Because omg12 is near pi, estimate work with omg12a = pi - omg12
			k := astroid(x, y)
			var omg12a float64
			if geod.f >= 0 {
				omg12a = lamscale * (-x * k / (1 + k))
			} else {
				omg12a = lamscale * (-y * (1 + k) / k)
			}
			somg12 = math.Sin(omg12a)
			comg12 = -math.Cos(omg12a)
			// Update spherical estimate of alp1 using omg12 instead of lam12
			salp1 = cbet2 * somg12
			calp1 = sbet12a - cbet2*sbet1*math.Pow(somg12, 2)/(1-comg12)
		}

	}
	// Sanity check on starting guess.  Backwards check allows NaN through.
	if !(salp1 <= 0) {
		salp1, calp1 = norm(salp1, calp1)
	} else {
		salp1 = 1
		calp1 = 0
	}
	return //sig12, salp1, calp1, salp2, calp2, dnm
}

func (geod *Geodesic) lambda12(sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
	slam120, clam120 float64, diffp bool, C1a, C2a, C3a []float64) (lam12, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps, domg12, dlam12 float64) {
	// Solve hybrid problem
	if sbet1 == 0 && calp1 == 0 {
		// Break degeneracy of equatorial line.  This case has already been
		// handled.
		calp1 = -geod.tiny
	}
	// sin(alp1) * cos(bet1) = sin(alp0)
	salp0 := salp1 * cbet1
	calp0 := math.Hypot(calp1, salp1*sbet1) // calp0 > 0
	// real somg1, comg1, somg2, comg2, lam12
	// tan(bet1) = tan(sig1) * cos(alp1)
	// tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
	ssig1 = sbet1
	somg1 := salp0 * sbet1
	csig1 = calp1 * cbet1
	comg1 := csig1
	ssig1, csig1 = norm(ssig1, csig1)
	// Math.norm(somg1, comg1); -- don't need to normalize!

	// Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
	// about this case, since this can yield singularities in the Newton
	// iteration.
	// sin(alp2) * cos(bet2) = sin(alp0)
	salp2 = salp1
	if cbet2 != cbet1 {
		salp2 = salp0 / cbet2
	}
	// calp2 = sqrt(1 - sq(salp2))
	//       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
	// and subst for calp0 and rearrange to give (choose positive sqrt
	// to give alp2 in [0, pi/2]).
	if cbet2 != cbet1 || math.Abs(sbet2) != -sbet1 {
		if cbet1 < -sbet1 {
			calp2 = math.Sqrt(math.Pow(calp1*cbet1, 2)+(cbet2-cbet1)*(cbet1+cbet2)) / cbet2
		} else {
			calp2 = math.Sqrt(math.Pow(calp1*cbet1, 2)+(sbet1-sbet2)*(sbet1+sbet2)) / cbet2
		}
	} else {
		calp2 = math.Abs(calp1)
	}
	// tan(bet2) = tan(sig2) * cos(alp2)
	// tan(omg2) = sin(alp0) * tan(sig2).
	ssig2 = sbet2
	somg2 := salp0 * sbet2
	csig2 = calp2 * cbet2
	comg2 := csig2
	ssig2, csig2 = norm(ssig2, csig2)
	// Math.norm(somg2, comg2); -- don't need to normalize!
	// sig12 = sig2 - sig1, limit to [0, pi]
	sig12 = math.Atan2(math.Max(0.0, csig1*ssig2-ssig1*csig2),
		csig1*csig2+ssig1*ssig2)
	// omg12 = omg2 - omg1, limit to [0, pi]
	somg12 := math.Max(0.0, comg1*somg2-somg1*comg2)
	comg12 := comg1*comg2 + somg1*somg2
	// eta = omg12 - lam120
	eta := math.Atan2(somg12*clam120-comg12*slam120,
		comg12*clam120+somg12*slam120)
	// real B312
	k2 := math.Pow(calp0, 2) * geod.ep2
	eps = k2 / (2*(1+math.Sqrt(1+k2)) + k2)
	geod.c3f(eps, C3a)
	B312 := sinCosSeries(true, ssig2, csig2, C3a) -
		sinCosSeries(true, ssig1, csig1, C3a)
	domg12 = -geod.f * geod.a3f(eps) * salp0 * (sig12 + B312)
	lam12 = eta + domg12
	if diffp {
		if calp2 == 0 {
			dlam12 = -2 * geod.f1 * dn1 / sbet1
		} else {
			_, dlam12, _, _, _ = geod.lengths(
				eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
				reducedLength, C1a, C2a)
			dlam12 *= geod.f1 / (calp2 * cbet2)
		}
	} else {
		dlam12 = math.NaN()
	}
	return
}

func (geod *Geodesic) genInverse(lat1, lon1, lat2, lon2 float64, outmask int) (a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12 float64) {
	// General version of the inverse problem
	a12 = math.NaN()
	s12 = math.NaN()
	m12 = math.NaN()
	M12 = math.NaN()
	M21 = math.NaN()
	S12 = math.NaN()

	outmask &= outMask
	// Compute longitude difference (AngDiff does this carefully).  Result is
	// in [-180, 180] but -180 is only for west-going geodesics.  180 is for
	// east-going and meridional geodesics.
	lon12, lon12s := angDiff(lon1, lon2)
	// Make longitude difference positive.
	lonsign := 1.0
	if lon12 < 0 {
		lonsign = -1
	}
	// If very close to being on the same half-meridian, then make it so.
	lon12 = lonsign * angRound(lon12)
	lon12s = angRound((180 - lon12) - lonsign*lon12s)
	lam12 := radians(lon12)
	var slam12, clam12 float64
	if lon12 > 90 {
		slam12, clam12 = sincosd(lon12s)
		clam12 = -clam12
	} else {
		slam12, clam12 = sincosd(lon12)
	}
	// If really close to the equator, treat as on equator.
	lat1 = angRound(latFix(lat1))
	lat2 = angRound(latFix(lat2))
	// Swap points so that point with higher (abs) latitude is point 1
	// If one latitude is a nan, then it becomes lat1.
	swapp := -1.0
	if math.Abs(lat1) >= math.Abs(lat2) {
		swapp = 1
	}
	if swapp < 0 {
		lonsign *= -1
		lat2, lat1 = lat1, lat2
	}
	// Make lat1 <= 0
	latsign := 1.0
	if lat1 >= 0 {
		latsign = -1
	}
	lat1 *= latsign
	lat2 *= latsign
	// Now we have
	//
	//     0 <= lon12 <= 180
	//     -90 <= lat1 <= 0
	//     lat1 <= lat2 <= -lat1
	//
	// longsign, swapp, latsign register the transformation to bring the
	// coordinates to this canonical form.  In all cases, 1 means no change was
	// made.  We make these transformations so that there are few cases to
	// check, e.g., on verifying quadrants in atan2.  In addition, this
	// enforces some symmetries in the results returned.

	// real phi, sbet1, cbet1, sbet2, cbet2, s12x, m12x
	sbet1, cbet1 := sincosd(lat1)
	sbet1 *= geod.f1
	// Ensure cbet1 = +epsilon at poles
	sbet1, cbet1 = norm(sbet1, cbet1)
	cbet1 = math.Max(geod.tiny, cbet1)
	sbet2, cbet2 := sincosd(lat2)
	sbet2 *= geod.f1
	// Ensure cbet2 = +epsilon at poles
	sbet2, cbet2 = norm(sbet2, cbet2)
	cbet2 = math.Max(geod.tiny, cbet2)
	// If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
	// |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
	// a better measure.  This logic is used in assigning calp2 in Lambda12.
	// Sometimes these quantities vanish and in that case we force bet2 = +/-
	// bet1 exactly.  An example where is is necessary is the inverse problem
	// 48.522876735459 0 -48.52287673545898293 179.599720456223079643
	// which failed with Visual Studio 10 (Release and Debug)
	if cbet1 < -sbet1 {
		if cbet2 == cbet1 {
			if sbet2 < 0 {
				sbet2 = sbet1
			} else {
				sbet2 = -sbet1
			}
		}
	} else {
		if math.Abs(sbet2) == -sbet1 {
			cbet2 = cbet1
		}
	}
	dn1 := math.Sqrt(1 + geod.ep2*math.Pow(sbet1, 2))
	dn2 := math.Sqrt(1 + geod.ep2*math.Pow(sbet2, 2))
	// real a12, sig12, calp1, salp1, calp2, salp2
	// index zero elements of these arrays are unused
	C1a := makeRange(0, nC1+1)
	C2a := makeRange(0, nC2+1)
	C3a := makeRange(0, nC3)
	var somg12, comg12, omg12, sig12, s12x, m12x float64

	meridian := lat1 == -90 || slam12 == 0
	if meridian {
		// Endpoints are on a single full meridian, so the geodesic might lie on
		// a meridian.
		calp1 = clam12
		salp1 = slam12 // Head to the target longitude
		calp2 = 1.0
		salp2 = 0.0 // At the target we're heading north
		// tan(bet) = tan(sig) * cos(alp)
		ssig1 := sbet1
		csig1 := calp1 * cbet1
		ssig2 := sbet2
		csig2 := calp2 * cbet2

		// sig12 = sig2 - sig1
		sig12 := math.Atan2(math.Max(0.0, csig1*ssig2-ssig1*csig2),
			csig1*csig2+ssig1*ssig2)

		s12x, m12x, _, M12, M21 = geod.lengths(
			geod.n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
			outmask|DISTANCE|reducedLength, C1a, C2a)
		// Add the check for sig12 since zero lengths geodesics might yield m12 <
		// 0.  Test case was
		//
		//    echo 20.001 0 20.001 0 | GeodSolve -i
		//
		// In fact, we will have sig12 > pi/2 for meridional geodesic which is
		// not a shortest path.
		if sig12 < 1 || m12x >= 0 {
			if sig12 < 3*geod.tiny {
				sig12 = 0
				m12x = 0
				s12x = 0
			}
			m12x *= geod.b
			s12x *= geod.b
			a12 = degrees(sig12)
		} else {
			// m12 < 0, i.e., prolate and too close to anti-podal
			meridian = false
		}
	}
	// somg12 > 1 marks that it needs to be calculated
	somg12 = 2.0
	comg12 = 0.0
	omg12 = 0.0
	if !meridian && sbet1 == 0 && (geod.f <= 0 || lon12s >= geod.f*180) {
		// Mimic the way Lambda12 works with calp1 = 0
		// Geodesic runs along equator
		calp1 = 0
		calp2 = 0.0
		salp1 = 1.0
		salp2 = 1.0
		s12x = geod.a * lam12
		sig12 = lam12 / geod.f1
		omg12 = sig12
		m12x = geod.b * math.Sin(sig12)
		if outmask&geodesicScale > 0 {
			M12 = math.Cos(sig12)
			M21 = M12
		}
		a12 = lon12 / geod.f1
	} else if !meridian {
		// Now point1 and point2 belong within a hemisphere bounded by a
		// meridian and geodesic is neither meridional or equatorial.

		// Figure a starting point for Newton's method
		var dnm float64
		sig12, salp1, calp1, salp2, calp2, dnm = geod.inverseStart(
			sbet1, cbet1, dn1, sbet2, cbet2, dn2, lam12, slam12, clam12, C1a, C2a)
		if sig12 >= 0 {
			// Short lines (InverseStart sets salp2, calp2, dnm)
			s12x = sig12 * geod.b * dnm
			m12x = math.Pow(dnm, 2) * geod.b * math.Sin(sig12/dnm)
			if outmask&geodesicScale > 0 {
				M12 = math.Cos(sig12 / dnm)
				M21 = M12
			}
			a12 = degrees(sig12)
			omg12 = lam12 / (geod.f1 * dnm)
		} else {
			// Newton's method.  This is a straightforward solution of f(alp1) =
			// lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
			// root in the interval (0, pi) and its derivative is positive at the
			// root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
			// alp1.  During the course of the iteration, a range (alp1a, alp1b) is
			// maintained which brackets the root and with each evaluation of f(alp)
			// the range is shrunk if possible.  Newton's method is restarted
			// whenever the derivative of f is negative (because the new value of
			// alp1 is then further from the solution) or if the new estimate of
			// alp1 lies outside (0,pi); in this case, the new starting guess is
			// taken to be (alp1a + alp1b) / 2.
			// real ssig1, csig1, ssig2, csig2, eps
			numit := 0
			tripn := false
			tripb := false
			// Bracketing range
			salp1a := geod.tiny
			calp1a := 1.0
			salp1b := geod.tiny
			calp1b := -1.0
			var v, ssig1, csig1, ssig2, csig2, eps, domg12, dv float64
			for numit < maxIt2 {
				// the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
				// WGS84 and random input: mean = 2.85, sd = 0.60
				v, salp2, calp2, sig12, ssig1, csig1, ssig2, csig2,
					eps, domg12, dv = geod.lambda12(
					sbet1, cbet1, dn1, sbet2, cbet2, dn2,
					salp1, calp1, slam12, clam12, numit < maxIt1,
					C1a, C2a, C3a)
				// 2 * tol0 is approximately 1 ulp for a number in [0, pi].
				// Reversed test to allow escape with NaNs
				mult := 1.0
				if tripn {
					mult = 8
				}
				if tripb || !(math.Abs(v) >= mult*geod.tol0) {
					break
				}
				// Update bracketing values
				if v > 0 && (numit > maxIt1 || calp1/salp1 > calp1b/salp1b) {
					salp1b = salp1
					calp1b = calp1
				} else if v < 0 && (numit > maxIt1 || calp1/salp1 < calp1a/salp1a) {
					salp1a = salp1
					calp1a = calp1
				}
				numit += 1
				if numit < maxIt1 && dv > 0 {
					dalp1 := -v / dv
					sdalp1 := math.Sin(dalp1)
					cdalp1 := math.Cos(dalp1)
					nsalp1 := salp1*cdalp1 + calp1*sdalp1
					if nsalp1 > 0 && math.Abs(dalp1) < math.Pi {
						calp1 = calp1*cdalp1 - salp1*sdalp1
						salp1 = nsalp1
						salp1, calp1 = norm(salp1, calp1)
						// In some regimes we don't get quadratic convergence because
						// slope -> 0.  So use convergence conditions based on epsilon
						// instead of sqrt(epsilon).
						tripn = math.Abs(v) <= 16*geod.tol0
						continue
					}
				}
				// Either dv was not positive or updated value was outside
				// legal range.  Use the midpoint of the bracket as the next
				// estimate.  This mechanism is not needed for the WGS84
				// ellipsoid, but it does catch problems with more eccentric
				// ellipsoids.  Its efficacy is such for
				// the WGS84 test set with the starting guess set to alp1 = 90deg:
				// the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
				// WGS84 and random input: mean = 4.74, sd = 0.99
				salp1 = (salp1a + salp1b) / 2
				calp1 = (calp1a + calp1b) / 2
				salp1, calp1 = norm(salp1, calp1)
				tripn = false
				tripb = (math.Abs(salp1a-salp1)+(calp1a-calp1) < geod.tolb) || (math.Abs(salp1-salp1b)+(calp1-calp1b) < geod.tolb)
			}
			var lengthMask int
			if outmask&(reducedLength|geodesicScale) > 0 {
				lengthMask = outmask | DISTANCE
			} else {
				lengthMask = outmask | empty
			}
			s12x, m12x, _, M12, M21 = geod.lengths(
				eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2, cbet1, cbet2,
				lengthMask, C1a, C2a)

			m12x *= geod.b
			s12x *= geod.b
			a12 = degrees(sig12)
			if outmask&area > 0 {
				// omg12 = lam12 - domg12
				sdomg12 := math.Sin(domg12)
				cdomg12 := math.Cos(domg12)
				somg12 = slam12*cdomg12 - clam12*sdomg12
				comg12 = clam12*cdomg12 + slam12*sdomg12
			}
		}
	}
	// end elif not meridian

	if outmask&DISTANCE > 0 {
		s12 = 0.0 + s12x // Convert -0 to 0
	}

	if outmask&reducedLength > 0 {
		m12 = 0.0 + m12x // Convert -0 to 0
	}

	if outmask&area > 0 {
		// From Lambda12: sin(alp1) * cos(bet1) = sin(alp0)
		salp0 := salp1 * cbet1
		calp0 := math.Hypot(calp1, salp1*sbet1) // calp0 > 0
		// real alp12
		if calp0 != 0 && salp0 != 0 {
			// From Lambda12: tan(bet) = tan(sig) * cos(alp)
			ssig1 := sbet1
			csig1 := calp1 * cbet1
			ssig2 := sbet2
			csig2 := calp2 * cbet2
			k2 := math.Pow(calp0, 2) * geod.ep2
			eps := k2 / (2*(1+math.Sqrt(1+k2)) + k2)
			// Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0).
			A4 := math.Pow(geod.a, 2) * calp0 * salp0 * geod.e2
			ssig1, csig1 = norm(ssig1, csig1)
			ssig2, csig2 = norm(ssig2, csig2)
			C4a := makeRange(0, nC4)
			geod.c4f(eps, C4a)
			B41 := sinCosSeries(false, ssig1, csig1, C4a)
			B42 := sinCosSeries(false, ssig2, csig2, C4a)
			S12 = A4 * (B42 - B41)
		} else {
			// Avoid problems with indeterminate sig1, sig2 on equator
			S12 = 0.0
		}
		if !meridian && somg12 > 1 {
			somg12 = math.Sin(omg12)
			comg12 = math.Cos(omg12)
		}
		var alp12 float64
		if !meridian && comg12 > -0.7071 && sbet2-sbet1 < 1.75 {
			// omg12 < 3/4 * pi and Long difference not too big and
			// and Lat difference not too big
			// Use tan(Gamma/2) = tan(omg12/2)
			// * (tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
			// with tan(x/2) = sin(x)/(1+cos(x))
			domg12 := 1 + comg12
			dbet1 := 1 + cbet1
			dbet2 := 1 + cbet2
			alp12 = 2 * math.Atan2(somg12*(sbet1*dbet2+sbet2*dbet1), domg12*(sbet1*sbet2+dbet1*dbet2))
		} else {
			// alp12 = alp2 - alp1, used in atan2 so no need to normalize
			salp12 := salp2*calp1 - calp2*salp1
			calp12 := calp2*calp1 + salp2*salp1
			// The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
			// salp12 = -0 and alp12 = -180.  However this depends on the sign
			// being attached to 0 correctly.  The following ensures the correct
			// behavior.
			if salp12 == 0 && calp12 < 0 {
				salp12 = geod.tiny * calp1
				calp12 = -1.0
			}
			alp12 = math.Atan2(salp12, calp12)
		}
		S12 += geod.c2 * alp12
		S12 *= swapp * lonsign * latsign
		// Convert -0 to 0
		S12 += 0.0
	}
	// Convert calp, salp to azimuth accounting for lonsign, swapp, latsign.
	if swapp < 0 {
		salp2, salp1 = salp1, salp2
		calp2, calp1 = calp1, calp2
		if outmask&geodesicScale > 0 {
			M21, M12 = M12, M21
		}
	}
	salp1 *= swapp * lonsign
	calp1 *= swapp * latsign
	salp2 *= swapp * lonsign
	calp2 *= swapp * latsign
	return a12, s12, salp1, calp1, salp2, calp2, m12, M12, M21, S12
}

func sinCosSeries(sinp bool, sinx, cosx float64, c []float64) float64 {
	// Evaluate a trig series using Clenshaw summation.
	// Evaluate
	// y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
	//            sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
	// using Clenshaw summation.  N.B. c[0] is unused for sin series
	// Approx operation count = (n + 5) mult and (2 * n + 2) add
	k := len(c) // Point to one beyond last element
	n := k
	if sinp {
		n -= 1
	}
	ar := 2 * (cosx - sinx) * (cosx + sinx) // 2 * cos(2 * x)
	y1 := 0.0                               // accumulators for sum
	y0 := 0.0
	if n&1 > 0 {
		// equivalent of n&1
		k -= 1
		y0 = c[k]
	}
	// Now n is even
	n = n / 2
	for n > 0 {
		n -= 1
		// Unroll loop x 2, so accumulators return to their original role
		k -= 1
		y1 = ar*y0 - y1 + c[k]
		k -= 1
		y0 = ar*y1 - y0 + c[k]
	}
	if sinp {
		// sin(2 * x) * y0
		return 2 * sinx * cosx * y0
	} else {
		// cos(x) * (y0 - y1)
		return cosx * (y0 - y1)
	}
}

func astroid(x, y float64) float64 {
	// solve astroid equation.
	// Solve k^4+2*k^3-(x^2+y^2-1)*k^2-2*y^2*k-y^2 = 0 for positive root k.
	// This solution is adapted from Geocentric::Reverse.
	p := math.Pow(x, 2)
	q := math.Pow(y, 2)
	r := (p + q - 1) / 6
	if !(q == 0 && r <= 0) {
		// Avoid possible division by zero when r = 0 by multiplying equations
		// for s and t by r^3 and r, resp.
		S := p * q / 4 // S = r^3 * s
		r2 := math.Pow(r, 2)
		r3 := r * r2
		// The discriminant of the quadratic equation for T3.  This is zero on
		// the evolute curve p^(1/3)+q^(1/3) = 1
		disc := S * (S + 2*r3)
		u := r
		if disc >= 0 {
			T3 := S + r3
			// Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
			// of precision due to cancellation.  The result is unchanged because
			// of the way the T is used in definition of u.
			if T3 < 0 {
				T3 += -math.Sqrt(disc)
			} else {
				T3 += math.Sqrt(disc) // T3 = (r * t)^3
			}
			// N.B. cbrt always returns the real root.  cbrt(-8) = -2.
			T := cbrt(T3) // T = r * t
			// T can be zero; but then r2 / T -> 0.
			u += T
			if T != 0 {
				u += r2 / T
			}
		} else {
			// T is complex, but the way u is defined the result is real.
			ang := math.Atan2(math.Sqrt(-disc), -(S + r3))
			// There are three possible cube roots.  We choose the root which
			// avoids cancellation.  Note that disc < 0 implies that r < 0.
			u += 2 * r * math.Cos(ang/3)
		}
		v := math.Sqrt(math.Pow(u, 2) + q) // guaranteed positive
		// Avoid loss of accuracy when u < 0.
		var uv float64
		if u < 0 {
			uv = q / (v - u)
		} else {
			uv = u + v // u+v, guaranteed positive
		}
		w := (uv - q) / (2 * v) // positive?
		// Rearrange expression for k to avoid loss of accuracy due to
		// subtraction.  Division by 0 not possible because uv > 0, w >= 0.
		return uv / (math.Sqrt(uv+math.Pow(w, 2)) + w) // guaranteed positive
	} else {
		// q == 0 && r <= 0
		// y = 0 with |x| <= 1.  Handle this case directly.
		// for y small, positive root is k = abs(y)/sqrt(1-x^2)
		return 0
	}
}

func a1m1f(eps float64) float64 {
	// return A1-1
	coeff := []float64{1, 4, 64, 0, 256}
	m := nA1 / 2
	t := polyval(m, coeff, 0, math.Pow(eps, 2)) / coeff[m+1]
	return (t + eps) / (1 - eps)
}

func c1f(eps float64, c []float64) {
	// return C1
	coeff := []float64{-1, 6, -16, 32, -9, 64, -128, 2048, 9, -16, 768, 3, -5, 512, -7, 1280, -7, 2048}
	eps2 := math.Pow(eps, 2)
	d := eps
	o := 0
	for l := 1; l < nC1+1; l++ { // l is index of C1p[l]
		m := (nC1 - l) / 2 // order of polynomial in eps^2
		c[l] = d * polyval(m, coeff, o, eps2) / coeff[o+m+1]
		o += m + 2
		d *= eps
	}
}

func a2m1f(eps float64) float64 {
	// return A2-1
	coeff := []float64{-11, -28, -192, 0, 256}
	m := nA2 / 2
	t := polyval(m, coeff, 0, math.Pow(eps, 2)) / coeff[m+1]
	return (t - eps) / (1 + eps)
}

func c2f(eps float64, c []float64) {
	// return C2
	coeff := []float64{1, 2, 16, 32, 35, 64, 384, 2048, 15, 80, 768, 7, 35, 512, 63, 1280, 77, 2048}
	eps2 := math.Pow(eps, 2)
	d := eps
	o := 0
	for l := 1; l < nC2+1; l++ { // l is index of C2[l]
		m := (nC2 - l) / 2 // order of polynomial in eps^2
		c[l] = d * polyval(m, coeff, o, eps2) / coeff[o+m+1]
		o += m + 2
		d *= eps
	}
}
