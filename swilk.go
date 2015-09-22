package onlinestats

import (
	"math"
	"sort"
	"strconv"
)

func SWilk(x []float64) (float64, float64, error) {

	data := make([]float64, len(x)+1)
	copy(data[1:], x)
	sort.Float64s(data[1:])
	data[0] = math.NaN()

	length := len(x)
	w, pw, err := swilkHelper(data, length, nil)
	return w, pw, err
}

// Calculate the Shapiro-Wilk W test and its significance level
// Based on the public domain code at https://joinup.ec.europa.eu/svn/sextante/soft/sextante_lib/trunk/algorithms/src/es/unex/sextante/tables/normalityTest/SWilk.java

/*
 * Constants and polynomial coefficients for swilk(). NOTE: FORTRAN counts the elements of the array x[length] as
 * x[1] through x[length], not x[0] through x[length-1]. To avoid making pervasive, subtle changes to the algorithm
 * (which would inevitably introduce pervasive, subtle bugs) the referenced arrays are padded with an unused 0th
 * element, and the algorithm is ported so as to continue accessing from [1] through [length].
 */
var C1 = []float64{math.NaN(), 0.0E0, 0.221157E0, -0.147981E0, -0.207119E1, 0.4434685E1, -0.2706056E1}
var C2 = []float64{math.NaN(), 0.0E0, 0.42981E-1, -0.293762E0, -0.1752461E1, 0.5682633E1, -0.3582633E1}
var C3 = []float64{math.NaN(), 0.5440E0, -0.39978E0, 0.25054E-1, -0.6714E-3}
var C4 = []float64{math.NaN(), 0.13822E1, -0.77857E0, 0.62767E-1, -0.20322E-2}
var C5 = []float64{math.NaN(), -0.15861E1, -0.31082E0, -0.83751E-1, 0.38915E-2}
var C6 = []float64{math.NaN(), -0.4803E0, -0.82676E-1, 0.30302E-2}
var C7 = []float64{math.NaN(), 0.164E0, 0.533E0}
var C8 = []float64{math.NaN(), 0.1736E0, 0.315E0}
var C9 = []float64{math.NaN(), 0.256E0, -0.635E-2}
var G = []float64{math.NaN(), -0.2273E1, 0.459E0}

const (
	Z90   = 0.12816E1
	Z95   = 0.16449E1
	Z99   = 0.23263E1
	ZM    = 0.17509E1
	ZSS   = 0.56268E0
	BF1   = 0.8378E0
	XX90  = 0.556E0
	XX95  = 0.622E0
	SQRTH = 0.70711E0
	TH    = 0.375E0
	SMALL = 1E-19
	PI6   = 0.1909859E1
	STQR  = 0.1047198E1
	UPPER = true
)

/**
* ALGORITHM AS R94 APPL. STATIST. (1995) VOL.44, NO.4
*
* Calculates Shapiro-Wilk normality test and P-value for sample sizes 3 <= n <= 5000 .
* Corrects AS 181, which was found to be inaccurate for n > 50.
*
* As described above with the constants, the data arrays x[] and a[] are referenced with a base element of 1 (like FORTRAN)
* instead of 0 (like Java) to avoid screwing up the algorithm. To pass in 100 data points, declare x[101] and fill elements
* x[1] through x[100] with data. x[0] will be ignored.
*
* @param x
*                Input; Data set to analyze; 100 points go in x[101] array from x[1] through x[100]
* @param n
*                Input; Number of data points in x
* @param a
*                Shapiro-Wilk coefficients.  Can be nil, or pre-computed by swilkCoeffs and passed in.
 */

type SwilkFault int

func (s SwilkFault) Error() string {
	return "swilk fault " + strconv.Itoa(int(s))
}

func swilkHelper(x []float64, n int, a []float64) (w float64, pw float64, err error) {

	if n > 5000 {
		return 0, 0, SwilkFault(2)
	}

	pw = 1.0
	if w >= 0.0 {
		w = 1.0
	}
	an := float64(n)
	if n < 3 {
		return 0, 0, SwilkFault(1)
	}

	if a == nil {
		a = SwilkCoeffs(n)
	}

	if n < 3 {
		return
	}

	// If W input as negative, calculate significance level of -W
	var w1, xx float64
	if w < 0.0 {
		w1 = 1.0 + w
	} else {

		// Check for zero range

		range_ := x[n] - x[1]
		if range_ < SMALL {
			return 0, 0, SwilkFault(6)
		}

		// Check for correct sort order on range - scaled X
		// TODO(dgryski): did the FORTRAN code puke on out-of-order X ? with ifault=7 ?
		xx = x[1] / range_
		sx := xx
		sa := -a[1]
		j := n - 1
		for i := 2; i <= n; i++ {
			xi := x[i] / range_
			// IF (XX-XI .GT. SMALL) PRINT *,' ANYTHING'
			sx += xi
			if i != j {
				sa += float64(sign(1, i-j)) * a[imin(i, j)]
			}
			xx = xi
			j--
		}

		// Calculate W statistic as squared correlation between data and coefficients
		sa /= float64(n)
		sx /= float64(n)
		ssa := 0.0
		ssx := 0.0
		sax := 0.0
		j = n
		var asa float64
		for i := 1; i <= n; i++ {
			if i != j {
				asa = float64(sign(1, i-j))*a[imin(i, j)] - sa
			} else {
				asa = -sa
			}
			xsx := x[i]/range_ - sx
			ssa += asa * asa
			ssx += xsx * xsx
			sax += asa * xsx
			j--
		}

		// W1 equals (1-W) calculated to avoid excessive rounding error
		// for W very near 1 (a potential problem in very large samples)

		ssassx := math.Sqrt(ssa * ssx)
		w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx)
	}
	w = 1.0 - w1

	// Calculate significance level for W (exact for N=3)

	if n == 3 {
		pw = PI6 * (math.Asin(math.Sqrt(w)) - STQR)
		return w, pw, nil
	}
	y := math.Log(w1)
	xx = math.Log(an)
	m := 0.0
	s := 1.0
	if n <= 11 {
		gamma := poly(G, 2, an)
		if y >= gamma {
			pw = SMALL
			return w, pw, nil
		}
		y = -math.Log(gamma - y)
		m = poly(C3, 4, an)
		s = math.Exp(poly(C4, 4, an))
	} else {
		m = poly(C5, 4, xx)
		s = math.Exp(poly(C6, 3, xx))
	}
	pw = alnorm((y-m)/s, UPPER)

	return w, pw, nil
}

// Precomputes the coefficients array a for SWilk
func SwilkCoeffs(n int) []float64 {

	a := make([]float64, n+1)

	an := float64(n)

	n2 := n / 2

	if n == 3 {
		a[1] = SQRTH
	} else {
		an25 := an + 0.25
		summ2 := 0.0
		for i := 1; i <= n2; i++ {
			a[i] = ppnd((float64(i) - TH) / an25)
			summ2 += a[i] * a[i]
		}
		summ2 *= 2.0
		ssumm2 := math.Sqrt(summ2)
		rsn := 1.0 / math.Sqrt(an)
		a1 := poly(C1, 6, rsn) - a[1]/ssumm2

		// Normalize coefficients

		var i1 int
		var fac float64
		if n > 5 {
			i1 = 3
			a2 := -a[2]/ssumm2 + poly(C2, 6, rsn)
			fac = math.Sqrt((summ2 - 2.0*a[1]*a[1] - 2.0*a[2]*a[2]) / (1.0 - 2.0*a1*a1 - 2.0*a2*a2))
			a[1] = a1
			a[2] = a2
		} else {
			i1 = 2
			fac = math.Sqrt((summ2 - 2.0*a[1]*a[1]) / (1.0 - 2.0*a1*a1))
			a[1] = a1
		}
		for i := i1; i <= n2; i++ {
			a[i] = -a[i] / fac
		}
	}

	return a
}

/**
 * Constructs an int with the absolute value of x and the sign of y
 *
 * @param x
 *                int to copy absolute value from
 * @param y
 *                int to copy sign from
 * @return int with absolute value of x and sign of y
 */
func sign(x int, y int) int {
	var result = x
	if x < 0 {
		result = -x
	}
	if y < 0 {
		result = -result
	}
	return result
}

// Constants & polynomial coefficients for ppnd(), slightly renamed to avoid conflicts. Could define
// them inside ppnd(), but static constants are more efficient.

// Coefficients for P close to 0.5
const (
	A0_p = 3.3871327179E+00
	A1_p = 5.0434271938E+01
	A2_p = 1.5929113202E+02
	A3_p = 5.9109374720E+01
	B1_p = 1.7895169469E+01
	B2_p = 7.8757757664E+01
	B3_p = 6.7187563600E+01

	// Coefficients for P not close to 0, 0.5 or 1 (names changed to avoid conflict with swilk())
	C0_p = 1.4234372777E+00
	C1_p = 2.7568153900E+00
	C2_p = 1.3067284816E+00
	C3_p = 1.7023821103E-01
	D1_p = 7.3700164250E-01
	D2_p = 1.2021132975E-01

	// Coefficients for P near 0 or 1.
	E0_p = 6.6579051150E+00
	E1_p = 3.0812263860E+00
	E2_p = 4.2868294337E-01
	E3_p = 1.7337203997E-02
	F1_p = 2.4197894225E-01
	F2_p = 1.2258202635E-02

	SPLIT1 = 0.425
	SPLIT2 = 5.0
	CONST1 = 0.180625
	CONST2 = 1.6
)

/**
 * ALGORITHM AS 241 APPL. STATIST. (1988) VOL. 37, NO. 3, 477-484.
 *
 * Produces the normal deviate Z corresponding to a given lower tail area of P; Z is accurate to about 1 part in 10**7.
 *
 * @param p
 * @return
 */
func ppnd(p float64) float64 {
	q := p - 0.5
	var r float64
	if math.Abs(q) <= SPLIT1 {
		r = CONST1 - q*q
		return q * (((A3_p*r+A2_p)*r+A1_p)*r + A0_p) / (((B3_p*r+B2_p)*r+B1_p)*r + 1.0)
	} else {
		if q < 0.0 {
			r = p
		} else {
			r = 1.0 - p
		}
		if r <= 0.0 {
			return 0.0
		}
		r = math.Sqrt(-math.Log(r))
		var normal_dev float64
		if r <= SPLIT2 {
			r -= CONST2
			normal_dev = (((C3_p*r+C2_p)*r+C1_p)*r + C0_p) / ((D2_p*r+D1_p)*r + 1.0)
		} else {
			r -= SPLIT2
			normal_dev = (((E3_p*r+E2_p)*r+E1_p)*r + E0_p) / ((F2_p*r+F1_p)*r + 1.0)
		}
		if q < 0.0 {
			normal_dev = -normal_dev
		}
		return normal_dev
	}
}

/**
 * Algorithm AS 181.2 Appl. Statist. (1982) Vol. 31, No. 2
 *
 * Calculates the algebraic polynomial of order nord-1 with array of coefficients c. Zero order coefficient is c[1]
 *
 * @param c
 * @param nord
 * @param x
 * @return
 */
func poly(c []float64, nord int, x float64) float64 {
	poly := c[1]
	if nord == 1 {
		return poly
	}
	p := x * c[nord]
	if nord != 2 {
		n2 := nord - 2
		j := n2 + 1
		for i := 1; i <= n2; i++ {
			p = (p + c[j]) * x
			j--
		}
	}
	poly += p
	return poly
}

// Constants & polynomial coefficients for alnorm(), slightly renamed to avoid conflicts.
const (
	CON_a    = 1.28
	LTONE_a  = 7.0
	UTZERO_a = 18.66
	P_a      = 0.398942280444
	Q_a      = 0.39990348504
	R_a      = 0.398942280385
	A1_a     = 5.75885480458
	A2_a     = 2.62433121679
	A3_a     = 5.92885724438

	B1_a = -29.8213557807
	B2_a = 48.6959930692

	C1_a = -3.8052E-8
	C2_a = 3.98064794E-4
	C3_a = -0.151679116635
	C4_a = 4.8385912808
	C5_a = 0.742380924027
	C6_a = 3.99019417011

	D1_a = 1.00000615302
	D2_a = 1.98615381364
	D3_a = 5.29330324926
	D4_a = -15.1508972451
	D5_a = 30.789933034
)

/**
 * Algorithm AS66 Applied Statistics (1973) vol.22, no.3
 *
 * Evaluates the tail area of the standardised normal curve from x to infinity if upper is true or from minus infinity to x if
 * upper is false.
 */
func alnorm(x float64, upper bool) float64 {
	up := upper
	z := x
	if z < 0.0 {
		up = !up
		z = -z
	}
	var fn_val float64
	if z > LTONE_a && (!up || z > UTZERO_a) {
		fn_val = 0.0
	} else {
		y := 0.5 * z * z
		if z <= CON_a {
			fn_val = 0.5 - z*(P_a-Q_a*y/(y+A1_a+B1_a/(y+A2_a+B2_a/(y+A3_a))))
		} else {
			fn_val = R_a * math.Exp(-y) / (z + C1_a + D1_a/(z+C2_a+D2_a/(z+C3_a+D3_a/(z+C4_a+D4_a/(z+C5_a+D5_a/(z+C6_a))))))
		}
	}
	if !up {
		fn_val = 1.0 - fn_val
	}
	return fn_val
}

func imin(i, j int) int {
	if i < j {
		return i
	}
	return j
}
