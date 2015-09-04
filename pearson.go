package onlinestats

import "math"

func Pearson(a, b []float64) float64 {

	if len(a) != len(b) {
		panic("len(a) != len(b)")
	}

	n := float64(len(a))

	var abar, bbar float64
	for i := range a {
		abar += a[i]
		bbar += b[i]
	}
	abar, bbar = abar/n, bbar/n

	var numerator float64
	var sumAA, sumBB float64

	for i := range a {
		numerator += (a[i] - abar) * (b[i] - bbar)
		sumAA += (a[i] - abar) * (a[i] - abar)
		sumBB += (b[i] - bbar) * (b[i] - bbar)
	}

	return numerator / (math.Sqrt(sumAA) * math.Sqrt(sumBB))
}
