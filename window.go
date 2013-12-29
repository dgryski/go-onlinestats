package onlinestats

// add http://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm

import "math"

type Windowed struct {
	data []float64
	head int

	length int
	sum    float64
	sumsq  float64
}

func NewWindowed(capacity int) *Windowed {
	return &Windowed{
		data: make([]float64, capacity),
	}
}

func (w *Windowed) Push(n float64) float64 {
	old := w.data[w.head]

	w.length++

	w.data[w.head] = n
	w.head++
	if w.head >= len(w.data) {
		w.head = 0
	}

	w.sum -= old
	w.sum += n

	w.sumsq -= old * old
	w.sumsq += n * n

	return old
}

func (w *Windowed) Len() int {
	if w.length < len(w.data) {
		return w.length
	}

	return len(w.data)
}

func (w *Windowed) Mean() float64 {
	return w.sum / float64(w.Len())
}

func (w *Windowed) Var() float64 {
	n := float64(w.Len())
	// population variance, not sample variance
	return (w.sumsq - (w.sum*w.sum)/n) / n
}

func (w *Windowed) Stddev() float64 {
	return math.Sqrt(w.Var())
}
