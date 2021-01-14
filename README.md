# FindPeaks1D
![CI](https://github.com/ymtoo/FindPeaks1D.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/ymtoo/FindPeaks1D.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ymtoo/FindPeaks1D.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://ymtoo.github.io/FindPeaks1D.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ymtoo.github.io/FindPeaks1D.jl/dev)

Finding peaks in a 1-D signal in Julia. The implementation is based on [`find_peaks`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html) in `SciPy`.

## Installation
```julia-repl
julia>]
pkg> add FindPeaks1D
```

## Usage
```julia
using FindPeaks1D, ImageFiltering, Plots

n = 48001
s1 = ImageFiltering.Kernel.gaussian((1000,), (n,))
s2 = ImageFiltering.Kernel.gaussian((500,), (n,))
s = s1.parent/maximum(s1.parent) + 0.5 * circshift(
    s2.parent/maximum(s2.parent), (10000,))
```
```julia
pkindices, properties = findpeaks1d(s; 
                                    height=0.1, 
                                    prominence=0.2, 
                                    width=1000.0, 
                                    relheight=0.9)

plot(s; color="black", label=false)
scatter!(pkindices, s[pkindices]; color="red", markersize=5, label="peaks")
vline!(properties["leftips"]; color="blue", width=2, label="peak edges")
vline!(properties["rightips"]; color="blue", width=2, label=false)
xlabel!("Sample")
```
![window](peaks-1.svg)
```julia
pkindices, properties = findpeaks1d(s; 
                                    height=0.1, 
                                    distance=12000, 
                                    prominence=0.2, 
                                    width=1000.0, 
                                    relheight=0.9)

plot(s; color="black", label=false)
scatter!(pkindices, s[pkindices]; color="red", markersize=5, label="peaks")
vline!(properties["leftips"]; color="blue", width=2, label="peak edges")
vline!(properties["rightips"]; color="blue", width=2, label=false)
xlabel!("Sample")
```
![window](peaks-2.svg)