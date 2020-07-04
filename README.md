# PlasmaDispersionFunctions.jl

![CI](https://github.com/jwscook/PlasmaDispersionFunctions.jl/workflows/CI/badge.svg)
[![codecov.io](http://codecov.io/github/jwscook/PlasmaDispersionFunctions.jl/coverage.svg?branch=master)](http://codecov.io/github/jwscook/PlasmaDispersionFunctions.jl?branch=master)

Julia implementation of the plasma dispersion function

<a href="https://www.codecogs.com/eqnedit.php?latex=\it{Z}_n(z)=\frac{1}{\sqrt{\pi}}\int_{\infty}^{\infty}\frac{x^n&space;e^{-x^2}}{x-z}dx" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\it{Z}_n(z)=\frac{1}{\sqrt{\pi}}\int_{\infty}^{\infty}\frac{x^n&space;e^{-x^2}}{x-z}dx" title="\it{Z}_n(z)=\frac{1}{\sqrt{\pi}}\int_{\infty}^{\infty}\frac{x^n e^{-x^2}}{x-z}dx" /></a>,

with a recurrence relation

<a href="https://www.codecogs.com/eqnedit.php?latex={\it&space;Z}_n(x)&space;=&space;x&space;{\it&space;Z}_{n-1}(x)&space;&plus;&space;\delta_{0,&space;n&space;\mathrm{mod}&space;2}&space;2^{\frac{1-n}{2}}&space;\Pi_{i=1}^{\frac{n-3}{2}}&space;\left(&space;2i&space;&plus;&space;1&space;\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?{\it&space;Z}_n(x)&space;=&space;x&space;{\it&space;Z}_{n-1}(x)&space;&plus;&space;\delta_{0,&space;n&space;\mathrm{mod}&space;2}&space;2^{\frac{1-n}{2}}&space;\Pi_{i=1}^{\frac{n-3}{2}}&space;\left(&space;2i&space;&plus;&space;1&space;\right)" title="{\it Z}_n(z) = z {\it Z}_{n-1}(z) + \delta_{0, n \mathrm{mod} 2} 2^{\frac{1-n}{2}} \Pi_{i=1}^{\frac{n-3}{2}} \left( 2i + 1 \right)" /></a>.

and

<a href="https://www.codecogs.com/eqnedit.php?latex={\it&space;Z}_0(z)&space;=&space;i&space;\sqrt{\pi}&space;\mathrm{erfcx}&space;\left(-i&space;z&space;\right)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?{\it&space;Z}_0(z)&space;=&space;i&space;\sqrt{\pi}&space;\mathrm{erfcx}&space;\left(-i&space;z&space;\right)" title="{\it Z}_0(z) = i \sqrt{\pi} \mathrm{erfcx} \left(-i z \right)" /></a>
