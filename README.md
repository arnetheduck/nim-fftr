# fftr

The fastest Fourier transform in the Rhein computes the complex [discrete Fourier transform](https://en.wikipedia.org/wiki/Discrete_Fourier_transform) (DFT) of a signal and back again with `O(N log N)` performance as usual.

Faster ones exist, but hey, it has decent speed, is less than 400loc and is written in pure Nim for minimal hassle.

It's also likely the only one developed from a viewpoint overlooking the Rhein river.

## Usage

In your library:

```nim
requires "https://github.com/arnetheduck/nim-fftr"
```

Compute frequencies in a signal:

```nim
import fftr, std/[math, sequtils]

let
  signal = (0..1023).mapIt(complex64(sin(TAU * 0.1 * float64(it))))

  # abs gets us back into real space
  # false for forward FFT, true for inverse (TODO flip it? separate names?)
  frequencies = fft(signal, false).mapIt(abs(it))
```

For performance, it's important to compile the application with [LTO](https://en.wikipedia.org/wiki/Interprocedural_optimization#WPO_and_LTO) enabled so that the `complex` module gets inlined properly - see [nim.cfg](./nim.cfg) for an example of how to do this.

## Benchmarks

Beats FFTW on at least one input size :)

[fftr](./benches/bench_fftr.nim):

```sh
nim c -d:release -r benches/bench_fftr

   min time    avg time  std dv   runs name
   0.001 ms    0.001 ms  ±0.000  x1000 pow2 - 64
   0.001 ms    0.001 ms  ±0.001  x1000 pow2 - 128
   0.003 ms    0.004 ms  ±0.002  x1000 pow2 - 256
   0.006 ms    0.007 ms  ±0.001  x1000 pow2 - 512
   0.013 ms    0.014 ms  ±0.001  x1000 pow2 - 1024
   0.029 ms    0.032 ms  ±0.001  x1000 pow2 - 2048
   0.067 ms    0.070 ms  ±0.002  x1000 pow2 - 4096
   0.309 ms    0.321 ms  ±0.004  x1000 pow2 - 16384
   1.423 ms    1.461 ms  ±0.017  x1000 pow2 - 65536
   0.001 ms    0.001 ms  ±0.000  x1000 prime - 5
   0.002 ms    0.002 ms  ±0.000  x1000 prime - 17
   0.023 ms    0.025 ms  ±0.001  x1000 prime - 149
   0.023 ms    0.025 ms  ±0.001  x1000 prime - 151
   0.025 ms    0.027 ms  ±0.001  x1000 prime - 251
   0.111 ms    0.116 ms  ±0.003  x1000 prime - 1009
   0.238 ms    0.246 ms  ±0.006  x1000 prime - 2017
   0.488 ms    0.505 ms  ±0.009  x1000 prime - 2879
   4.934 ms    5.051 ms  ±0.057   x989 prime - 32767
  10.866 ms   11.053 ms  ±0.102   x452 prime - 65521
  21.629 ms   21.957 ms  ±0.165   x228 prime - 65537
 244.610 ms  249.728 ms  ±6.073    x20 prime - 746483
 244.258 ms  246.235 ms  ±1.584    x21 prime - 746497
  10.468 ms   10.688 ms  ±0.091   x468 prime-power - 44521
  47.279 ms   47.748 ms  ±0.308   x105 prime-power - 160801
   1.309 ms    1.361 ms  ±0.021  x1000 mult-of-power-of-2 - 24576
   1.915 ms    1.959 ms  ±0.023  x1000 mult-of-power-of-2 - 20736
   4.171 ms    4.251 ms  ±0.037  x1000 small-comp-large-prime - 30270
   0.002 ms    0.002 ms  ±0.000  x1000 small-comp - 18
   0.028 ms    0.029 ms  ±0.001  x1000 small-comp - 360
   6.936 ms    7.046 ms  ±0.046   x710 small-comp - 44100
   4.653 ms    4.759 ms  ±0.041  x1000 small-comp - 48000
   4.919 ms    5.049 ms  ±0.045   x991 small-comp - 46656
  11.567 ms   11.741 ms  ±0.067   x426 small-comp - 100000
```

[fftw](./benches/bench_fftw.nim):

```sh
nim c -d:release -r benches/bench_fftw

   min time    avg time  std dv   runs name
   0.003 ms    0.005 ms  ±0.002  x1000 pow2 - 64
   0.003 ms    0.005 ms  ±0.001  x1000 pow2 - 128
   0.004 ms    0.006 ms  ±0.002  x1000 pow2 - 256
   0.005 ms    0.006 ms  ±0.002  x1000 pow2 - 512
   0.006 ms    0.007 ms  ±0.001  x1000 pow2 - 1024
   0.008 ms    0.010 ms  ±0.002  x1000 pow2 - 2048
   0.016 ms    0.024 ms  ±0.008  x1000 pow2 - 4096
   0.060 ms    0.086 ms  ±0.037  x1000 pow2 - 16384
   0.285 ms    0.471 ms  ±0.348  x1000 pow2 - 65536
   0.002 ms    0.002 ms  ±0.001  x1000 prime - 5
   0.001 ms    0.002 ms  ±0.001  x1000 prime - 17
   0.029 ms    0.040 ms  ±0.011  x1000 prime - 149
   0.038 ms    0.045 ms  ±0.008  x1000 prime - 151
   0.039 ms    0.048 ms  ±0.009  x1000 prime - 251
   0.073 ms    0.108 ms  ±0.035  x1000 prime - 1009
   0.139 ms    0.205 ms  ±0.080  x1000 prime - 2017
   0.200 ms    0.312 ms  ±0.140  x1000 prime - 2879
   0.858 ms    0.955 ms  ±0.103  x1000 prime - 32767
   4.152 ms    4.952 ms  ±1.324   x936 prime - 65521
   1.920 ms    2.050 ms  ±0.119  x1000 prime - 65537
  88.072 ms  101.897 ms ±12.072    x50 prime - 746483
  55.563 ms   63.118 ms  ±8.796    x79 prime - 746497
   1.701 ms    2.138 ms  ±0.660  x1000 prime-power - 44521
   5.052 ms    5.345 ms  ±0.471   x904 prime-power - 160801
   0.105 ms    0.129 ms  ±0.037  x1000 mult-of-power-of-2 - 24576
   0.114 ms    0.137 ms  ±0.033  x1000 mult-of-power-of-2 - 20736
   0.711 ms    0.754 ms  ±0.025  x1000 small-comp-large-prime - 30270
   0.004 ms    0.005 ms  ±0.001  x1000 small-comp - 18
   0.009 ms    0.011 ms  ±0.004  x1000 small-comp - 360
   0.230 ms    0.251 ms  ±0.028  x1000 small-comp - 44100
   0.222 ms    0.241 ms  ±0.029  x1000 small-comp - 48000
   0.270 ms    0.292 ms  ±0.022  x1000 small-comp - 46656
   0.529 ms    0.560 ms  ±0.033  x1000 small-comp - 100000
```

## Guts and references

Stuff that seemed useful while writing the code:

* [Stockham](http://wwwa.pikara.ne.jp/okojisan/otfft-en/stockham1.html)
* [Radix2](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#The_radix-2_DIT_case) FFT for power-of-2-length transforms
* [Bluestein](https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein.27s_algorithm) for prime length transforms
* [Mixed radix](https://community.arm.com/arm-community-blogs/b/graphics-gaming-and-vr-blog/posts/speeding-up-fast-fourier-transform-mixed-radix-on-mobile-arm-mali-gpu-by-means-of-opencl---part-1) of Radix-2 and Bluestein for the rest
* Butterflies for 1, 2 and 4-length transforms

More can be done of course - in particular, Radix-3,5,7 are high on the wishlist.

## Missing stuff

* planning
* multi-dimensional FFT:s
* more butterflies
* more algorithms (Rader, Good-Thomas etc)
* iterative / in-place variants
* SIMD
* etc...

## Resources

* [nimfftw3](https://github.com/SciNim/nimfftw3) wrapper for Nim - faster, less convenient to use

## Music of choice

Dooz Kawa - Etioles de Sol
