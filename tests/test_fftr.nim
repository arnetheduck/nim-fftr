import std/[math, sequtils], unittest2, fftr

proc makeBuffer(size: int): seq[Complex64] =
  var tmp = newSeq[Complex64](size)
  for j in 0..tmp.high:
    tmp[j].re = sin(0.1 * TAU * float64(j))
  tmp

suite "FFT":
  for i in 1..512:
    test "fft vs dft " & $i:
      let
        signal = makeBuffer(i)
        aa = dft(signal, false)
        bb = fft(signal, false)

        af = foldl(aa, a + b, complex64(0.0))
        bf = foldl(bb, a + b, complex64(0.0))
      checkpoint($(aa))
      checkpoint($(bb))
      checkpoint($(af))
      checkpoint($(bf))
      check:
        bf.re - af.re < 1e-10
        bf.im - af.im < 1e-10

    test "fft vs dft inverse " & $i:
      let
        signal = makeBuffer(i)
        aa = dft(signal, true)
        bb = fft(signal, true)

        af = foldl(aa, a + b, complex64(0.0))
        bf = foldl(bb, a + b, complex64(0.0))
      checkpoint($(aa))
      checkpoint($(bb))
      checkpoint($(af))
      checkpoint($(bf))
      check:
        bf.re - af.re < 1e-10
        bf.im - af.im < 1e-10

  test "normalized ifft":
    let
      signal = makeBuffer(512)
      aa = dft(dft(signal, false), true, normalize=true)
      bb = fft(fft(signal, false), true, normalize=true)

      sf = foldl(signal, a + b, complex64(0.0))
      af = foldl(aa, a + b, complex64(0.0))
      bf = foldl(bb, a + b, complex64(0.0))
    checkpoint($(aa))
    checkpoint($(bb))
    checkpoint($(af))
    checkpoint($(bf))
    check:
      abs(bf.re - af.re) < 1e-10
      abs(bf.im - af.im) < 1e-10
      abs(sf.re - bf.re) < 1e-10
      abs(sf.im - bf.im) < 1e-10

# These are broken on nim 1.6:
# static:
#   for i in [3, 4, 8, 9, 16, 17]:
#     let
#       signal = makeBuffer(i)
#       a = dft(signal, false)
#       b = fft(signal, false)
#     echo a
#     echo b
#     doAssert almostEqual(a, b)
