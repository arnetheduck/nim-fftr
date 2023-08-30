import std/[math, sequtils], unittest2, fftr

proc makeBuffer(size: int): seq[Complex64] =
  var tmp = newSeq[Complex64](size)
  for j in 0..tmp.high:
    tmp[j].re = sin(0.1 * TAU * float64(j))
  tmp

proc almostEqual[T](a, b: openArray[T]): bool =
  if a.len != b.len:
    return false
  for i in 0..<a.len:
    if not almostEqual(a[i].re, b[i].re):
      return false
  true

suite "FFT":
  for i in 1..512:
    test "fft vs dft " & $i:
      let
        signal = makeBuffer(i)
        aa = dft(signal, false)
        bb = fft(signal, false)

        af = foldl(aa, a + abs(b), 0.0)
        bf = foldl(bb, a + abs(b), 0.0)
      checkpoint($(aa))
      checkpoint($(bb))
      check:
        almostEqual(af, bf, 7)

    test "fft vs dft inverse " & $i:
      let
        signal = makeBuffer(i)
        aa = dft(signal, true)
        bb = fft(signal, true)

        af = foldl(aa, a + abs(b), 0.0)
        bf = foldl(bb, a + abs(b), 0.0)
      checkpoint($(aa))
      checkpoint($(bb))
      check:
        almostEqual(af, bf, 7)

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
