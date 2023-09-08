## Discrete Fourier transform implementation(s)
##
## Usage:
##
## ```nim
## import fftr, std/math, std/sequtils
##
## let
##   signal = (0..1023).mapIt(
##     complex64(sin(TAU * 0.1 * float64(it))))
##
##   freqs = fft(tmp, false) # False for forward, true for inverse FFT
## ```

import
  std/[bitops, complex, math]

export complex

func twiddleTerm(fftLen: int, inverse: bool): float64 =
  # the constant term of the twiddle computation
  if inverse:
    TAU / float64(fftLen)
  else:
    -TAU / float64(fftLen)

func twiddle(i: int, term: float64): Complex64 =
  let theta = float64(i) * term
  complex64(cos(theta), sin(theta))

func rotate90(value: Complex64, inverse: bool): Complex64 =
  if not inverse:
    complex64(value.im, -value.re)
  else:
    complex64(-value.im, value.re)

func computeTwiddles(len: int, inverse: bool): seq[Complex64] =
  let term = twiddleTerm(len, inverse)

  result = newSeq[Complex64](len)
  for i in 0..<len:
    result[i] = twiddle(i, term)

func butterfly1(input: openArray[Complex64], output: var openArray[Complex64]) =
  assert input.len == 1
  assert input.len == output.len
  output[0] = input[0]

func butterfly2(i0, i1: Complex64, o0, o1: var Complex64) =
  o0 = i0 + i1
  o1 = i0 - i1

func butterfly2(input: openArray[Complex64], output: var openArray[Complex64]) =
  assert input.len == 2
  assert input.len == output.len
  butterfly2(input[0], input[1], output[0], output[1])

func butterfly4(
    input: openArray[Complex64], output: var openArray[Complex64],
    inverse: bool) =
  assert input.len == 4
  assert input.len == output.len

  var a0, a1, a2, a3: Complex64
  butterfly2(input[0], input[2], a0, a2)
  butterfly2(input[1], input[3], a1, a3)

  a3 = rotate90(a3, inverse)

  butterfly2(a0, a1, output[0], output[2])
  butterfly2(a2, a3, output[1], output[3])

template butterflies: bool =
  case input.len
  of 0:
    true
  of 1:
    butterfly1(input, output)
    true
  of 2:
    butterfly2(input, output)
    true
  of 4:
    butterfly4(input, output, ctx.inverse)
    true
  else:
    false

type
  Radix2* = object
    twiddles: seq[Complex64]
    scratch: seq[Complex64]
    inverse: bool

func radix2(
    ctx: var Radix2, strideShift: int,
    input: openArray[Complex64], output: var openArray[Complex64]) =
  # Simple divide-and-conquer FFT for powers of 2
  assert isPowerOfTwo(input.len)
  assert input.len == output.len

  if butterflies():
    return

  let
    len = input.len
    half = len div 2

  # Avoid overwriting scratch space during recursion
  # TODO iterative / inplace algo instead?
  for i in 0..<half:
    ctx.scratch[i + half] = input[i * 2]
  radix2(
    ctx, strideShift + 1, ctx.scratch.toOpenArray(half, len - 1),
    output.toOpenArray(0, half-1))

  for i in 0..<half:
    ctx.scratch[i + half] = input[i * 2 + 1]

  radix2(
    ctx, strideShift + 1, ctx.scratch.toOpenArray(half, len - 1),
    output.toOpenArray(half, output.high))

  for i in 0..<half:
    let
      p = output[i]
      q = ctx.twiddles[i shl strideShift] * output[i + half]
    output[i] = p + q
    output[i + half] = p - q

func process*(
    ctx: var Radix2, input: openArray[Complex64],
    output: var openArray[Complex64]) =
  radix2(ctx, 0, input, output)

func init(T: type Radix2, len: int, inverse: bool): T =
  T(
    twiddles: computeTwiddles(len, inverse),
    scratch: newSeq[Complex64](len),
    inverse: inverse)

type
  Bluestein* = object
    len: int
    bk: seq[Complex64]
    multiplier: seq[Complex64]
    inner: Radix2
    scratch, scratch2: seq[Complex64]

func init*(T: type Bluestein, len: int, inverse: bool): T =
  let
    twiddles = computeTwiddles(len * 2, not inverse)
  var
    bk = newSeq[Complex64](len)
  bk[0].re = 1.0
  var coeff = 0
  for m in 1..<len:
    coeff += 2*m - 1
    if coeff >= 2*len: coeff -= 2*len
    bk[m] = twiddles[coeff]

  let
    innerLen = nextPowerOfTwo(2 * len - 1)
    innerScale = 1.0 / float64(innerLen)

  var
    scratch = newSeq[Complex64](innerLen)

  scratch[0] = bk[0] * innerScale
  for i in 1..<len:
    let t = bk[i] * innerScale
    scratch[i] = t
    scratch[scratch.len - i] = t

  var inner = Radix2.init(scratch.len, inverse)
  var multiplier = newSeq[Complex64](scratch.len)
  process(inner, scratch, multiplier)

  T(
    len: len,
    bk: bk,
    multiplier: multiplier,
    inner: inner,
    scratch: scratch,
    scratch2: newSeq[Complex64](scratch.len)
  )

func process*(
    ctx: var Bluestein, input: openArray[Complex64],
    output: var openArray[Complex64]) =
  # Less efficient FFT for any length
  assert input.len == output.len

  for i in 0..<input.len:
    ctx.scratch[i] = input[i] * conjugate(ctx.bk[i])

  for i in input.len..<ctx.scratch.len:
    ctx.scratch[i].reset()

  process(ctx.inner, ctx.scratch, ctx.scratch2)

  for i in 0..<ctx.scratch2.len:
    # The conjugate inverts the inner fft below
    ctx.scratch2[i] = conjugate(ctx.scratch2[i] * ctx.multiplier[i])

  process(ctx.inner, ctx.scratch2, ctx.scratch)

  for i in 0..<input.len:
    output[i] = ctx.scratch[i] * ctx.bk[i]

type
  RadixMixed* = object
    len: int
    twiddles: seq[Complex64]
    scratch: seq[Complex64]

    inner1: Radix2
    inner2: Bluestein

func init*(T: type RadixMixed, len: int, inverse: bool): T =
  let
    height = 1 shl len.countTrailingZeroBits()
    width = len div height
    term = twiddleTerm(len, inverse)
    twiddles = block:
      var tmp = newSeq[Complex64](len)
      for x in 0..<width:
        for y in 0..<height:
          tmp[x * height + y] = twiddle(x*y, term)
      tmp
  T(
    len: len,
    twiddles: twiddles,
    scratch: newSeq[Complex64](len),
    inner1: Radix2.init(height, inverse),
    inner2: Bluestein.init(width, inverse)
  )

func transpose[T](
    input: openArray[T], output: var openArray[T], width, height: int) =
  assert input.len == output.len

  for x in 0..<width:
    for y in 0..<height:
      output[y + x * height] = input[x + y * width]

func process*(
    ctx: var RadixMixed, input: openArray[Complex64],
    output: var openArray[Complex64]) =
  # Mixed-radix FFT
  assert input.len == output.len

  let
    height = 1 shl input.len.countTrailingZeroBits()
    width = input.len div height

  transpose(input, output, width, height)

  block:
    for i in 0..<width:
      process(
        ctx.inner1,
        output.toOpenArray(i*height, (i+1)*height - 1),
        ctx.scratch.toOpenArray(i*height, (i+1)*height - 1))

  for i in 0..<ctx.scratch.len:
    ctx.scratch[i] = ctx.scratch[i] * ctx.twiddles[i]

  transpose(ctx.scratch, output, height, width)

  block:
    for i in 0..<height:
      process(
        ctx.inner2,
        output.toOpenArray(i*width, (i+1)*width - 1),
        ctx.scratch.toOpenArray(i*width, (i+1)*width - 1))

  transpose(ctx.scratch, output, width, height)

func dft*(input: openArray[Complex64], inverse: bool): seq[Complex64] =
  # Slow DFT - useful for testing
  result.setLen(input.len)

  let twiddles = computeTwiddles(input.len, inverse)
  for k in 0..<result.len:
    for n in 0..<result.len:
      result[k] += input[n] * twiddles[k * n mod twiddles.len]

func fft*(input: openArray[Complex64], inverse: bool): seq[Complex64] =
  # Convenience entry point that chooses FFT based on input length
  result.setLen(input.len)

  if isPowerOfTwo(input.len):
    var ctx = Radix2.init(input.len, inverse)
    process(ctx, input, result)
  elif input.len.countTrailingZeroBits == 0:
    var ctx = Bluestein.init(input.len, inverse)
    process(ctx, input, result)
  else:
    var ctx = RadixMixed.init(input.len, inverse)
    process(ctx, input, result)
