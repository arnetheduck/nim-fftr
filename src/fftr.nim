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

const SQRT_2 = 1.41421356237309504880

func rotate45(value: Complex64, inverse: bool): Complex64 =
  if not inverse:
    complex64((value.re + value.im) / SQRT_2, (value.im - value.re) / SQRT_2)
  else:
    complex64((value.re - value.im) / SQRT_2, (value.im + value.re) / SQRT_2)

func rotate90(value: Complex64, inverse: bool): Complex64 =
  if not inverse:
    complex64(value.im, -value.re)
  else:
    complex64(-value.im, value.re)

func rotate180(value: Complex64): Complex64 =
  complex64(-value.re, -value.im)

func twiddleTerm(fftLen: int, inverse: bool): float64 =
  # the constant term of the twiddle computation
  if inverse:
    TAU / float64(fftLen)
  else:
    -TAU / float64(fftLen)

func twiddle(i: int, term: float64): Complex64 =
  let theta = float64(i) * term
  complex64(cos(theta), sin(theta))

func computeTwiddles1(len: int, inverse: bool): seq[Complex64] =
  let term = twiddleTerm(len, inverse)

  result = newSeq[Complex64](len)
  for i in 0..<len:
    result[i] = twiddle(i, term)

func computeTwiddles2(len: int, inverse: bool): seq[Complex64] =
  let term = twiddleTerm(len, inverse)

  result = newSeq[Complex64](len)
  for i in 0..<(len div 2):
    let tmp = twiddle(i, term)
    result[i] = tmp
    result[i + len div 2] = rotate180(tmp)

func computeTwiddles4(len: int, inverse: bool): seq[Complex64] =
  let term = twiddleTerm(len, inverse)

  result = newSeq[Complex64](len)
  for i in 0..<(len div 4):
    let tmp = twiddle(i, term)
    result[i] = tmp
    result[i + len div 4] = rotate90(tmp, inverse)
    result[i + len div 2] = rotate180(tmp)
    result[i + 3 * len div 4] = rotate90(tmp, not inverse)

func computeTwiddles8(len: int, inverse: bool): seq[Complex64] =
  let term = twiddleTerm(len, inverse)

  result = newSeq[Complex64](len)
  for i in 0..<(len div 8):
    let tmp = twiddle(i, term)
    let tmp45 = rotate45(tmp, inverse)
    result[i] = tmp
    result[i + len div 4] = rotate90(tmp, inverse)
    result[i + len div 2] = rotate180(tmp)
    result[i + 3 * len div 4] = rotate90(tmp, not inverse)
    let i45 = i + len div 8
    result[i45] = tmp45
    result[i45 + len div 4] = rotate90(tmp45, inverse)
    result[i45 + len div 2] = rotate180(tmp45)
    result[i45 + 3 * len div 4] = rotate90(tmp45, not inverse)

func computeTwiddles(len: int, inverse: bool): seq[Complex64] =
  if len mod 8 == 0: computeTwiddles8(len, inverse)
  elif len mod 4 == 0: computeTwiddles4(len, inverse)
  elif len mod 2 == 0: computeTwiddles2(len, inverse)
  else: computeTwiddles1(len, inverse)

func butterfly1(i0: Complex64, o0: var Complex64) =
  o0 = i0

func butterfly2(i0, i1: Complex64, o0, o1: var Complex64) =
  o0 = i0 + i1
  o1 = i0 - i1

func butterfly4(
    i0, i1, i2, i3: Complex64, o0, o1, o2, o3: var Complex64,
    inverse: bool) =
  var a0, a1, a2, a3: Complex64
  butterfly2(i0, i2, a0, a2)
  butterfly2(i1, i3, a1, a3)

  a3 = rotate90(a3, inverse)

  butterfly2(a0, a1, o0, o2)
  butterfly2(a2, a3, o1, o3)

type
  Stockham2* = object
    twiddles: seq[Complex64]
    scratch: seq[Complex64]
    inverse: bool

func stockham2(
    ctx: var Stockham2, stride, n: int,
    input: ptr UncheckedArray[Complex64],
    output: ptr UncheckedArray[Complex64],
    tmp: ptr UncheckedArray[Complex64]) =
  # Stockham FFT for powers of two
  assert isPowerOfTwo(n)

  case n
  of 0: discard
  of 1: butterfly1(input[0], output[0])
  of 2:
    for q in 0..<stride:
      butterfly2(input[q], input[q + stride], output[q], output[q + stride])
  of 4:
    for q in 0..<stride:
      butterfly4(
          input[q], input[q + stride], input[q + 2 * stride], input[q + 3 * stride],
          output[q], output[q + stride], output[q + 2 * stride], output[q + 3 * stride],
          ctx.inverse)
  else:
    let m = n div 2

    block: # twiddle is a noop on the first round
      for q in 0..<stride:
        const p = 0
        let
          a = input[q + stride * p]
          b = input[q + stride * (p + m)]
        output[q + stride * (2 * p)] = a + b
        output[q + stride * ((2 * p) + 1)] = (a - b) # * wp

    for p in 1..<m:
      let wp = ctx.twiddles[p * stride]

      for q in 0..<stride:
        let
          a = input[q + stride * p]
          b = input[q + stride * (p + m)]
        output[q + stride * (2 * p)] = a + b
        output[q + stride * ((2 * p) + 1)] = (a - b) * wp

    stockham2(ctx, 2 * stride, m, output, tmp, output);

func process*(
    ctx: var Stockham2, input: openArray[Complex64],
    output: var openArray[Complex64]) =
  let
    direct = input.len.countTrailingZeroBits mod 2 == 0 or input.len <= 4
    o = cast[ptr UncheckedArray[Complex64]](addr output[0])
    tmp = cast[ptr UncheckedArray[Complex64]](addr ctx.scratch[0])
  # stockham2 will flip between o and tmp for every iteration so we want the final
  # result to end up in the output
  stockham2(
    ctx, 1, input.len,
    cast[ptr UncheckedArray[Complex64]](unsafeAddr input[0]),
    if direct: o else: tmp,
    if direct: tmp else: o
    )

func init(T: type Stockham2, len: int, inverse: bool): T =
  T(
    twiddles: if len >= 4: computeTwiddles(len, inverse) else: @[],
    scratch: newSeq[Complex64](len),
    inverse: inverse)

type
  Bluestein* = object
    len: int
    bk: seq[Complex64]
    multiplier: seq[Complex64]
    inner: Stockham2
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

  var inner = Stockham2.init(scratch.len, inverse)
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
    output[i] = conjugate(ctx.scratch[i]) * conjugate(ctx.bk[i])

type
  RadixMixed* = object
    len: int
    twiddles: seq[Complex64]
    scratch: seq[Complex64]

    inner1: Stockham2
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
    inner1: Stockham2.init(height, inverse),
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
  # Mixed-radix 6-step FFT
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

type
  FftNormalization* = enum
    ## Normalization modes for the FFT and the IFFT
    backward,  ## Only the IFFT (i.e. the backward transform) is scaled by `1/n`.
    ortho,     ## Both the FFT and the IFFT are scaled by `1/sqrt(n)`.
    forward,   ## Only the FFT (i.e. the forward transform) is scaled by `1/n`.
    disabled   ## No normalization is applied to any transform.

func dftImpl(input: openArray[Complex64], inverse: bool = false, norm: FftNormalization): seq[Complex64] =
  ## Slow DFT implementation (useful for testing)
  let size = input.len
  result.setLen(size)

  let twiddles = computeTwiddles(size, inverse)
  for k in 0..<result.len:
    for n in 0..<result.len:
      result[k] += input[n] * twiddles[k * n mod twiddles.len]

  if norm == disabled or (inverse and norm == forward) or (not inverse and norm == backward):
    return

  let scalingFactor = if norm == ortho:
      sqrt(float64(size))
    else:
      float64(size)
  for i in 0..<result.len:
    result[i] = result[i] / scalingFactor

func dft*(input: openArray[Complex64], inverse: bool = false, norm: FftNormalization = backward, n = 0): seq[Complex64] =
  ## Calculates the DFT or the IDFT of an input signal using a slow, straight-forward implementation
  ##
  ## This function is only meant for testing, since it is a straight-forward but slow implementaiton of the DFT
  ##
  ## Inputs:
  ## - input: The input signal.
  ## - inverse: If true, calculates the inverse FFT, otherwise calculates the FFT (which is the default).
  ## - norm: Normalization mode. Default is `backward`.
  ##         For norm=`backward` no normalization is applied to the FFT
  ##         (i.e. to the forward transform) and scaling by `1/n` is applied to the IFFT
  ##         (i.e. to the backward transform).
  ##         For norm=`forward` a `1/n` factor is only applied to the FFT (i.e. to the forward tranform).
  ##         For norm=`ortho`, both transforms are scaled by `1/sqrt(n)`.
  ##         For norm=`disabled`, no normalization is applied to any transform.
  ## - n: Length of the the output.
  ##      If `n` is not given, the length of the input is automatically used.
  ##      If `n` is smaller than the length of the input, the input is cropped.
  ##      If `n` is larger, the input is padded with zeros.
  let size = if n <= 0: input.len else: n

  if size > input.len:
    # Extend the input with zeros as much as needed
    var extendedInput = newSeq[Complex64](size)
    extendedInput[0..input.high] = input
    dftImpl(extendedInput, inverse=inverse, norm=norm)
  elif size < input.len:
    # Crop the input as much as needed
    dftImpl(input[0..<size], inverse=inverse, norm=norm)
  else:
    dftImpl(input, inverse=inverse, norm=norm)

func fftImpl(input: openArray[Complex64], inverse: bool = false, norm: FftNormalization): seq[Complex64] =
  ## Calculates the FFT or the IFFT of an input signal using the best method given the input length
  ##
  ## Inputs:
  ## - input: The input signal.
  ## - inverse: if true, calculates the inverse FFT, otherwise calculates the FFT.
  ## - norm: Normalization mode. Default is `backward`.
  ##         For norm=`backward` no normalization is applied to the FFT
  ##         (i.e. to the forward transform) and scaling by `1/n` is applied to the IFFT
  ##         (i.e. to the backward transform).
  ##         For norm=`forward` a `1/n` factor is only applied to the FFT (i.e. to the forward tranform).
  ##         For norm=`ortho`, both transforms are scaled by `1/sqrt(n)`.
  ##         For norm=`disabled`, no normalization is applied to any transform.
  let size = input.len
  result.setLen(size)

  if isPowerOfTwo(size):
    var ctx = Stockham2.init(size, inverse)
    process(ctx, input, result)
  elif size.countTrailingZeroBits == 0:
    var ctx = Bluestein.init(size, inverse)
    process(ctx, input, result)
  else:
    var ctx = RadixMixed.init(size, inverse)
    process(ctx, input, result)

  if norm == disabled or (inverse and norm == forward) or (not inverse and norm == backward):
    return

  let scalingFactor = if norm == ortho:
      sqrt(float64(size))
    else:
      float64(size)
  for i in 0..<result.len:
    result[i] = result[i] / scalingFactor

func fft*(input: openArray[Complex64], inverse: bool = false, norm: FftNormalization = backward, n = 0): seq[Complex64] =
  ## Calculates the FFT or the IFFT of an input signal
  ##
  ## This function automatically selects the best calculation method given the desired FFT size
  ##
  ## Inputs:
  ## - input: The input signal.
  ## - inverse: If true, calculates the inverse FFT, otherwise calculates the FFT (which is the default).
  ## - norm: Normalization mode. Default is `backward`.
  ##         For norm=`backward` no normalization is applied to the FFT
  ##         (i.e. to the forward transform) and scaling by `1/n` is applied to the IFFT
  ##         (i.e. to the backward transform).
  ##         For norm=`forward` a `1/n` factor is only applied to the FFT (i.e. to the forward tranform).
  ##         For norm=`ortho`, both transforms are scaled by `1/sqrt(n)`.
  ##         For norm=`disabled`, no normalization is applied to any transform.
  ## - n: Length of the the output.
  ##      If `n` is not given, the length of the input is automatically used.
  ##      If `n` is smaller than the length of the input, the input is cropped.
  ##      If `n` is larger, the input is padded with zeros.
  let size = if n <= 0: input.len else: n

  if size > input.len:
    # Extend the input with zeros as much as needed
    var extendedInput = newSeq[Complex64](size)
    extendedInput[0..input.high] = input
    fftImpl(extendedInput, inverse=inverse, norm=norm)
  elif size < input.len:
    # Crop the input as much as needed
    fftImpl(input[0..<size], inverse=inverse, norm=norm)
  else:
    fftImpl(input, inverse=inverse, norm=norm)

func ifft*(input: openArray[Complex64], n = 0, norm: FftNormalization = backward): seq[Complex64] {.inline.} =
  ## Calculates the IFFT of an input signal using the best method given the input length
  ##
  ## This is a convenience function that simply calls `fft` with the inverse argument set to true.
  ## 
  ## Inputs:
  ## - input: The input signal.
  ## - norm: Normalization mode. Default is `backward`.
  ##         For norm=`backward` no normalization is applied to the FFT
  ##         (i.e. to the forward transform) and scaling by `1/n` is applied to the IFFT
  ##         (i.e. to the backward transform).
  ##         For norm=`forward` a `1/n` factor is only applied to the FFT (i.e. to the forward tranform).
  ##         For norm=`ortho`, both transforms are scaled by `1/sqrt(n)`.
  ##         For norm=`disabled`, no normalization is applied to any transform.
  ## - n: Length of the the output.
  ##      If `n` is not given, the length of the input is automatically used.
  ##      If `n` is smaller than the length of the input, the input is cropped.
  ##      If `n` is larger, the input is padded with zeros.
  fft(input, inverse=true, norm=norm, n=n)
