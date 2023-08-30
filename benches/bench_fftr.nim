import math, benchy, fftr

proc benchFFT(name: string, size: int) =
  GC_fullCollect()
  var v = newSeq[Complex64](size)
  timeIt name & " - " & $size:
    discard fft(v, false)

for size in [64, 128, 256, 512, 1024, 2048, 4096, 16384, 65536]:
  benchFFT("pow2", size)

for size in [5, 17, 149, 151, 251, 1009, 2017, 2879, 32767, 65521, 65537, 746483, 746497]:
  benchFFT "prime", size

for size in [211^2, 401^2]:
  benchFFT "prime-power", size

for size in [24576, 20736]:
  benchFFT "mult-of-power-of-2", size

for size in [30270]:
  benchFFT "small-comp-large-prime", size

for size in [18, 360, 44100, 48000, 46656, 100000]:
  benchFFT "small-comp", size
