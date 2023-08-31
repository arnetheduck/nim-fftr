# Package

version       = "0.1.0"
author        = "Jacek Sieka"
description   = "The fastest Fourier transform in the Rhein (so far)"
license       = "MIT"
srcDir        = "src"
installExt    = @["nim"]

requires "nim >= 1.6.0"
requires "unittest2"

taskRequires "bench", "benchy", "fftw3"

let nimc = getEnv("NIMC", "nim") # Which nim compiler to use
let lang = getEnv("NIMLANG", "c") # Which backend (c/cpp/js)
let flags = getEnv("NIMFLAGS", "") # Extra flags for the compiler
let verbose = getEnv("V", "") notin ["", "0"]
let cfg =
  " --styleCheck:usages --styleCheck:error" &
  (if verbose: "" else: " --verbosity:0 --hints:off") &
  " --outdir:build -f"

proc build(args, path: string) =
  exec nimc & " " & getPathsClause() & " " & lang & " " & cfg & " " & flags & " " & args & " " & path

proc run(args, path: string) =
  build args & " -r", path

task test, "Run tests":
  run "", "tests/test_fftr.nim"

task bench, "Benchmarks":
  run "-d:release", "benches/bench_fftr.nim"
