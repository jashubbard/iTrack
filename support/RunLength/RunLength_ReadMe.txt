ReadMe: RunLength Readme

The function RunLength applies a run-length en- and decoding:
  [b,n] = RunLength(x)
  x     = RunLength(b,n)
With b is the vector of values, n is the number of their repetitions.

You can find a lot of RLE tools in the FileExchange already. This C-Mex is
about 5 times faster than good vectorized M-versions.
The M-file RunLength_M contains vectorized and loop M-code for education.

The MEX file RunLength.c must be compiled before using:
A compiler must be installed at first (see "doc mex"). Then several
altenatives can be used for compiling:
Implicit:      Running the M-file RunLength starts a compilation.
Installer:     InstallMex('RunLength.c', 'uTest_RunLength')
Manual Win:    mex -O RunLength.c
Manual Linux:  mex -O CFLAGS="\$CFLAGS -std=c99" RunLength.c
Download:      http://www.n-simon.de/mex

Running uTest_RunLength tests validity and speed.

Jan Simon, 17-Sep-2013

