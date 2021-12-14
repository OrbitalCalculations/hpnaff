# hpNaff
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5780013.svg)](https://doi.org/10.5281/zenodo.5780013)

## Numerical Analysis of Fundamental Frequencies (NAFF)
A little command line tool (macOS, Linux, Windows[^1]) to perform high-resolution frequency analysis of time series following the Laskar method.

Pre-compiled binaries are available [here](https://paloz.marum.de/bitbucket/projects/ESSP/repos/hpnaff/browse/Binaries)

For macOS: after download, copy 'hpnaff' to a location in your search path, and  add execute permissions by 'chmod a+x ./hpnaff; cp hpnaff /usr/local/bin' 
Then try it on an [example file](https://paloz.marum.de/confluence/download/attachments/30048264/etp_1kyr_36Ma.txt?version=1&modificationDate=1551703946466&api=v2) that contains a weighted mix of Earth's orbital eccentricity, obliquity and climatic precession.

[^1]: Windows version requires installation of [Windows Swift Toolchain](https://www.swift.org/download/)

## Usage:
```
% ./hpnaff -v
Welcome to hpNAFF Analysis. ©2015-2021 Heiko Pälike. 3rd party licenses: see --licenses option
Please provide exactly one file name argument
Usage: hpnaff [options] file
file format: tab separated columns, data must be equidistantly spaced. Instead of file can specify '-' for stdin.
  -h, --help:
      Displays help message.
  -r, --fracRMSChangeLimit:
      fracRMSChangeLimit. Default: 0.000002
  -f, --fracFreqAccuracyLimit:
      fracFreqAccuracyLimit. Default: 0.00000001
  -m, --maxFrequencies:
      maxFrequencies. Default: 30
  -c, --freqCycleLimit:
      freqCycleLimit. Default: 100
  -l, --lowerFreqLimit:
      lowerFreqLimit. Default: 0.0
  -u, --upperFreqLimit:
      upperFreqLimit. Default: 0.5 (NYQUIST)
  -d, --detrendingOrder:
      Order for polynomial detrending. 0 removes DC only. >0 removes a polynomial fit. Default: not set
  -dC, --detrendChunks:
      Option to apply detrending (-d) to individual chunks. Default: false
  -fAL, --filterAbscissaLower:
      filter on Abscissa (x-axis) data greater or equal than value specified.
  -fAU, --filterAbscissaUpper:
      filter on Abscissa (x-axis) data lesser or equal than value specified.
  -s, --chunkSize:
      Split the input file into chunks of this size. Default: not set, no chunking
  -w, --window_offsets:
      Specify windowing offset for evolutive analysis. Default: chunkSize if set, or None.
  -z, --columnSeparator:
      Use this as columns separtor for input file instead of default tab \t. Default: \t
  -k, --skipHeaderLines:
      Specify number of (header)lines to be skipped. Default: 0
  -n, --decimate:
      Specify decimation of input data. Default: 1 (take every value)
  -a, --abscissaColumn:
      Specify column number for abscissa (time or depth). Default: 1
  -o, --ordinateColumn:
      Specify column number for ordinate (data). Default: 2
  -C, --complex:
      Compute spectra for complex input. Default: false. If true, by default takes values from abscissa column + 1.
  -Co, --ComplexOrdinateColumn:
      Specify column number for ordinate (data), complex part. Default: 3 (if C option used)
  -t, --deltaT:
      Specify time/depth offset between data. Default: 1.0
  -x, --ratios:
      Compute and output frequency ratio matrix. Default: false.
  --ratio_list:
      If ratio_option is set, provides comma separated list of frequency ratios to highlight in output.
  --output_detrended_data:
      Output detrended data to this filename. Not implemented yet.
  --version:
      Shows current version and build of this software
  -v, --verbose:
      Print verbose messages. Specify multiple times to increase verbosity.
  --references:
      Shows references to papers and sources.
  --licenses:
      Shows licenses of code components used in this software.
```

## Example: Try it on a sample file: etp_1kyr_36Ma.txt
```
user$ ./hpnaff etp_1kyr_36Ma.txt -d 1 -m 75
```


```
===================================================================================================================================================
    # | chunk |       t0       | points |    frequency |     amplitude |        phase | significance |         period |    arcsec/yr |  phase (deg)
    0 |     0 |         0.0000 |  14001 |  0.000000000 |    0.00000000 |  0.000000000 |  0.000000000 |            inf |    0.0000000 |    0.0000000
    1 |     0 |         0.0000 |  14001 |  0.002472456 |    0.75721258 |  3.740563127 |  0.795224850 |   404.45610756 |    3.2043032 |  214.3184802
    2 |     0 |         0.0000 |  14001 |  0.010522542 |    0.59639514 | -0.083030087 |  0.840326031 |    95.03406718 |   13.6372149 |   -4.7572736
    3 |     0 |         0.0000 |  14001 |  0.024999241 |    0.53877772 |  5.657761056 |  0.827374133 |    40.00121445 |   32.3990163 |  324.1658300
    4 |     0 |         0.0000 |  14001 |  0.008040668 |    0.50803290 |  2.884638325 |  0.832167805 |   124.36778406 |   10.4207051 |  165.2776014
    5 |     0 |         0.0000 |  14001 |  0.010102115 |    0.39766702 |  4.191460276 |  0.868275179 |    98.98917382 |   13.0923408 |  240.1529838
    6 |     0 |         0.0000 |  14001 |  0.042828239 |    0.31115152 |  1.466354530 |  0.901608023 |    23.34908056 |   55.5053976 |   84.0159259
    7 |     0 |         0.0000 |  14001 |  0.007627856 |    0.30706911 |  0.558756463 |  0.900168886 |   131.09844523 |    9.8857008 |   32.0143871
    8 |     0 |         0.0000 |  14001 |  0.000407279 |    0.27097231 |  2.625060814 |  0.917986628 |  2455.31930796 |    0.5278336 |  150.4049056
    9 |     0 |         0.0000 |  14001 |  0.045302069 |    0.26148887 |  5.142663822 |  0.908826103 |    22.07404698 |   58.7114815 |  294.6529325
   10 |     0 |         0.0000 |  14001 |  0.025832617 |    0.23907704 |  3.749692849 |  0.917331368 |    38.71075060 |   33.4790718 |  214.8415747
   11 |     0 |         0.0000 |  14001 |  0.001034652 |    0.22887920 |  0.316263467 |  0.926607014 |   966.50814962 |    1.3409095 |   18.1205619
   12 |     0 |         0.0000 |  14001 |  0.053350844 |    0.23311864 |  1.355155680 |  0.912974999 |    18.74384580 |   69.1426943 |   77.6447010
   13 |     0 |         0.0000 |  14001 |  0.009911387 |    0.18079485 |  5.692179885 |  0.941223401 |   100.89404923 |   12.8451580 |  326.1378837
   14 |     0 |         0.0000 |  14001 |  0.019210470 |    0.15673509 | -0.475810951 |  0.949106113 |    52.05494759 |   24.8967689 |  -27.2619593
   15 |     0 |         0.0000 |  14001 |  0.052929481 |    0.14841497 |  5.755575872 |  0.947991689 |    18.89306257 |   68.5966076 |  329.7702061
   16 |     0 |         0.0000 |  14001 |  0.009488271 |    0.15112041 | -0.166596764 |  0.960411285 |   105.39327551 |   12.2967997 |   -9.5452915
   17 |     0 |         0.0000 |  14001 |  0.025098798 |    0.17311913 |  2.629896209 |  0.958954086 |    39.84254538 |   32.5280423 |  150.6819533
   18 |     0 |         0.0000 |  14001 |  0.025414340 |    0.14355863 |  4.818018103 |  0.948646430 |    39.34786366 |   32.9369851 |  276.0521029
   19 |     0 |         0.0000 |  14001 |  0.024888808 |    0.14098312 |  5.710170222 |  0.964800625 |    40.17870229 |   32.2558949 |  327.1686540
   20 |     0 |         0.0000 |  14001 |  0.010323185 |    0.13699989 | -1.010239654 |  0.950068676 |    96.86932826 |   13.3788478 |  -57.8824685
   21 |     0 |         0.0000 |  14001 |  0.007440236 |    0.12986997 |  2.040413354 |  0.952028440 |   134.40433245 |    9.6425463 |  116.9070736
   22 |     0 |         0.0000 |  14001 |  0.001455303 |    0.12708837 |  2.315320953 |  0.947284449 |   687.14223688 |    1.8860724 |  132.6581188
   23 |     0 |         0.0000 |  14001 |  0.010925100 |    0.11817448 | -0.359553604 |  0.955088752 |    91.53234018 |   14.1589300 |  -20.6009040
   24 |     0 |         0.0000 |  14001 |  0.009667251 |    0.11874251 |  2.622499903 |  0.954076811 |   103.44202802 |   12.5287567 |  150.2581762
   25 |     0 |         0.0000 |  14001 |  0.010413463 |    0.09715315 |  2.916722459 |  0.980490597 |    96.02953115 |   13.4958485 |  167.1158869
   26 |     0 |         0.0000 |  14001 |  0.000901134 |    0.09569726 |  1.601650029 |  0.968540900 |  1109.71245232 |    1.1678701 |   91.7677869
   27 |     0 |         0.0000 |  14001 |  0.012994167 |    0.10086003 |  0.611296497 |  0.967944111 |    76.95760564 |   16.8404408 |   35.0247093
   28 |     0 |         0.0000 |  14001 |  0.042724400 |    0.09486676 |  1.049712893 |  0.980633322 |    23.40582906 |   55.3708222 |   60.1441185
   29 |     0 |         0.0000 |  14001 |  0.000231340 |    0.09813349 |  0.120421080 |  0.967044518 |  4322.64694648 |    0.2998163 |    6.8996196
   30 |     0 |         0.0000 |  14001 |  0.005563540 |    0.10009924 |  2.397732816 |  0.967425007 |   179.74168484 |    7.2103475 |  137.3799708
   31 |     0 |         0.0000 |  14001 |  0.010638958 |    0.09688648 | -0.508419978 |  0.964555303 |    93.99417100 |   13.7880891 |  -29.1303190
   32 |     0 |         0.0000 |  14001 |  0.008918727 |    0.08741944 |  2.593049312 |  0.963917464 |   112.12362149 |   11.5586705 |  148.5707817
   33 |     0 |         0.0000 |  14001 |  0.007853668 |    0.09431921 |  1.158980960 |  0.966015093 |   127.32903632 |   10.1783540 |   66.4047175
   34 |     0 |         0.0000 |  14001 |  0.024592220 |    0.08828171 |  3.104291584 |  0.962948350 |    40.66326600 |   31.8715176 |  177.8628061
   35 |     0 |         0.0000 |  14001 |  0.008385587 |    0.08641115 |  2.992612020 |  0.964466604 |   119.25223732 |   10.8677206 |  171.4640385
   36 |     0 |         0.0000 |  14001 |  0.018145511 |    0.08763660 |  3.810401603 |  0.965100179 |    55.11004768 |   23.5165828 |  218.3199301
   37 |     0 |         0.0000 |  14001 |  0.024810831 |    0.08240291 |  4.521072773 |  0.982190253 |    40.30497880 |   32.1548364 |  259.0383888
   38 |     0 |         0.0000 |  14001 |  0.042929640 |    0.09470460 | -1.714770886 |  0.974168391 |    23.29392958 |   55.6368128 |  -98.2491346
   39 |     0 |         0.0000 |  14001 |  0.001312871 |    0.08074684 |  0.777678607 |  0.975863805 |   761.68923898 |    1.7014813 |   44.5577020
   40 |     0 |         0.0000 |  14001 |  0.002884917 |    0.08065197 |  2.947032331 |  0.969528197 |   346.63039200 |    3.7388528 |  168.8525146
   41 |     0 |         0.0000 |  14001 |  0.009282092 |    0.07995912 | -0.902914836 |  0.959472393 |   107.73432964 |   12.0295917 |  -51.7332094
   42 |     0 |         0.0000 |  14001 |  0.002172887 |    0.07898195 |  1.985388299 |  0.966843584 |   460.21714558 |    2.8160620 |  113.7543702
   43 |     0 |         0.0000 |  14001 |  0.025182288 |    0.07351399 |  1.270869502 |  0.989395181 |    39.71045117 |   32.6362447 |   72.8154588
   44 |     0 |         0.0000 |  14001 |  0.002006049 |    0.07122474 |  1.007509645 |  0.970121706 |   498.49230512 |    2.5998395 |   57.7260505
   45 |     0 |         0.0000 |  14001 |  0.008233354 |    0.07018462 |  4.902075131 |  0.962605563 |   121.45718541 |   10.6704267 |  280.8682159
   46 |     0 |         0.0000 |  14001 |  0.043838032 |    0.06715352 |  3.232504867 |  0.965810782 |    22.81124276 |   56.8140900 |  185.2088861
   47 |     0 |         0.0000 |  14001 |  0.007224019 |    0.06531452 |  4.276908703 |  0.978785703 |   138.42710917 |    9.3623280 |  245.0488180
   48 |     0 |         0.0000 |  14001 |  0.045402169 |    0.07436495 |  2.063108445 |  0.977076727 |    22.02537938 |   58.8412112 |  118.2074066
   49 |     0 |         0.0000 |  14001 |  0.002344771 |    0.06402237 |  1.101819579 |  0.974716549 |   426.48086343 |    3.0388233 |   63.1296117
   50 |     0 |         0.0000 |  14001 |  0.034012471 |    0.06632261 |  0.311286623 |  0.964879805 |    29.40098080 |   44.0801621 |   17.8354097
   51 |     0 |         0.0000 |  14001 |  0.000615717 |    0.06955959 |  0.324652566 |  0.966238480 |  1624.12260860 |    0.7979693 |   18.6012218
   52 |     0 |         0.0000 |  14001 |  0.024999132 |    0.06125458 |  3.871230190 |  0.987944595 |    40.00138890 |   32.3988750 |  221.8051514
   53 |     0 |         0.0000 |  14001 |  0.012584554 |    0.06352103 |  4.354136361 |  0.970616108 |    79.46248859 |   16.3095823 |  249.4736369
   54 |     0 |         0.0000 |  14001 |  0.007042302 |    0.05954369 |  3.178070923 |  0.967975779 |   141.99901671 |    9.1268238 |  182.0900509
   55 |     0 |         0.0000 |  14001 |  0.053168977 |    0.06126971 | -0.361831203 |  0.974607869 |    18.80796025 |   68.9069938 |  -20.7314008
   56 |     0 |         0.0000 |  14001 |  0.018476857 |    0.05675901 |  1.358226583 |  0.977498040 |    54.12175935 |   23.9460065 |   77.8206508
   57 |     0 |         0.0000 |  14001 |  0.025289765 |    0.06621217 | -0.908124254 |  0.973848255 |    39.54168736 |   32.7755361 |  -52.0316871
   58 |     0 |         0.0000 |  14001 |  0.045202487 |    0.05979910 |  4.539755712 |  0.986931113 |    22.12267664 |   58.5824230 |  260.1088423
   59 |     0 |         0.0000 |  14001 |  0.001552070 |    0.05459046 |  3.384103161 |  0.981459103 |   644.30080647 |    2.0114828 |  193.8948286
   60 |     0 |         0.0000 |  14001 |  0.010774958 |    0.06013873 | -0.920994981 |  0.970525666 |    92.80778711 |   13.9643455 |  -52.7691253
   61 |     0 |         0.0000 |  14001 |  0.035180014 |    0.05381168 |  6.108664303 |  0.971837152 |    28.42523015 |   45.5932984 |  350.0006830
   62 |     0 |         0.0000 |  14001 |  0.025932854 |    0.05497186 |  0.528441367 |  0.981080756 |    38.56112370 |   33.6089791 |   30.2774600
   63 |     0 |         0.0000 |  14001 |  0.019098347 |    0.05289607 |  5.666427635 |  0.977106362 |    52.36055271 |   24.7514576 |  324.6623884
   64 |     0 |         0.0000 |  14001 |  0.052829237 |    0.05592637 |  4.622260192 |  0.978276065 |    18.92891225 |   68.4666917 |  264.8360008
   65 |     0 |         0.0000 |  14001 |  0.011547573 |    0.05030718 |  3.661455364 |  0.969207741 |    86.59828253 |   14.9656548 |  209.7859392
   66 |     0 |         0.0000 |  14001 |  0.052667555 |    0.04771644 |  4.565272429 |  0.975236135 |    18.98702155 |   68.2571512 |  261.5708425
   67 |     0 |         0.0000 |  14001 |  0.008514638 |    0.04804064 |  4.938804770 |  0.978348985 |   117.44480713 |   11.0349707 |  282.9726692
   68 |     0 |         0.0000 |  14001 |  0.017723966 |    0.04660213 |  1.785804511 |  0.976900743 |    56.42078067 |   22.9702600 |  102.3190615
   69 |     0 |         0.0000 |  14001 |  0.013408861 |    0.04526629 | -0.218064195 |  0.976731996 |    74.57754915 |   17.3778840 |  -12.4941580
   70 |     0 |         0.0000 |  14001 |  0.004944881 |    0.04401069 |  3.792380702 |  0.977679719 |   202.22932980 |    6.4085660 |  217.2874085
   71 |     0 |         0.0000 |  14001 |  0.024477668 |    0.04221972 |  2.795577253 |  0.982047896 |    40.85356574 |   31.7230571 |  160.1747779
   72 |     0 |         0.0000 |  14001 |  0.025733760 |    0.04584982 |  3.159868039 |  0.992264403 |    38.85945885 |   33.3509534 |  181.0471025
   73 |     0 |         0.0000 |  14001 |  0.044937430 |    0.04272239 |  0.214967789 |  0.979100871 |    22.25316390 |   58.2389096 |   12.3167471
   74 |     0 |         0.0000 |  14001 |  0.018373301 |    0.04170326 |  4.267372537 |  0.976418419 |    54.42680232 |   23.8117976 |  244.5024360
   75 |     0 |         0.0000 |  14001 |  0.016041231 |    0.04508949 |  2.028721000 |  0.980103195 |    62.33935561 |   20.7894353 |  116.2371511
   ```

## Notes
"Evolutive" outputs can be obtained with the -s option, e.g. -s 1000 computes an analysis for every 1000 points.

Source code available at repository below.
Note: There are three branches, "main" contains a new swift package manager version, which can also be compiled for Linux/Windows, and "Accelerate" whichh uses the LinearAlgebra system library. The Linux version is linked statically. The "Windows" branch is modified to enable compilation with the Windows Toolchain [^1]

https://paloz.marum.de/bitbucket/scm/essp/hpnaff.git


Notes: 
This is a Swift implementation of 
https://ops.aps.anl.gov/manuals/SDDStoolkit/SDDStoolkitsu61.html
Related efforts:
https://github.com/MichaelEhrlichman/FortNAFF
https://github.com/adrn/SuperFreq
https://github.com/nkarast/PyNAFF

## References
References:
1. Laskar, J., 1990, The chaotic motion of the Solar System. A numerical estimate of the size of the chaotic zones, Icarus, 88, 266-291.
2. Laskar, J., 1993, Frequency analysis for multi-dimensional systems. Global dynamics and diffusion, Physica D, 67, 257-281.
3. Dumas, S., Laskar, J., 1993, Global Dynamics and Long-Time Stability in Hamiltonian Systems Via Numerical Frequency Analysis, Phys. Rev. Letters, 70 (20), 2975-2979.
4. Laskar, J. : 1999, Introduction to frequency map analysis, in proc. of NATO ASI 533 3DHAM95, S'Agaro, Spain, 134150.
5. Papaphilippou, Y., Frequency maps for LHC models, PAC99.
6. Papahilippou, Y. Zimmermann, F., Weak-strong beam-beam simulations for the Large Hadron Collider, Phys. Rev. ST Accel.
Beams 2, 104001 (1999).
7. Robin, D., Steir, C., Laskar, J., Nadolski, L. : 2000, Global dynamics of the ALS revealed through experimental Frequency Map Analysis, Phys. Rev. Let., 85, pp. 558-561.
8. Laskar, J., Frequency map analysis and quasiperiodic decompositions. preprint (https://arxiv.org/pdf/math/0305364.pdf) (2003).
9. [Valluri & Merritt, 1998](http://iopscience.iop.org/article/10.1086/306269/fulltext/37764.text.html#sc2)
