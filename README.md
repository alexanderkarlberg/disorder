[![Build Status](https://img.shields.io/github/actions/workflow/status/alexanderkarlberg/disorder/cmake-single-platform.yml?label=build&logo=github&style=flat-square)](https://github.com/alexanderkarlberg/disorder/actions/workflows/cmake-single-platform.yml)

Welcome to `disorder` - a program to compute fixed order predictions
at high order in deep inelastic scattering! Below are some
installation and run instructions. More detailed help can be found in
the `docs` folder.

Although the code is released under GPLv3 or later any scientific use of the code should result in a citation of

Alexander Karlberg, arXiv:2401.XXXXX.

More detailed information about the use of the code can be found in the docs folder.

Installation
============

To compile `disorder`, you will need

* <span style="font-variant:small-caps;">Hoppet</span>, v1.3.0 or newer (https://github.com/hoppet-code/hoppet) [currently also needs the 2024-01-n3lo-splittings-functions branch]
* LHAPDF6 (http://lhapdf.hepforge.org/) [tested with v6.5.4]

Optionally you may need fastjet installed as well
(https://fastjet.fr/).

Once all dependencies are installed on your machine, disorder can be
compiled using `cmake`:

	mkdir build && cd build
  	cmake ..

  	make [-j]

from the main directory. This will create an executable `disorder` along
with two auxiliary executables, `mergedata` and `getpdfuncert`.

`disorder` can also be installed in the user's path by invoking

	make install

This will install disorder in the default location (/usr/local/bin
typically). The user can change the path by specifying

	cmake -DCMAKE_INSTALL_PREFIX=/path/to/subdir ..

If `hoppet-config` or `lhapdf-config` are not in the user's path, the full
path can be specified manually through

	cmake -DHOPPET_CONFIG=/path/to/hoppet-config -DLHAPDF_CONFIG=/path/to/lhapdf-config ..

where the path should include the config itself
(i.e. `/usr/local/bin/hoppet-config`)

By default fastjet is not linked and only a skeleton analysis
(`analysis/simple_analysis.f`) is compiled. To link fastjet run

	cmake -DNEEDS_FASTJET=ON [-DFASTJET_CONFIG=/path/to/fastjet-config] ..

where the path to `fastjet-config` only needs to be specified if it is
not in the user's `$PATH`.

To compile a different analysis the user should first put it in the
analysis directory (here we assume it to be called `my_analysis.f`),
and then pass it to `cmake` through

	cmake -DANALYSIS=my_analysis.f ..

Usage
=====

To run `disorder`, use the disorder executable, and pass command
line arguments to specify inputs. Example:

	./disorder -pdf MSTW2008nlo68cl -n3lo -Q 12.0 -x 0.01 -scaleuncert


To get a list of parameters which can be specified run

	./disorder -help
   
The program will create a file `xsct_[...].dat` containing the total
cross section and Monte Carlo error in addition to any scale or PDF
uncertainties included.

You can include your favourite analysis in the analysis directory and
pass it to `cmake` as outlined above. The output of the analysis is stored
in `disorder_[...].dat` files.

The full list of possible options can be obtained from
`src/mod_parameters.f90`.

For a more detailed usage description please look in the manual which
can be found in the `docs` directory.

Third party code
================

Besides the dependencies listed above `disorder` incorporates code from the following sources:

* The POWHEG-BOX (https://powhegbox.mib.infn.it/) under GPLv2. Specifically the analysis framework and the `mergedata` programs are adapted from there.
* The command line tools (io_utils.f90 and lcl_dec.f90) are written by Gavin Salam and are under GPLv3. 
* Some of the code is adapted from proVBFH (https://github.com/fdreyer/proVBFH/) under GPLv3.
* The code relies heavily on disent, written by Mike Seymour. 

Bugs
====

There are no known bugs.
