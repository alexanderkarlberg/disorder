[![Build Status](https://img.shields.io/github/actions/workflow/status/alexanderkarlberg/disorder/cmake-single-platform.yml?label=build&logo=github&style=flat-square)](https://github.com/alexanderkarlberg/disorder/actions/workflows/cmake-single-platform.yml)

Welcome to `disorder` - a program to compute fixed order predictions
at high order in deep inelastic scattering! Below are some
installation and run instructions. More detailed help can be found in
the `docs` folder.

Although the code is released under GPLv3 or later any scientific use of the code should result in a citation of

Karlberg, A., disorder: Deep inelastic scattering at high orders. [SciPost Phys. Codebases 32 (2024)](https://doi.org/10.21468/scipostphyscodeb.32), [arXiv:2401.16964](https://arxiv.org/abs/2401.16964).

More detailed information about the use of the code can be found in the docs folder.

Installation
============

To compile `disorder`, you will need

* [<span style="font-variant:small-caps;">Hoppet</span>](https://github.com/hoppet-code/hoppet), v1.3.0 or newer [currently also needs the 2024-01-n3lo-splittings-functions branch]
* [LHAPDF6](http://lhapdf.hepforge.org/) [tested with v6.5.4]

Optionally you may need [fastjet](https://fastjet.fr/) installed as well.

Once all dependencies are installed on your machine, disorder can be
compiled using `cmake`:

	mkdir build && cd build
  	cmake ..

  	make [-j]

from the main directory. This will create an executable `disorder` along
with two auxiliary executables, `mergedata` and `getpdfuncert`.

`disorder` can also be installed in the user's path by invoking

	make install

This will install `disorder` in the default location (/usr/local/bin
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

* The [POWHEG-BOX](https://powhegbox.mib.infn.it/) under GPLv2. Specifically the analysis framework and the `mergedata` programs are adapted from there.
* The command line tools (io_utils.f90 and lcl_dec.f90) are written by Gavin Salam and are under GPLv3. 
* Some of the code is adapted from [proVBFH](https://github.com/fdreyer/proVBFH/) under GPLv3.
* The code relies heavily on disent, written by Mike Seymour. The version included here is based on the one included with v1.0.5 of [dispatch](https://github.com/gavinsalam/dispatch), but has received significant additional modifications.

Citation policy
===============

It is difficult to give an exhaustive list of references for all use-cases of `disorder`. The following should therefore be seen as a minimal set:

Any use of the code should always result in a citation of

* Karlberg, A. (2024). disorder: Deep inelastic scattering at high orders. SciPost Physics Codebases. https://doi.org/10.21468/scipostphyscodeb.32, [arXiv:2401.16964](https://arxiv.org/abs/2401.16964).

<span style="font-variant:small-caps;">Hoppet</span> should also always be cited

* G.P. Salam, J. Rojo, A Higher Order Perturbative Parton Evolution Toolkit (HOPPET), [arXiv:0804.3755](https://arxiv.org/abs/0804.3755).

along with references for the DIS coefficient functions at the appropriate order (see Refs. [3]-[15] in the paper).  

Whenever the code is run in the P2B mode the following paper should also be cited

* S. Catani, M.H. Seymour, A General algorithm for calculating jet cross-sections in NLO QCD, [arXiv:hep-ph/9605323](https://arxiv.org/abs/hep-ph/9605323).

Bugs
====

There are no known bugs at the moment.

In the past there have been bugs, discovered by myself and others. I would in particular like to thank

* Melissa van Beekveld
* Silvia Ferrario Ravasio

for reporting a number of them!
