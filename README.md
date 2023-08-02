# disorder
----------------------------------------------------------------------
Installation
---

To compile disorder, you will need
* hoppet, the struct-func-devel branch (http://hoppet.hepforge.org/)
* LHAPDF (http://lhapdf.hepforge.org/)

The struct-func-devel branch of hoppet can be downloaded using

  svn checkout https://svn.hepforge.org/hoppetsvn/branches/struct-func-devel/

Once all dependencies are installed on your machine, disorder
can be compiled using:

  ./configure [--with-LHAPDF=DIR --with-hoppet=DIR]
  make

in the main directory. This will create an executable "disorder".

No arguments need to bne passed to configure if lhapdf=config and
hoppet-config are both in the user PATH.

To compile disorder with a local installation of hoppet or
LHAPDF, change the "DIR" above in the configure step to the folder
containing the hoppet-config and lhapdf-config executables
respectively.


----------------------------------------------------------------------
Usage
---

To run disorder, use the disorder executable, and pass command
line arguments to specify inputs. Example:

  ./disorder -pdf MSTW2008nlo68cl -nlo -Q 12.0 -x 0.1

The program will create a file "xsct_nlo.dat" containing the total
cross section and Monte Carlo error in addition to standard
hiostograms of Q, x and y.

A list of possible options can be obtained from parameters.f90.


----------------------------------------------------------------------
