----------------------------------------------------------------------
 Installation
---

To compile disorder, you will need
* hoppet, v 1.3.0 or newer (https://github.com/gavinsalam/hoppet)
* LHAPDF (http://lhapdf.hepforge.org/)

Once all dependencies are installed on your machine, disorder can be
compiled using:

  ./configure [--with-LHAPDF=DIR --with-hoppet=DIR]

  make [-j]

in the main directory. This will create an executable "disorder".

No arguments need to be passed to configure if lhapdf-config and
hoppet-config are both in the user's $PATH.

To compile disorder with a local installation of hoppet or LHAPDF,
change the "DIR" above in the configure step to the folder containing
the hoppet-config and lhapdf-config executables respectively.

Two auxiliary executable, mergedata and getpdfuncert, can be compiled
with

  make aux [-j]

----------------------------------------------------------------------
Usage
---

To run disorder, use the disorder executable, and pass command
line arguments to specify inputs. Example:

  ./disorder -pdf MSTW2008nlo68cl -nlo -Q 12.0 -x 0.1


To get a list of parameters which can be specified run

   ./disorder -help
   
The program will create a file "xsct_[...].dat" containing the total
cross section and Monte Carlo error in addition to any scale or PDF
uncertainties included.

You can include your favourite analysis in the analysis directory and
modify the Makefile accordingly. The output of the analysis is stored
in 'disorder_[...].dat' files.

The full list of possible options can be obtained from
src/mod_parameters.f90.

For a more detailed usage description please look in the manual which
can be found in the docs directory.


----------------------------------------------------------------------
