# bespoke_amber_restraints

command-line programs for using X-ray crystallography data to optimize harmonic restraints to the bare minimum required to keep the structure from leaving the electron density

sorry the documentation is currently poor.  Start with the master script optimize_weights_runme.com.  It expects to have certain files available.
example data here:
https://bl831.als.lbl.gov/~jamesh/amber/1aho/example_starter4.tgz 

if you want to use the automatic hydration features you will need to compile the two short C programs included in this repo:

### [float_add](docs/float_add.md)

> Add, subtract, scale and offset raw floating-point flat files with arbitrary headers.

### [float_func](docs/float_func.md)

> Perform any C function on one or two raw floating-point flat files with arbitrary headers.



## Author
<ADDRESS><A HREF="mailto:JMHolton@lbl.gov">James Holton <JMHolton@lbl.gov> </A></ADDRESS>
