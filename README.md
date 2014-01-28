MESH2MCNP
=========

Mesh conversion tool using the PENMSHXP/PENTRAN geometry format

MESH2MCNP matches discrete geometry so that you cannot blame the Sn meshing in your Monte Carlo analysis.

Keeping your local MESH2MCNP up-to-date
---------------------------------------

```
cd mesh2mcnp/
git pull # this will work if you have not made changes to targeted files in the local clone
```

If `git pull` does not work it is because there are conflicting files, if you do not care for this, just issue `git stash` prior to `git pull`, and that will tuck away the problem files (you can recover stashed files too).

When weird things happen, just delete the `mesh2mcnp` folder, and re-clone with `git clone https://github.com/critcel/mesh2mcnp.git mesh2mcnp`

MESH2MCNP should also have a home at `/path/to/software/mesh2mcnp` and system paths will be updated to find this location for using MESH2MCNP.

Help
----
Type `mesh2mcnp` (on HPC clusters) to see the options.
Dissecting the options, run `python mesh2mcnp.py` and you'll get the options listed.

- `f`: Supply the `.pen` file (MANDATORY)
- `mcm`: Give the GMIX mcm file and it will insert the mcm file in place of the multigroup based material cards.  Also, this option replaces the densities in the cell cards accordingly.
- `nofm`: Turn off the FMESH tallies.
- `lp`: Print the log to stdout; if you don't use this m2mc.log is created and you can look there.

Other options not explored by the examples:

- `f4mul`: tally multiplier for f4 universe mesh tallies (`-f4mul 112.5`)
- `utmat`: filter tallies against specific materials of interest (e.g. only materials 8,2,11 `-utmat 8 2 19`)
- `unit`: assign unit volumes to material universe cells

Example DRCC Problem
--------------------

### Continuous Option
To run DRCC example, and generate a `.mc` file (continuous)
```
cd case3-drcc
python ../mesh2mcnp.py -f drcc_out.pen -mcm drcc.mcm -nofm -lp
cp drcc_out.mc drcc_out_continuous.mc
```
### Multigroup Option
To run DRCC example, and generate a `.mc` file (xsmcnp/multigroup)
```
cd case3-drcc
python ../mesh2mcnp.py -f drcc_out.pen -mcm drcc.mcm -nofm -lp
cp drcc_out.mc drcc_out_multigroup.mc
```

