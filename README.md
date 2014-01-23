MESH2MCNP
=========

Mesh conversion tool using the PENMSHXP/PENTRAN geometry format

MESH2MCNP matches discrete geometry 

Keeping your local MESH2MCNP up-to-date
---------------------------------------

```
cd mesh2mcnp/
git pull # this will work if you have not migrated anything
```

If `git pull` does not work it is because there are conflicting files, if you do not care for this, just issue `git stash` prior to `git pull`, and that will tuck away the problem files (you can recover stashed files too).

When weird things happen, just delete the `mesh2mcnp` folder, and re-clone with `git clone git@github.com:critcel/mesh2mcnp`

MESH2MCNP should also have a home at `/path/to/software/mesh2mcnp` and system paths will be updated to find this location for using MESH2MCNP.

Example DRCC Problem
--------------------

To run DRCC example, and generate a `.mc` file
```
cd case3-drcc
python ../mesh2mcnp.py -f drcc_out.pen -mcm drcc.mcm -nofm -lp
```

Dissecting the options, run `python mesh2mcnp.py` and you'll get the arguments.

- `f`: Supply the `.pen` file (MANDATORY)
- `mcm`: Give the GMIX mcm file and it will insert the mcm file in place of the multigroup based material cards.  Also, this option replaces the densities in the cell cards accordingly.
- `nofm`: Turn off the FMESH tallies.
- `lp`: Print the log to stdout; if you don't use this m2mc.log is created and you can look there.
