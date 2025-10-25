# AMM75
AMM75 is a NEMO ocean circulation configuration. It is a physics only shelf sea model of the Northwest European Shelf, that is based off of the AMM15 configuration (Graham et al 2018).
AMM75 specifically is the same as AMM15 but coarsened to 7.5 km horizontal resolution. This has been created to investigate sensitivities to model resolution in a shelf sea model. The AMM15 reference repository is [https://github.com/NOC-MSM/CO_AMM15_CHAMFER](https://github.com/NOC-MSM/CO_AMM15_CHAMFER/tree/v0.1.0). This configuration is based on a version of NEMO (v4.0.4) that has functioning momentum trend diagnostics
```
svn -r 15194 co https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/branches/UKMO/NEMO_4.0.4_momentum_trends nemo_4.0.4_trd).
```
The NEMO directiory of this repository includes the AMM15/AMM75 specific source code (SRC), namelists (NAMELISTS), boundary coordinate files, input/output files (XIOS) as well as various High Performance Computing files which are specific to ARCHER2. The domain file is available from https://gws-access.jasmin.ac.uk/public/jmmp/AMM7/AMM75/INPUTS/domain_cfg.nc

Forcing datasets are those used to force the AMM15 configuration, specifically:
The model has been forced using ERA5 for the surface forcing, GloASea6 for the lateral boundaries, FES2014 for the tides and a River climatology [https://gws-access.jasmin.ac.uk/public/jmmp/AMM7/AMM75/RIV/](https://gws-access.jasmin.ac.uk/public/jmmp/AMM7/AMM75/RIV/AMM75_River_Climatology.nc) . Model simulations are initialised using GloSea6 temperature and salinity.


Refs.
Graham, J. A., O’Dea, E., Holt, J., Polton, J., Hewitt, H. T., Furner, R., . . . others
(2018). Amm15: a new high-resolution nemo configuration for operational
simulation of the european north-west shelf. Geoscientific Model Development,
11 (2), 681–696. 
