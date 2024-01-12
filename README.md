# CDFvidar (Regridding initial conditions for CCAM)

CDFvidar is used to prepare initial conditions and mesonest atmospheric data
for downscaling with the Conformal Cubic Atmospheric Model (CCAM).  Typically
CDFvidar is used to interpolate analyses, reanalyses and General Circulation
Models to the conformal cubic grid used by CCAM.

## Website

For documentation, see our website at

[https://confluence.csiro.au/display/CCAM/CCAM]

## Dependencies

CDFvidar requires the NetCDF C library.

## Building CDFvidar

To build CDFvidar with intel, gnu and cray fortran compiler use

```
make
make GFORTRAN=yes
make CRAY=yes
```

Debugging is also enabled with

```
make TEST=yes
```
