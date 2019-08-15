# FUSE

## Description and credits

This is a source code repository for the **Framework for Understanding Structural Errors** or **FUSE**. The initial implementation is described in [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735). This implementation involves four main additional features:

1. a snow module described in [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736),
2. a calibration mode using [Duan et al. (WRR, 1992)](http://dx.doi.org/10.1029/91WR02985) shuffled complex evolution alogrithm (SCE),
3. a distributed mode enabling to run FUSE on a grid with spatially distributed parameter values, and
4. all the input, output and parameter files are now NetCDF files.

## FUSE modes and case studies

FUSE can be run in five complementary modes:

1. `run_def` runs FUSE with default parameter values,
2. `calib_sce` runs FUSE in an SCE-calibration mode,
3. `run_best` runs FUSE using the best (lowest RMSE) parameter set found by SCE.
4. `run_pre_catch` runs FUSE with parameter sets from previously calibrated catchment(s),
5. `run_pre_dist` runs FUSE with spatially distributed parameter values.

To get you started with FUSE, we provide files for two case studies involing modeling at different spatial scales:

* catchment scale: forcing and streamflow data for the [USGS 09066300 MIDDLE CREEK NEAR MINTURN, CO.](https://waterdata.usgs.gov/nwis/inventory/?site_no=09066300&agency_cd=USGS&) catchment - available [here [0.5MB]](
https://dl.dropboxusercontent.com/s/f6omcgz8hsirlr0/fuse_catch.zip?dl=0) for download,  
* grid scale: forcing from a climate model on a 1/8th degree grid for a 58 x 28 grid cells domain - available [here [32MB]](
https://dl.dropboxusercontent.com/s/vt1vcyqkld8zqyc/fuse_grid.zip?dl=0) for download.

We recommend that you follow the steps A to F below to run FUSE for the catchment case study, and then follow step G to run FUSE on a grid using spatially distributed parameter values.

## A. Fork this repository and compile FUSE
1. Fork this repository to a directory `$(MASTER)` on your machine (see the [SUMMA manual](http://summa.readthedocs.io/en/latest/development/SUMMA_and_git/) for a step-by-step procedure)
1. Change directory to `$(MASTER)/build/` and edit the `Makefile`, by:
   1. defining the name of the master directory (line 10),
   2. defining the fortran compiler (line 202),
   3. defining the path to the NetCDF libraries (lines 204-225, note that the NetCDF libraries should be compiled using the same compiler you are using to compile FUSE).
 1. Compile the SCE code (see Section H below).
 1. Compile the FUSE code (type `make`).
 1. Change to `$(MASTER)/bin/` and try running FUSE by typing `./fuse.exe`. If the output is `1st command-line argument is missing (fileManager)`, you have probably compiled FUSE correctly. 
 
## B. Populate the bin directory
FUSE relies on a `FILEMANAGER`, which defines the paths to the FUSE `settings`, `input`, `output` directories, as well as other settings essential to run FUSE (described in Sections C and D). 
1. Move the file `fm_catch.txt` provided for the catchment case study to `$(MASTER)/bin/`. This is the `FILEMANAGER` for the catchment case study.
1. Update the lines 4 to 6 of `fm_catch.txt` using the path of the `fuse_catch` directory on your machine.
   
## C. Populate the setup directory
The `setup` directory must contain the following files (provided for the catchment case study):

   1. The file `M_DECISIONS` (called `fuse_zDecisions_902.txt` in the case studies) describes the different options available in the FUSE modeling framework. These modeling decisions are described in detail by [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735), except decision 9 described in [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736).
   2. The file `CONSTRAINTS` (called `fuse_zConstraints_snow.txt` in the case studies) defines in particular the default parameter values and lower and upper parameter bounds. The list of parameters corresponds to those described in [Clark et al. (WRR, 2008)](http://dx.doi.org/10.1029/2007WR006735) and [Henn et al. (WRR, 2015)](http://dx.doi.org/10.1002/2014WR016736). 
   3. The file `MOD_NUMERIX` (called `fuse_zNumerix.txt` in the case studies) defines decisions regarding the numerical solution technique. Examples of the impact of these decisions are described by [Clark and Kavetski (WRR 2010)](http://dx.doi.org/10.1029/2009WR008894) and [Kavetski and Clark (WRR 2010)](http://dx.doi.org/10.1029/2009WR008896).
   4. The file `FORCINGINFO` (called `input_info.txt` in the case studies) provides metadata for the NetCDF input file. It defines the name and units of the variables in the input file.

## D. Populate the input directory
The `input` directory must contain the following files (provided for the catchment case study):

   1. The file `forcefile` (called `us_09066300_input.nc` in the catchment case study) contains the input data in a 2D (resp. 3D) arrays for modeling at the catchment (resp. grid) scale. The name of this file is made by appending the `suffix_forcing` defined in the `FILEMANAGER` (see B) to the basin ID (see E `$2`).
   2. The file `BFILE` (called `us_09066300_elev_bands.nc` in the catchment case study) describe the elevation bands required when the snow module is on. The dimensions of this file must match that of `forcefile`. The name of this file is made by appending the `suffix_elev_bands` defined in the `FILEMANAGER` (see B) to the basin ID (see E `$2`).
   
Note that the dimension of the NetCDF files will determine if FUSE is run at the catchment or grid-scale. FUSE will look for the variables `lat` and `lon` and if they are arrays, it will run on the grid they define. This means that NetCDF input files for a single catchment must also include the variables `lat` and `lon`.
   
## E. Execute FUSE

Run FUSE using default parameter values at the catchment scale:
```
./fuse.exe fm_catch.txt us_09066300 run_def
```

then calibrate it:

```
./fuse.exe fm_catch.txt us_09066300 calib_sce
```

then run it with the best SCE parameter set:

```
./fuse.exe fm_catch.txt us_09066300 run_best
```

The first arugment is the file manager, the second the basin ID and the third is the FUSE mode. Now you can increase the value of `MAXN` in `fm_catch.txt` to enable SCE to use more iterations and improve the model fit.

## F. Content of the output directory
Running FUSE in its different modes will create the following files in the `output` directory (provided for the catchment case study for comparison purposes):
* the files whose name contains `runs` contain the simulations, 
* the files whose name contains `para` contain the parameter values,
* the last element of the file name indicates which FUSE mode was used.

## G. Run FUSE for the grid case study

Download the data for the grid scale case study (see `FUSE modes and case studies` section) and modify `fm_grid.txt` following thesteps described in section C for `fm_catch.txt`. 

Run FUSE using default parameter values over the grid:

```
./fuse.exe fm_grid.txt cesm1-cam5 run_def
```

then using spatially distributed parameter values over the grid:

```
./fuse.exe fm_grid.txt cesm1-cam5 run_pre_dist cesm1-cam5_902_para_dist.nc
```

The file `cesm1-cam5_902_para_dist.nc` contains the distributed parameter values and must be located in the `output` folder. Note that with the current code i) the output parameter files (e.g., `*para_def.nc`) are empty and ii) it is not checked whether parameter values are within the prescribed realistic intervals.
   
## H. Compile SCE
The code of the shuffled complex evolution method (`$(MASTER)/build/FUSE_SRC/FUSE_SCE/sce.f`) was written in F77, so it must be compiled separately. If you use `ifort`, try the following flags:
  ```
  ifort -O2 -c -fixed sce_16plus.f
  ```

If necessary, rename the compiled file, so that it can be found by the `Makefile`, which by default will be looking for a file named `sce_16plus.o`.

### License
FUSE is distributed under the GNU Public License Version 3. For details see the file `LICENSE` in the FUSE root directory or visit the [online version](http://www.gnu.org/licenses/gpl-3.0.html).
