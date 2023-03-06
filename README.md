# upp_cape_sandbox

no_upp: 
- this directory contains code to read a text file profile and calculate CAPE and CIN from UPP
- to build/run on hera: `module load intel/2022.2.0; cd no_upp; make; ./cape_driver.exe`

use_upp:
- this directory contains code to read a text file profile and calculate CAPE and CIN from UPP
- this code uses UPP code directly, except for UPP_PHYSICS.f
- to build/run on hera:
  - `module load intel/2022.2.0 netcdf/4.7.0 ncview`
  - `cd use_upp`
  - optional: edit `user_build_config` with path to UPP source
  - `make`
  - `cd run`
  - `./cape_driver.exe`
  - `ncview cape_output.nc` to view output
