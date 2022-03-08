# BARDATA

BARDATA is a modeling component that may be used to supply 2D "baroclinic" variables to ocean tide and surge models. These variables are computed from a global circulation model (currently GOFS 3.1) 3D variables (pressure, salinity, temperature and flow velocities if needed). The system contains the "ogcmdl" that performs the data download and the computation of the 2D variables.

BARDATA contains all the libraries required to build "ogcmdl" and the [NUOPC](https://earthsystemmodeling.org/nuopc/) files to create the "nuopc cap" library to used to couple BARDATA with any tide/surge model.

The software is developed by the [Coastal Marine modeling Branch(CMMB)](https://coastaloceanmodels.noaa.gov/) / [Office of Coast Survey](https://nauticalcharts.noaa.gov/), [National Ocean Service (NOS)](https://oceanservice.noaa.gov/) at [NOAA](https://www.noaa.gov/)

### Developer and Contact Information
Panagiotis Velissariou <Panagiotis.Velissariou@noaa.gov>

## Getting started

First, download the code by cloning this repo:

```
git clone https://github.com/noaa-ocs-modeling/BARDATA.git BARDATA
```

or by fetching the [zip archive](https://github.com/noaa-ocs-modeling/BARDATA/archive/refs/heads/master.zip) (unzip the archive to get the source code contained in the BARDATA folder).

You can then build BARDATA using the top level `build.sh` script.


## Build

You may run the `build.sh` script as:
```
build.sh --help
```
to get the "help" screen with all available options you may pass to the script.

