# pyfocs

version==version==0.3.1.

pyfocs has been known by btmm_process (obscure non-pythonic name) and pyfox (an unmaintained package on PyPi) resulting in the new name for the library.

# Getting Started

## Installation

This installation assumes you have the anaconda distribution of python. If you do not have anaconda see the basic [troubleshooting section](#Troubleshooting).

### Using a package manager
pyfocs can be installed by using:

`pip install pyfocs`

which installs pyfocs plus all dependencies. This install method has caused problems for Windows OS. If you encounter errors when running pyfocs using this method, we instead recommend following the below method. If installing through pip and using an anaconda environment, then this step should be performed after installing all necessary conda packages, preferably in a new anaconda environment. [`pip` and `conda` can create conflicts in the packages, especially if used interchangeably multiple times.](https://www.anaconda.com/using-pip-in-a-conda-environment/).

### From source
Alternatively you can download the source code from this repository (green button with "Clone or Download"), extract the package, navigate to the directory containing it, and run:

`python setup.py install`

Note that Windows users will need to use anaconda power prompt or a similar python environment.

Both methods should result in the `PyFOX.py` being callable from the command line.

### Dependency issues
Please see the [troubleshooting section](#Troubleshooting). For issues not covered there please raise an issue and include details on your own system, your python version, and the versions of the packages pyfocs requires.

## Example

Download the data in the `example` directory. Within that directory is an example configuration file in yaml format. Adjust the `dir_pre` and `external` paths to be those of the example folder. Then, you should be able to run

`PyFOX.py path/to/example_configuration.yml`

Alternatively, providing no path to the yaml file will open a file browser for selecting the configuration file.

# Overview

The Bayreuth Micrometeorology python library for processing Fiber Optic Distributed Sensing (FODS) data. The library consists of a family of simple functions and a master script (`PyFOX`) that can be used to process output from a Silixa Distribute Temperature Sensing (DTS) device, such as an Ultima or XT, from the original `*.xml` files to calibrated temperatures with physical labels. This library is built around the [xarray](http://xarray.pydata.org) package for handling n-dimensional data, especially in a netcdf format.

## Other libraries

Other similar libraries exist, such as the [one developed at Delft University](https://github.com/bdestombe/python-geotechnical-profile), which can be more useful for some applications, especially those with double-ended configurations.

# PyFOX Steps

Data and the surrounding directory structure is assumed to follow this outline.
![](data_structure_scheme.jpg).

Each Subdirectory corresponds to a particular step in the processing.

1) Archives original `.xml` files into specified time interval.

2) Creates netcdfs of the raw data, including the instrument reported temperature, stokes intensity, and anti-stokes intensity. Dimensions of Length Along the Fiber, `LAF`, and time.

3) Labels the data, integrates external data streams and other reference data, performs step-loss corrections, performs single ended calibration based on Hausner et al., (2011). Splits multicore data into individual cores. Reports instrument reported temperature, calibrated temperature, log-power ratio of stoke and anti-stokes intensities, stokes intensity, anti-stokes intensities, and all data labels. Dimensions are `LAF` and `time`. New coordinates specified by location type in the location library can be used to label the data along with a `number of labels` by `number of LAF` coordinate.

4) Converts data labels with physical coordinates. Drops the LAF label and only includes the physical location (`xyz`) and `time`. Each `core` dimension is saved as a separate netcdf. Cores do not share the `xyz` dimension and must be aligned with each other. They do share the `time` dimension.

## Example jupyter notebook

For space reasons we only include the data for following steps 2-4 in the example notebook. The example notebook walks through the iterative approach for processing FODS data.

### References

Hausner, M. B., Suárez, F., Glander, K. E., & Giesen, N. Van De. (2011). Calibrating Single-Ended Fiber-Optic Raman Spectra Distributed Temperature Sensing Data. Sensors, 11, 10859–10879. https://doi.org/10.3390/s111110859

### Muppet Archiver

Batch script for scheduled archiving of `.xml` files on the Silixa DTS devices. Why muppet? Unviersity of Bayreuth Micrometeorology names their Silixa devices after muppet characters. Requires an anaconda 3.* distribution of python. Task scheduler must point to the `.bat` script and not the python script.

## <a name="Troubleshooting"></a>Troubleshooting/First time with python

### Install the anaconda version of python
https://www.anaconda.com/distribution/#download-section

You should be prompted to install python for your particular OS. Install version 3.7*.

## PyYAML updating
An error in an old version of pip stops the updating of PyYAML. If you get an error related to this library you can solve it with:

`pip install --ignore-installed PyYAML`

### Windows specific stuff
You will always need to use the anaconda power prompt or a similar environment. You will find it difficult/impossible to run python outside of the anaconda environment.

From https://www.correlatedsolutions.com/support/index.php?/Knowledgebase/Article/View/85/1/running-python-scripts-from-anywhere-under-windows

```
For the first time, windows asks what application to use to run *.py files.
Select python in Anaconda directory: make sure check box "Always use this app to open .py files" is checked
Click "More Apps"
Click "Look for another app on this PC"
Find path to anaconda python, e.g. "C:\Users\Your Name\Anaconda3\python"
Click "Open"
Click "OK"
```

Other issues may be related to running the anaconda prompt without administrator issues.
