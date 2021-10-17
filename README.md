Introduction
-------
This Matlab package implements an optimization algorithm to find the best fuel surrogate matching a set of experimental properties, such as molecular weight, H/C ratio, cetane number, TSI (Threshold Soot Index), and distillation curves.

The methodology is inspired by the following papers:
* K. Narayanaswamy, H. Pitsch, P. Pepiot, Combustion and Flame, 165, p. 288-309 (2016)
* K. Narayanaswamy, P. Pepiot, Combustion Theory and Modelling, 22(5), p. 883-897 (2018)

It is written purely in Matlab language. It is self-contained. There is no external dependency.

Note: this package requires Matlab **R2016b** or later, since it utilizes a new series of Matlab features.

Installation
-------
1. Download the package to a local folder by running: 
```console
git clone https://github.com/acuoci/OpenSMOKEpp_SurrogateOptimizer.git
```
2. Run Matlab and navigate to the `OpenSMOKEpp_SurrogateOptimizer` folder, then run the `SurrogateOptimizer.m` script.
```console
SurrogateOptimizer
```

Notes
-------
This Matlab package implements an optimization algorithm to find the best fuel surrogate matching a set of experimental properties of a fuel.
The main script is ./src/OptimizeSurrogate.m

Properties currently implemented are:
- Molecular weight (MW)
- Hydrogen/Carbon ratio (HC)
- Density (rho)
- Viscosity (mu)
- Cetane number (CN)
- Yield sooting index (YSI)
- Distillation curve (DC)*
- Ignition delay times curve (IDT)*
- Laminar burning velocities curve (LBV)*

The user can set these properties as target (or any subset, just set the relative weight to zero).
*Curves used as target are not set in the main script, but on external files in the ./data folder

The Matlab package uses as input a database "./data/Database.xml" containing some chemical species and their properties.
You can expand it according to your available data.
Note that IDT and LBV are not computed exactly, but using models derived for some mixtures. For this reason, IDT and LBV curves can be used only with a limited palette of species.

The package can exploit parallel computing.
The optimization is carried out on two levels (be sure to have Optimization Toolbox available on your Matlab installation)
- a first rough scan with the genetic algorithm (ga)
- a series of deeper trials with patternsearch and/or fmincon
Few plots visualizations are available, as well as the possibility to save results to an output file

This tool is the result of Simone Pertesana's thesis work, under the supervision of Marco Mehl and Alberto Cuoci (Department of chemistry, materials and chemical engineering - Politecnico di Milano)


FeedBack
-------
If you find any bug or have any suggestion, please do file issues. We are graceful for any feedback and will do my best to improve this package.

Contributing
-------
1. Fork repository
2. Create a new branch, e.g. `checkout -b my-stuff`
3. Commit and push your changes to that branch
4. Make sure that the test works (!) (see known errors)
5. Create a pull request
6. I accept your pull request

Don't bunch up changes, e.g. if you have bug-fixes, new features and style changes, rather make 3 separate pull requests. Ensure that you introduce tests/examples for any new functionality

License
-------
Released under GPLv2 license.

Contact
-------
alberto.cuoci at polimi dot it
