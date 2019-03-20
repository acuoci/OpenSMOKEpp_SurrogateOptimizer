Introduction
-------
This Matlab package implements an optimization algorithm to find the best fuel surrogate matching a set of experimental properties, such as molecular weight, H/C ratio, cetane number, TSI (Threshold Soot Index), and distillation curves.

The methodology is inspired by the following papers:
* K. Narayanaswamy, H. Pitsch, P. Pepiot, Combustion and Flame, 165, p. 288-309 (2016)
* K. Narayanaswamy, P. Pepiot, Combustion Theory and Modelling, 22(5), p. 883-897 (2018)

It is written purely in Matlab language. It is self-contained. There is no external dependency.

Note: this package requires Matlab **R2016b** or latter, since it utilizes a new series of Matlab features.

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

Don't bunch up changes, e.g. if you have bug-fixes, new features and style changes, rather make 3 seperate pull requests. Ensure that you introduce tests/examples for any new functionality

License
-------
Released under GPLv2 license.

Contact
-------
alberto.cuoci at polimi dot it
