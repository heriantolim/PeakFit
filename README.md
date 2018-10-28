# PeakFit
**PeakFit** provides a tool to fit spectral data with a linear combination of symmetric peak functions such as Gaussian or Lorentzian. The background component of the spectra can be optionally taken into account by specifying the related parameter. In which case, the algorithm will attempt to fit the background with a polynomial function of the given order.

The fit model is given by:

f(*x*) ~ peak1(*x*) + peak2(*x*) + ... + peakN(*x*) + polynomial(*x*)

Each peak function is characterized by three fit coefficients: `Center`, `Width`, and either one of `Area` and `Height`. The polynomial function is characterized by n+1 fit coefficients, where n is the order of the polynomial.

## Licensing
This software is licensed under the GNU General Public License (version 3).

## Tested On
- MATLAB R2013b - R2018a

## Requirements
- MATLAB Curve Fitting Toolbox
- [MatCommon](https://github.com/heriantolim/MatCommon)
- [MatGraphics](https://github.com/heriantolim/MatGraphics) - required for doing the examples
- [PhysConst](https://github.com/heriantolim/PhysConst) - required for the `convertunit` method

## Setting Up
1. Download or git-clone this repository and other repositories listed in the [Requirements](https://github.com/heriantolim/PeakFit#requirements).
2. Add the repositories to the MATLAB's search path via `addpath(genpath( ... ))` OR this [version control system](https://github.com/heriantolim/MatlabVerCon).

## Usage
Construct the PeakFit object in the following ways, and the fit results will be populated in the object's [public properties](https://github.com/heriantolim/PeakFit#public-properties).

```MATLAB
obj=PeakFit(Data, ... Name-Value ...)
```

OR

```MATLAB
obj=PeakFit(XData, YData, ... Name-Value ...)
```

OR

```MATLAB
% Create an empty PeakFit object.
obj=PeakFit();

% Specify the data points and peak-fit settings via property assignments.
obj.XData= ... ;
obj.YData= ... ;
obj.Property1Name=Value1;
obj.Property2Name=Value2;
...

% Perform the peak fitting by reinstantiating the PeakFit object.
obj=PeakFit(obj);
```

`Data` must be specified as a two-column (or two-row) matrix where the first column (or first row) is the X data points and the second column (or second row) is the Y data points. In the alternative syntax, `XData` and `YData` are respectively the X and the Y data points, specified as vectors, of the curve to be fitted.

The peak-fit settings can be specified after the mandatory arguments in Name-Value syntax, e.g. `PeakFit(Data, 'Window', [100,900], 'NumPeaks', 3)`. If no settings are specified, the algorithm will attempt to fit all peaks in the data using the default settings. The default behavior may not be optimal for noisy data. See [Best Practices](https://github.com/heriantolim/PeakFit#best-practices) for recommendations in making an optimal fit.

### Name-Value Pair Arguments
Any of the [public properties](https://github.com/heriantolim/PeakFit#public-properties) can be specified as arguments in Name-Value syntax during the object construction. Specifying other things as Name-Value arguments will return an error.

## Examples
### Default Behavior
If the `PeakFit` is called without specifying the number of peaks, start points, or lower or upper bounds, then the algorithm will attempt to fit all peaks that it can it can guess. The following image is a photoluminescence spectrum of Erbium in Y<sub>2</sub>SiO<sub>5</sub> at near-liquid temperature. The spectrum was fitted using the command:
```MATLAB
Fit=PeakFit(Data,'PeakShape','Lorentzian');
```
The full code is given in [Examples/Er_PL_in_YSO.m](/Examples/Er_PL_in_YSO.m). It can be seen that many of the peaks were not resolved properly, and only the tallest peaks were correctly identified.

![Er PL in YSO](/Examples/Er_PL_in_YSO.png)

### Best Practices

## Public Properties
- `Data`, `XData`, `YData`: The data points of the curve to be fitted. Please ensure that the Y data points are all positive, otherwise the peak fitting may not work properly.
- `Window`: A vector of length two [*a*, *b*] that limits the fitting to only the data points whose X coordinates lies within [*a*, *b*].
- `NumPeaks`: The number of peaks wished to be fitted. When the fitting peak shapes, start points, lower, or upper bounds are set with vectors of length greater than `NumPeaks`, then `NumPeaks` will be incremented to adjust to the maximum length of these vectors. When the maximum length of these vectors is less, then these vectors will be expanded and filled with the default values. When `NumPeaks`=0 and all the start point, lower, and upper are not set, then the algorithm will attempt to fit all peaks that it can guess in the curve. Defaults to 0.
- `PeakShape`: A string vector that specifies the peak shape of each peak. The choices of `PeakShape` currently are: 'Lorentzian' (1) and 'Gaussian' (2). `PeakShape` may be set with an integer, the initial of these names, e.g. 'L' or 'G', or the full name, e.g. 'Lorentzian'. When the length of `PeakShape` is less than `NumPeaks`, the remaining peaks will be set with the default `PeakShape`, which is 'Lorentzian'. If PeakShape contains only one element, then the default value is the same as that element.
- `BaselinePolyOrder`: An integer that specifies the order of the polynomial function used to fit the background of the spectrum. Defaults to 0, which means a constant polynomial. Set this to a negative value to exclude the polynomial from the fit model.

### Fit Results
- (Read-only) `Peak`: A struct containing the fit results for each peak.
- (Read-only) `Base`: A struct containing the fit results for the baseline.
- (Read-only) `Area`, `Center`, `Height`, `Width`: A 3-by-`NumPeaks` matrix that stores the fit results for the area, center, height, and width, respectively; with each column correspond to the fit results for a particular peak. The first row holds the values at convergence, and the second (third) row holds the 95% CI lower (upper) bounds. Note that `Area` (or `Height`) is a redundant variable to the fit model, as it can be computed from `Width` and `Height` (or `Area`). Hence, it would be enough to specify the start points, lower or upper bounds for the `Area` or `Height` alone, as one can be computed from the other for any given `Width`.
- (Read-only) `Baseline`: A 3-by-(n+1) matrix, where n=`BaselinePolyOrder`, that stores the fit results for the baseline. The elements in the i-th column correspond to the fit results for the x^(n-i+1) polynomial coefficient. The first row holds the values at convergence, and the second (third) row holds the 95% CI lower (upper) bounds.
- (Read-only) `RelStDev`: The relative standard deviation of the fit results.
- (Read-only) `CoeffDeterm`: The coefficient of determination of the fit results.
- (Read-only) `AdjCoeffDeterm`: The degree-of-freedom adjusted coefficient of determination of the fit results.
- (Read-only) `NumFunEvals`: The number of function evaluations.
- (Read-only) `NumIters`: The number of iterations.
- (Read-only) `ExitFlag`: The exit condition of the algorithm. Positive flags indicate convergence, within tolerances. Zero flags indicate that the maximum number of function evaluations or iterations was exceeded. Negative flags indicate that the algorithm did not converge to a solution.

### Fitting Start Points
- `AreaStart`, `CenterStart`, `WidthStart`, `HeightStart`, `BaselineStart`: A vector of initial values for the area, center, width, height, and baseline coefficients, respectively. The default values are determined heuristically. To make certain properties default, set their values to NaN. They will be then replaced with the default values upon fitting.

### Fitting Lower Bounds
- `AreaLow`, `CenterLow`, `WidthLow`, `HeightLow`, `BaselineLow`: A vector of lower bounds for the area, center, width, height, and baseline coefficients, respectively. The default values are determined heuristically. To make certain properties default, set their values to -Inf. They will be then replaced with the default values upon fitting.

### Fitting Upper Bounds
- `AreaUp`, `CenterUp`, `WidthUp`, `HeightUp`, `BaselineUp`: A vector of upper bounds for the area, center, width, height, and baseline coefficients, respectively. The default values are determined heuristically. To make certain properties default, set their values to Inf. They will be then replaced with the default values upon fitting.

### Algorithm Parameters
- (Read-only) `Method` = 'NonLinearLeastSquares'; The method used for the fitting.
- `Robust`: The type of the least-squares method to be used in the fitting. Avaliable values are 'off', 'LAR' (least absolute residual method), and 'Bisquare' (bisquare weight method). Defaults to 'off'.
- `Algorithm`: The algorithm to be used in the fitting. Available values are 'Lavenberg-Marquardt', 'Gauss-Newton', or 'Trust-Region'. Defaults to 'Trust-Region'.
- `MovMeanWidth`: The window width of the moving average used for smoothing the curve in order to filter out the noise before finding the maximas. This parameter is used only when CenterStart is not given. The value of this can be set as a positive integer which specifies the width in terms of the number of data points, OR a real scalar between 0 and 1 which specifies the width as a fraction of the total number of data points. Defaults to 0.02.
- `DiffMaxChange`: The maximum change in coefficients for finite difference gradients. Defaults to 0.1.
- `DiffMinChange`: The minimum change in coefficients for finite difference gradients. Defaults to 1e-8.
- `MaxFunEvals`: The allowed maximum number of evaluations of the model. Defaults to 1e5.
- `MaxIters`: The maximum number of iterations allowed for the fit. Defaults to 1e3.
- `TolFun`: The termination tolerance on the model value. Defaults to 1e-6.
- `TolX`: The termination tolerance on the coefficient values. Defaults to 1e-6.

## Public Methods
- `disp`: Displays the options, results, error and performance of the fitting.
- `model`: Returns the reconstructed data points (model) using the fit results. See [@PeakFit/model.m](/@PeakFit/model.m) for more info.
- `convertunit`: Converts the units of the data points and the fit results. Available units to be converted from or to are: 'eV' (electron volt), 'percm'(wavenumber in cm^-1), 'Ramanshift' (Raman shift in cm^-1). When converting to or from Raman shift, an additional argument is required that specifies the excitation wavelength in nanometer. See [@PeakFit/convertunit.m](/@PeakFit/convertunit.m) for more info.

### Static Methods
- `fnlorentzian`: The Lorentzian function. See [@PeakFit/fnlorentzian.m](/@PeakFit/fnlorentzian.m) for more info.
- `fngaussian`: The Gaussian function. See [@PeakFit/fngaussian.m](/@PeakFit/fngaussian.m) for more info.

## See Also
- [BiErfFit](https://github.com/heriantolim/BiErfFit)
