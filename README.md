# PeakFit
A set of peak fitting tools based on MATLAB for spectroscopic data analysis.

## Requirements
* MATLAB Curve Fitting Toolbox
* [PhysConst](https://github.com/heriantolim/PhysConst) - required when invoking the `convertunit` method


## Usage
Construct the PeakFit object in the following ways and the fit results will be populated in the [public properties](https://github.com/heriantolim/PeakFit#public-properties) of the object.

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

The peak-fit settings can be specified after the mandatory arguments in Name-Value syntax, e.g. `PeakFit(Data, Window, [100,900], NumPeaks, 3)`. If no settings are specified, the algorithm will attempt to fit all peaks in the data using the default settings. The default behavior may not be optimal for noisy data. See [Best Practices](https://github.com/heriantolim/PeakFit#best-practices) for recommendations in making an optimal fit.

## Name-Value Pair Arguments


## Public Properties

## Public Methods

## Best Practices

## See Also

