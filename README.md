# PeakFit
A set of peak fitting tools based on MATLAB for spectroscopic data analysis.

## Usage
Construct the PeakFit object in the following ways and the fit results will be populated in the readable object properties.

```MATLAB
obj=PeakFit(Data, ... Name-Value ...)
```

OR

```MATLAB
obj=PeakFit(XData, YData, ... Name-Value ...)
```

`Data` must be specified as a two-column (or two-row) matrix where the first column (or first row) is the X data points and the second column (or second row) is the Y data points. In the alternative syntax, `XData` and `YData` are respectively the X and the Y data points, specified as vectors, of the curve to be fitted.

The peak-fit settings can be specified after the mandatory arguments in Name-Value syntax, e.g. `PeakFit(Data, Window, [100,900], NumPeaks, 3)`. If no settings are specified, the algorithm will attempt to fit all peaks in the data using the default settings. The default behavior may not be optimal for noisy data. See [Best Practices](https://github.com/heriantolim/PeakFit#best-practices) for recommendations in making an optimal fit.

## Name-Value Pair Arguments

## Readable Properties

## Static Methods

## Best Practices
