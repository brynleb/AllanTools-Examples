# AllanTools-Examples
Examples of using the AllanTools Python package for time series stability analysis.

AllanTools is a python library for calculating Allan deviation and related time & frequency statistics written by Anders Wallin.

## List of files:
- `ADevTools.py` - A module containing useful AllanTools functions, including:
  - Computing common types of Allan deviation ('Overlapping', 'Modified', 'Total', etc.).
  - Calculating the statistical uncertainty of these Allan deviations.
  - Generating various types of colored noise ('Brownian', 'Violet', 'Pink', etc.).
- `ADevTools_Examples.py` - Example implementations of ADevTools functions, including:
  - Testing the response of different Allan deivations to simulated noise.
  - Computing the power spectral density (PSD) of time series data.
  - Understanding the PSD and Allan deviation of simulated data recorded at different sampling rates.

## Note about uncertainties
  - As of AllanTools release 2024.6, uncertainties and confidence intervals are not computed correctly for each type of Allan deviation. This code provides basic functions for computing asymmetric confidence intervals and the associated uncertainties for the 'Total' and 'Overlapped' deviations.
  - Note that this now seems to be a feature of **AllanToolKit**---a fork of the original AllanTools package.
 
## Additional resources:
- [AllanTools Github](https://github.com/aewallin/allantools)
- [AllanTools Documentation](https://allantools.readthedocs.io/en/latest/)
- [AllanToolKit GitLab](https://gitlab.com/amv213/allantoolkit)
- [A. Wallins blog post about confidence intervals](https://www.anderswallin.net/tag/allantools/)
- [NIST Stable32 software](http://www.stable32.com/)
- [NIST Handbook of Frequency Stability Analysis](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=50505)
