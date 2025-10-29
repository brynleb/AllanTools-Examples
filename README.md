# AllanTools-Examples
Examples of using the AllanTools Python package for time-series stability analysis with uncertainties.

[AllanTools](https://allantools.readthedocs.io/en/latest/) is a python library for calculating Allan deviation and related time & frequency statistics written by Anders Wallin.

## List of files:
- `ADevTools.py` - A module containing useful functions, including:
  - Computing common variants of Allan deviation ('ADev', 'Overlapping', 'Modified', 'Total', 'Modified Total').
  - Correct calculations of confidence intervals and uncertainties for these Allan deviations. 
  - Generating various types of colored noise ('Brownian', 'Violet', 'Pink', etc.).
- `ADevTools_Examples.py` - Example implementations of ADevTools functions, including:
  - Testing the response of different Allan deivations to simulated noise.
  - Computing the power spectral density (PSD) of time series data.
  - Understanding the PSD and Allan deviation of simulated data recorded at different sampling rates.

## Note about uncertainties
  - As of AllanTools release 2024.6, uncertainties and confidence intervals are not computed correctly for each type of Allan deviation and noise distribution. Presently, their code outputs standard errors assuming an idealized number of degrees of freedom ($N-1$).
  - This issue appears to be planned for a future release of AllanTools, but may require a significant change to their API (https://github.com/aewallin/allantools/issues/47)
  - Correct uncertainties are now a feature of **AllanToolKit**---a fork of the original AllanTools package.
 
## Additional resources:
- [AllanTools Github](https://github.com/aewallin/allantools)
- [AllanTools Documentation](https://allantools.readthedocs.io/en/latest/)
- [AllanToolKit GitLab](https://gitlab.com/amv213/allantoolkit)
- [A. Wallins blog post about confidence intervals](https://www.anderswallin.net/tag/allantools/)
- [NIST Stable32 software](http://www.stable32.com/)
- [NIST Handbook of Frequency Stability Analysis](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=50505)
