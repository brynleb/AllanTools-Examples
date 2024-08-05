# AllanTools-Examples
Examples for using the AllanTools Python package.

AllanTools is a python library for calculating Allan deviation and related time & frequency statistics written by Anders Wallin.



## List of examples:
- `ADev_Uncertainty.py`:
  - As of AllanTools release 2019.9, uncertainties and confidence intervals are not computed correctly for each type of Allan deviation. This code provides basic functions for computing asymmetric confidence intervals and the associated uncertainties for the 'Total' and 'Overlapped' deviations.
  - Note that this now seems to be a feature of **AllanToolKit**---a fork of the original AllanTools package.
 
## Additional resources:
- [AllanTools Github](https://github.com/aewallin/allantools)
- [AllanTools Documentation](https://allantools.readthedocs.io/en/latest/)
- [AllanToolKit GitLab](https://gitlab.com/amv213/allantoolkit)
- [A. Wallins blog post about confidence intervals](https://www.anderswallin.net/tag/allantools/)
- [NIST Stable32 software](http://www.stable32.com/)
