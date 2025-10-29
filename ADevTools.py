import numpy as np
import allantools as allan
from scipy.stats import chi2

## alpha = exponent of f in the frequency PSD
## NoiseType		Description					alpha
##---------------------------------------------------
## 'White PM'		White phase noise           +2
## 'Flicker PM'		Flicker phase noise         +1
## 'White FM'		White frequency noise        0
## 'Flicker FM'		Flicker frequency noise     -1
## 'Random Walk FM'	Random walk frequency noise	-2

ALPHAS = {'White PM': 2, 'Flicker PM': 1, 'White FM': 0, 'Flicker FM': -1, 'Random Walk FM': -2}

#####################################################################

def ComputeADev(data, taus='octave', rate=1., ADevType='Total', tauMax=0.4, ComputeErrors=True, NoiseType='White FM'):
	"""Compute Allan deviation of data.
	ARGUMENTS:
	data      (np.array) - Input data. Assumed to be a unitless fractional-frequency-like data.
	taus      (np.array) - Values of averaging time tau (s) at which to compute ADev. 
						   One can also use values of 'all', 'octave', 'decade', or 'log10'. Defaults to 'octave'.
	rate         (float) - Sampling rate of input data (Hz). Defaults to 1 Hz.
	ADevType       (str) - Control for type of Allan deviation to compute ('ADev', 'Overlapping', 'Modified', 'Total', 'Modified Total').
						   Defaults to 'Total'.
	tauMax       (float) - Maximum tau to plot (in units of the half data length when ADevType = 'Overlapping',
						   and the full data length when ADevType = 'Total'). Defaults to 0.4.
	ComputeErrors (bool) - Flag for computing ADev uncertainties based on Chi2 distribution (expensive for large data sets).
	NoiseType      (str) - Type of noise to assume for uncertainties calculation.
	RETURN FORMAT: (taus2, ADevs, ADevLowerErrors, ADevUpperErrors)
	taus2      (np.array) - Averaging times in seconds at which Allan deviations were computed.
	ADevs      (np.array) - Allan deviations.
	ADevErrors (np.array) - Statistical uncertainties on ADevs (lower, upper).
	"""

	if ADevType == 'ADev':
		## Classic - use only if required - relatively poor confidence.
		(taus2, ADevs, ADevErrors, ADevNs) = allan.adev(data, rate=rate, data_type='freq', taus=taus)
	elif ADevType == 'Overlapping':
		## General purpose - most widely used.
		(taus2, ADevs, ADevErrors, ADevNs) = allan.oadev(data, rate=rate, data_type='freq', taus=taus)
	elif ADevType == 'Modified':
		## Used to distinguish between White and Flicker Phase Modulation.
		(taus2, ADevs, ADevErrors, ADevNs) = allan.mdev(data, rate=rate, data_type='freq', taus=taus)
	elif ADevType == 'Total':
		## Better confidence at long averages for Allan deviation.
		(taus2, ADevs, ADevErrors, ADevNs) = allan.totdev(data, rate=rate, data_type='freq', taus=taus)
	elif ADevType == 'Modified Total':
		## Better confidence at long averages for modified Allan
		(taus2, ADevs, ADevErrors, ADevNs) = allan.mtotdev(data, rate=rate, data_type='freq', taus=taus)

	## Ensure tauMax is less than or equal to 1
	tauMax = min(abs(tauMax), 1.0)
	## Total number of samples
	ndata = len(data)
	## Number of averaging times requested
	ntau = len(taus2)

	## Keep only a subset of Allan deviation
	if taus == 'all':
		ntau = int(round(tauMax*ntau))
	elif taus == 'octave' or taus == 'decade' or taus == 'log10':
		itau = 0
		while taus2[itau] <= tauMax*taus2[-1]:
			itau += 1
		ntau = itau+1

	taus2  = taus2[:ntau]
	ADevs  = ADevs[:ntau]
	ADevNs = ADevNs[:ntau]
	ADevErrors = ADevErrors[:ntau]

	if ComputeErrors:
		ADevLowerErrors, ADevUpperErrors = ComputeADevErrors(ndata, ntau, rate, taus2, ADevs, ADevErrors, ADevType, NoiseType)[:2]
		ADevErrors = np.vstack((ADevLowerErrors, ADevUpperErrors))
	else:
		ADevErrors = np.vstack((ADevErrors, ADevErrors))

	return (taus2, ADevs, ADevErrors)

######################## End of AllanDev() ##########################
#####################################################################

def ComputeADevErrors(ndata, ntau, rate, taus, ADevs, ADevErrors, ADevType='Total', NoiseType='White FM', ModelType='chi2'):
	"""Compute 1-sigma confidence intervals and one- or two-sided uncertainties for Allan deviations.
	ARGUMENTS:
	ndata           (int) - Total number of samples
	ntau            (int) - Number of averaging times contained in 'tau'
	rate          (float) - Sampling rate of input data (Hz)
	taus       (np.array) - Averaging times
	ADevs      (np.array) - Allan deviation
	ADevErrors (np.array) - One-sided uncertainty output by AllanTools
	ADevType        (str) - Type of Allan deviation
	NoiseType       (str) - Type of noise present in data ('White PM', 'Flicker PM', 'White FM', 'Flicker FM', 'Random Walk FM')
	ModelType       (str) - Type of model to use to compute uncertainty ('chi2', otherwise AllanTools default)
	RETURN FORMAT: (ADevLowerErrors, ADevUpperErrors, ADevLowerCIs, ADevUpperCIs)
	ADevLowerErrors (np.array) - Lower bound of ADev uncertainties
	ADevUpperErrors (np.array) - Upper bound of ADev uncertainties
	ADevLowerCIs    (np.array) - Lower bound of ADev confidence intervals
	ADevUpperCIs    (np.array) - Upper bound of ADev confidence intervals
	"""

	alpha = ALPHAS[NoiseType]

	if ModelType == 'chi2':
		## Compute Allan deviation confidence intervals based on chi2-distribution
		ADevLowerErrors = np.zeros(ntau)
		ADevUpperErrors = np.zeros(ntau)
		ADevLowerCIs    = np.zeros(ntau)
		ADevUpperCIs    = np.zeros(ntau)

		for i in range(ntau):
			## Averaging factor tau = m*tau0 = m/rate
			m = rate*taus[i]
			N = float(ndata)

			if ADevType == 'ADev' or ADevType == 'Overlapping' or ADevType == 'Modified':
				## Equivalent degrees of freedom for overlapping Allan deviation
				## See NIST Handbook of Frequency Stability Analysis, Table 5, p. 40
				if alpha == 2:
					## White phase noise
					edf = (N+1)*(N-2*m)/(2*(N-m))
				elif alpha == 1:
					## Flicker phase noise
					a = np.log((N-1.)/(2.*m))
					b = np.log((2.*m+1.)*(N-1.)/4.)
					edf = np.sqrt(np.exp(a*b))
				elif alpha == 0:
					## White frequency noise
					a = 3.*(N-1.)/(2.*m) - 2.*(N-2.)/N
					b = 1. + 5./(4.*m**2)
					edf = a/b
				elif alpha == -1:
					## Flicker frequency noise
					if m == 1:
						edf = 2.*(N-2.)**2/(2.3*N - 4.9)
					else:
						edf = 5.*N**2/(4*m*N+12.*m**2)
				elif alpha == -2:
					## Random walk frequency noise
					a = (N-2.)/m
					b = ((N-1.)**2 - 3*m*(N-1.) + 4.*m**2)/(N-3.)**2
					edf = a*b

				# edf = allan.edf_simple(N, m, alpha)

				# if ADevType == 'ADev':
				# 	overlapping, modified = False, False
				# elif ADevType == 'Overlapping':
				# 	overlapping, modified = True, False
				# else:
				# 	overlapping, modified = False, True
				# edf = allan.edf_greenhall(alpha, d=2, m=m, N=N, overlapping=overlapping, modified=modified)

			elif ADevType == 'Total':
				## Equivalent degrees of freedom for total deviation
				## See NIST Handbook of Frequency Stability Analysis, Table 7, p. 41
				## Note that this table only includes alpha = 0, -1, -2 cases... 
				if alpha == 0:
					## White frequency noise
					# edf = 1.5*float(ndata)/m
					b, c = 1.50, 0.
				elif alpha == -1:
					## Flicker frequency noise
					b, c = 1.17, 0.22
				elif alpha == -2:
					## Random walk frequency noise
					b, c = 0.93, 0.36
				edf = b*N/m - c
				# edf = allan.edf_totdev(N, m, alpha)

			elif ADevType == 'Modified Total':
				## Equivalent degrees of freedom for modified total deviation
				## See NIST Handbook of Frequency Stability Analysis, Table 8, p. 41
				if alpha == 2:
					## White phase noise case
					b, c = 1.90, 2.10
				elif alpha == 1:
					## Flicker phase noise case
					b, c = 1.20, 1.40
				elif alpha == 0:
					## Random walk frequency noise case
					b, c = 1.10, 1.20
				elif alpha == -1:
					## Flicker frequency noise case
					b, c = 0.85, 0.50
				elif alpha == -2:
					## Random walk frequency noise case
					b, c = 0.75, 0.31
				edf = b*N/m - c
				# edf = allan.edf_mtotdev(N, m, alpha)

			## Chi-Squared values corresponding to +/- sigma
			(chi2Lower, chi2Upper) = ChiSquaredModel(edf)
			## Estimate two-sided confidence interval
			ADevLowerCIs[i]    = np.sqrt(edf/chi2Upper)*ADevs[i]
			ADevUpperCIs[i]    = np.sqrt(edf/chi2Lower)*ADevs[i]
			ADevLowerErrors[i] = np.abs(ADevLowerCIs[i] - ADevs[i])
			ADevUpperErrors[i] = np.abs(ADevUpperCIs[i] - ADevs[i])
	else:
		## Simple one-sided confidence interval output by allantools (typically an underestimate of true confidence interval)
		##   ADevErr = ADev/sqrt(n), where n is the number of data pairs used to compute each ADev
		ADevLowerErrors = ADevErrors
		ADevUpperErrors = ADevErrors
		ADevLowerCIs    = ADevs - ADevErrors
		ADevUpperCIs    = ADevs + ADevErrors

	return (ADevLowerErrors, ADevUpperErrors, ADevLowerCIs, ADevUpperCIs)

################### End of ComputeADevErrors() ######################
#####################################################################

def ChiSquaredModel(dof, p=0.682689):
	"""Compute chi-squared distribution for a given number of degrees of freedom (dof) at confidence level p.
	ARGUMENTS:
	dof   (float) - Number of degrees of freedom
	p     (float) - Value between 0 and 1 corresponding to confidence level.
					Defaults to erf(1/np.sqrt(2)) = 0.682689, corresponding to 1 sigma.
	RETURN FORMAT: (chi2Lower, chi2Upper)
	chi2Lower (float) - Lower limit of chi2 distribution
	chi2Upper (float) - Upper limit of chi2 distribution
	"""

	chi2Lower = chi2.ppf(0.5*(1-p), dof) ## Percent point function (inverse of cdf)
	chi2Upper = chi2.ppf(0.5*(1+p), dof)

	return (chi2Lower, chi2Upper)

#################### End of ChiSquaredModel() #######################
#####################################################################

def GenerateNoise(nData, Rate, NoiseColor='White', NoisePar=1.):
	"""Generate an array of colored noise.
	ARGUMENTS:
	nData      (int) - Number of samples
	Rate     (float) - Sampling rate (Hz)
	NoiseColor (str) - Color of noise to generate ('White', 'Brown', 'Violet', 'Pink')
					   Defaults to 'White'.
	NoisePar (float) - Noise parameter specific to type. Defaults to 1.
 	RETURN FORMAT: Noise
	Noise (np.array) - Generated noise values.
	"""

	if NoiseColor == 'White':
		## Generate time series with white noise that has constant PSD = b0, up to the nyquist frequency fs/2.
		b0 = NoisePar ## Power spectral density of white noise [X^2/Hz]
		Noise = allan.noise.white(nData, b0, Rate)
	elif NoiseColor == 'Brown':
		## Brownian or random walk (diffusion) noise with 1/f^2 PSD.
		## Not really a color, rather Brownian or random-walk. Obtained by integrating white noise.
		## The desired power-spectral density is b_minus2*f^-2.
		b_minus2 = NoisePar
		Noise = allan.noise.brown(nData, b_minus2, Rate)
	elif NoiseColor == 'Violet':
		## Violet noise with f^2 PSD. Obtained by differentiating white noise
		## The Desired power-spectral density is b2*f^2
		b_2 = NoisePar
		Noise = allan.noise.violet(nData, b_2, Rate)
	elif NoiseColor == 'Pink':
		## Pink noise (approximation) with 1/f PSD
		## depth: number of iterations for each point. High numbers are slower but generates a more correct spectrum at low frequencies.
		depth = NoisePar
		Noise = allan.noise.pink(nData, depth)

	return Noise

##################### End of GenerateNoise() ########################
#####################################################################