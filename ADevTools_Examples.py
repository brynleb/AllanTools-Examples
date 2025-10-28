## Import packages
import allantools as allan
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.signal import welch, periodogram

## Import custom modules
import ADevTools

## Set default plot options
plt.style.use('default')
mpl.rcParams['figure.dpi'] = 150
mpl.rcParams['savefig.dpi'] = 150
mpl.rcParams['font.size'] = 10
mpl.rcParams['legend.fontsize'] = 10

## Define sampling parameters

rates  = [5.0E3, 1.0E3, 2.0E2]	## Sample rates [Hz]
nrates = len(rates)

ndata  = [10001, 10001, 10001]	## Number of samples acquired at each sample rate
times  = [np.linspace(0., (ndata[i]-1)/rates[i], num=ndata[i]) for i in range(nrates)]

# ndata1 = 50001	## Number of samples at sample rate 1
# ndata2 = int((ndata1-1)*rate2/rate1) + 1	## Number of samples at sample rate 2
# ndata3 = int((ndata1-1)*rate3/rate1) + 1	## Number of samples at sample rate 3

# times1 = np.linspace(0., (ndata1-1)/rate1, num=ndata1)
# times2 = np.linspace(0., (ndata2-1)/rate2, num=ndata2)
# times3 = np.linspace(0., (ndata3-1)/rate3, num=ndata3)

S = lambda t, f, A, B: A*np.sin(2*np.pi*f*t) + B

rng    = np.random.default_rng(seed=1234)	## For repeatability
sigma  = 1.E-3 	## Standard deviation of Gaussian noise [V]
noises = [rng.normal(0., sigma, ndata[i]) for i in range(nrates)]
# noise = ADevTools.GenerateNoise(ndata1, rate1, NoiseType='White', NoisePar=1.E-3)

kwargs = {'f': 50.E0, 'A': 2.E-4, 'B': 1.E-3}
data   = [S(times[i], **kwargs) + noises[i] for i in range(nrates)]

## Compute PSDs
psd_type = 'welch'		## PSD type 'welch' or 'standard'
scaling  = 'spectrum' 	## Scaling of welch periodogram: 'density' for PSD, or 'spectrum' for squared magnitude spectrum
nsegs    = 8			## Number of segments to include in welch periodogram

freqs = [np.empty(1) for _ in range(nrates)]
asds  = [np.empty(1) for _ in range(nrates)]

for i in range(nrates):
	if psd_type == 'welch':
		freqs[i], psd = welch(data[i], rates[i], nperseg=int(ndata[i]/nsegs), scaling=scaling)
	else:
		freqs[i], psd = periodogram(data[i], rates[i], scaling=scaling)
	asds[i] = np.sqrt(psd)

## Compute Allan deviations
taus       = ['log10', 'log10', 'log10']	## Average times: 'all', 'octave', 'decade', 'log10'
adev_types = ['Total', 'Total', 'Total']	## Allan deviation type: 'ADev', 'Overlapping', 'Modified', 'Total', 'Modified Total'
comp_errs  = True

adev_taus = [np.empty(1) for _ in range(nrates)]
adev_devs = [np.empty(1) for _ in range(nrates)]
adev_errs = [np.empty(1) for _ in range(nrates)]

for i in range(nrates):
	(adev_taus[i], adev_devs[i], adev_errs[i]) = ADevTools.ComputeADev(data[i], taus=taus[i], rate=rates[i], ADevType=adev_types[i], ComputeErrors=comp_errs)

## Print signal statistics
# print('Time series:')
# print('Mean = {:.3e} V'.format(np.mean(noise1)))
# print('SDev = {:.3e} V'.format(np.std(noise1)))
# print('RMS  = {:.3e} V'.format(np.sqrt(np.mean(noise1**2))))
# if welch_scaling == 'density':
# 	print('Amplitude spectral density:')
# 	print('Mean = {:.3e} V/sqrt(Hz)'.format(np.mean(asd1)))
# 	print('SDev = {:.3e} V/sqrt(Hz)'.format(np.std(asd1)))
# 	print('RMS  = {:.3e} V/sqrt(Hz)'.format(np.sqrt(np.mean(asd1**2))))
# if welch_scaling == 'spectrum':
# 	print('Magnitude spectrum:')
# 	print('Mean = {:.3e} V'.format(np.mean(asd1)))
# 	print('SDev = {:.3e} V'.format(np.std(asd1)))
# 	print('RMS  = {:.3e} V'.format(np.sqrt(np.mean(asd1**2))))

## Plot results
fig = plt.figure(figsize=(2*4,2*2.5), layout="constrained")
gs  = GridSpec(2, 2, figure=fig, width_ratios=[1, 1], height_ratios=[1, 1])
ax0 = fig.add_subplot(gs[0, 0])
ax1 = fig.add_subplot(gs[1, 0])
ax2 = fig.add_subplot(gs[:, 1])
axs = [ax0, ax1, ax2]

colors = ['crimson', 'forestgreen', 'royalblue']

for i in range(nrates):
	axs[0].plot(times[i], data[i], color=colors[i], marker='', linestyle='-', alpha=0.5)
	axs[1].plot(freqs[i][1:], asds[i][1:], color=colors[i], marker='', linestyle='-', alpha=0.75)
	axs[2].errorbar(adev_taus[i], adev_devs[i], yerr=adev_errs[i], fmt='', ecolor=colors[i], color=colors[i], marker='.', linestyle='-', label=f'{rates[i]:.0f} Hz')

axs[0].set_xlabel('Time  (s)')
axs[0].set_ylabel('Signal  (V)')

axs[1].set_xlabel(r'$f$  (Hz)')
if scaling == 'density':
	axs[1].set_ylabel(r'ASD  (V/$\sqrt{Hz}$)')
if scaling == 'spectrum':
	axs[1].set_ylabel(r'Spectral Magnitude  (V)')
axs[1].set_xscale('log')
axs[1].set_yscale('log')
# axs[1].legend()

axs[2].set_xlabel(r'$\tau$  (s)')
axs[2].set_ylabel(r'$\sigma_A$  (V)')
axs[2].set_xscale('log')
axs[2].set_yscale('log')
axs[2].legend()

for ax in axs:
	ax.grid('both')

plt.show()