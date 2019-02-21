/*
 * libmfcc.c - Implementation of the standard MFCC algorithm in C using FFTW library for computing FFT
 * Author - Akarsh Prabhakara	
 *
 * MIT License
 * Copyright (c) 2017 Akarsh Prabhakara
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "fftw3.h"
#include "libmfcc.h"

void mfcc(double *mfcc_result, double *data, int samplingRate, int nfilt, int numcep, 
				int nfft, int ceplifter, int appendEnergy, int window, int numentries)

	// mfcc_result 		- stores the output of the MFCC computation - array of length numcep
	// data 			- input data of length numentries
	// samplingRate 	- the rate at which audio was sampled in Hz
	// nfilt 			- number of filters desired between 0 and samplingRate/2
	// numcep 			- number of MFCC coefficients
	// nfft 			- N in N-point FFT
	// ceplifter 		- 0 indicates no lifter 
	// window 			- 0/1/2, 0 - rectangular, 1 - Hanning, 2 - Hamming

{
	// Get frame data into temp_in and pad with zeros if numentries < nfft
	double *temp_in = (double *)malloc(nfft * sizeof(double));
	memset(temp_in, 0, nfft * sizeof(double));
	memcpy(temp_in, data, sizeof(double) * numentries);

	// Windowing
	if (window != 0)
		windowing(temp_in, numentries, window);

	// Run FFT
	fftw_plan p;
	fftw_complex *temp_out = (fftw_complex *)malloc((nfft / 2 + 1) * sizeof(fftw_complex));
	p = fftw_plan_dft_r2c_1d(nfft, temp_in, temp_out, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);

	free(temp_in);

	// Compute power spectra
	double *powspectrum = (double *)malloc((nfft / 2 + 1) * sizeof(double));
	for (int i = 0; i < nfft / 2 + 1; i++)
		powspectrum[i] = (1.0 / nfft) * (pow(temp_out[i][0], 2) + pow(temp_out[i][1], 2));

	free(temp_out);

	// Compute the filterbank parameters
	double **fbank = (double **)malloc(nfilt * sizeof(double*));
	for (int i = 0; i < nfilt; i++) 
		fbank[i] = (double *)malloc((nfft / 2 + 1) * sizeof(double));
	get_filterbank_parameters(fbank, nfilt, samplingRate, nfft);

	double specenergy = 0.0;
	for (int i = 0; i < nfft / 2 + 1; i++)
		specenergy += powspectrum[i];

	if (specenergy <= 0.0)
		specenergy = DBL_EPSILON;
	
	// Get filter bank output
	double *feat = (double *)malloc(nfilt * sizeof(double));
	for (int l = 0; l < nfilt; l++)
	{
		feat[l] = 0.0;
		for (int k = 0; k < nfft / 2 + 1; k++)
			feat[l] += powspectrum[k] * fbank[l][k];

		if (feat[l] > 0.0)
			feat[l] = log(feat[l]); 
		else
			feat[l] = DBL_EPSILON;	
	}

	free(powspectrum);

	for (int i = 0; i < nfilt; i++) 
		free(fbank[i]);	
	free(fbank);

	for (int i = 0; i < numcep; i++)
	{
		// DCT - II of filter bank output
		mfcc_result[i] = 0.0;
		for (int j = 0; j < nfilt; j++)
			mfcc_result[i] += feat[j] * cos((i * PI / nfilt) * (j + 0.5f));

		// Orthogonalization of DCT output
		if (i == 0)
			mfcc_result[i] *= sqrt(1.0 / nfilt);
		else
			mfcc_result[i] *= sqrt(2.0 / nfilt);

		// Ceplifter
		if (ceplifter != 0)
			mfcc_result[i] *= 1.0 + (ceplifter / 2.0) * sin(PI * i / ceplifter);
	}

	free(feat);

	// Append Energy
	if (appendEnergy == 1)
		mfcc_result[0] = log(specenergy);
}

void get_filterbank_parameters(double **fbank, int nfilt, int samplingRate, int nfft)
{
	double lowmel = hztomel(0.0);
	double highmel = hztomel(samplingRate / 2.0);

	// Generate nfilt center frequencies linearly spaced in the mel scale
	double *bin = (double *)malloc((nfilt + 2) * sizeof(double));
	for (int i = 0; i <= nfilt + 1; i++)
		bin[i] = floor(meltohz(i * (highmel - lowmel) / (nfilt + 1) + lowmel) * (nfft + 1) / samplingRate);

	// Triangular Filter Banks 
	for (int i = 0; i < nfilt; i++) 
	{
		memset(fbank[i], 0, (nfft / 2 + 1)*sizeof(double));
		for (int j = (int)bin[i]; j < (int)bin[i + 1]; j++)
			fbank[i][j] = (j - bin[i]) / (bin[i + 1] - bin[i]);
		for (int j = (int)bin[i + 1]; j < (int)bin[i + 2]; j++)
			fbank[i][j] = (bin[i + 2] - j) / (bin[i + 2] - bin[i + 1]);
	}
	free(bin);
}

double hztomel(double hz)
{
	return 2595 * log10(1 + hz / 700.0);
}

double meltohz(double mel)
{
	return 700 * (pow(10, mel / 2595.0) - 1);
}

void windowing(double *temp_in, int numentries, int window) 
{
	// Apply respective window before FFT
	for (int i = 0; i < numentries; i++)
	{
		if (window == 1)
			temp_in[i] *= (0.5 - 0.5 * cos(2 * PI * i / (numentries - 1)));
		else if (window == 2)
			temp_in[i] *= (0.54 - 0.46 * cos(2 * PI * i / (numentries - 1)));
	}
}