/*
 * libmfcc.h - Header file for libmfcc.c
 * Author - Akarsh Prabhakara	
 *
 * MIT License
 * Copyright (c) 2017 Akarsh Prabhakara
 */

#pragma once

#define DBL_EPSILON 2.2204460492503131e-16
#define PI 3.14159265358979323846264338327

void mfcc(double *mfcc_result, double *data, int samplingRate, int nfilt, int numcep, int nfft, int ceplifter, int appendEnergy, int window, int numentries);
void get_filterbank_parameters(double **fbank, int nfilt, int samplingRate, int nfft);
double hztomel(double hz);
double meltohz(double mel);
void windowing(double *temp_in, int numentries, int window);