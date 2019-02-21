/*
 * example.c - Driver code to demonstrate libmfcc working in real time
 * Author - Akarsh Prabhakara
 *
 * Parts of this code has been taken from 
 * https://github.com/EddieRingle/portaudio/tree/master/examples/paex_record.c for demo puroses only	
 */

#include <stdio.h>
#include <string.h> 
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include "portaudio.h"
#include "../src/libmfcc.h"

// MFCC parameters

#define nfilt 26			// number of filters 
#define numcep 13 			// number of cepstral coefficients per frame
#define winlen 0.025		// window / frame length in ms 
#define winstep 0.01 		// window / frame offset in ms
#define nfft 512			// N in N-point FFT
#define ceplifter 22		// 0 indicates no lifter 
#define preemph 0.97		// 0-1, 0 indicates no pre-emphasis
#define appendEnergy 1		// 0/1, 1 indicates replace first coefficient with the log of the energy in the frame
#define window 0			// 0/1/2, 0 - rectangular, 1 - Hanning, 2 - Hamming

// Audio input parameters

#define SAMPLE_RATE 16000
#define SAMPLES_PER_BUFFER 1024
#define NUM_SECONDS 1.5
#define WRITE_TO_FILE 1

// Select audio sample format

#if 0
#define PA_SAMPLE_TYPE  paFloat32
typedef float SAMPLE;
#define SAMPLE_SILENCE  (0.0f)
#define PRINTF_S_FORMAT "%.8f"
#elif 1
#define PA_SAMPLE_TYPE  paInt16
typedef short SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#elif 0
#define PA_SAMPLE_TYPE  paInt8
typedef char SAMPLE;
#define SAMPLE_SILENCE  (0)
#define PRINTF_S_FORMAT "%d"
#else
#define PA_SAMPLE_TYPE  paUInt8
typedef unsigned char SAMPLE;
#define SAMPLE_SILENCE  (128)
#define PRINTF_S_FORMAT "%d"
#endif

typedef struct
{
	int sampleIndex;
	int maxSampleIndex;
	SAMPLE *recordedSamples;
}
paTestData;

static int recordCallback(const void *inputBuffer, void *outputBuffer,
							unsigned long samplesPerBuffer,
							const PaStreamCallbackTimeInfo* timeInfo,
							PaStreamCallbackFlags statusFlags,
							void *userData)
{
	paTestData *data = (paTestData*)userData;
	const SAMPLE *rptr = (const SAMPLE*)inputBuffer;
	SAMPLE *wptr = &data->recordedSamples[data->sampleIndex];
	long samplesToCalc;
	int finished;
	unsigned long samplesLeft = data->maxSampleIndex - data->sampleIndex;

	(void) outputBuffer; 
	(void) timeInfo;
	(void) statusFlags;
	(void) userData;

	if(samplesLeft < samplesPerBuffer)
	{
		samplesToCalc = samplesLeft;
		finished = paComplete;
	}
	else
	{
		samplesToCalc = samplesPerBuffer;
		finished = paContinue;
	}

	if(inputBuffer == NULL)
		for(int i = 0; i < samplesToCalc; i++)
			*wptr++ = SAMPLE_SILENCE;
	else
		for(int i = 0; i < samplesToCalc; i++)
			*wptr++ = *rptr++;

	data->sampleIndex += samplesToCalc;
	return finished;
}

void calc_mfcc(SAMPLE *data, double **mfcc_result, int currentFrameID, int numframes, 
											int frame_size, int frame_offset, int numSamples)
{
	int startIndex = (currentFrameID - 1) * frame_offset;
	int endIndex = (currentFrameID != numframes) 
					? startIndex + frame_size - 1 
					: startIndex + numSamples - (currentFrameID - 1) * frame_offset - 1;
	
	double *local_data = (double *)malloc(sizeof(double) * (endIndex - startIndex + 1));

	// Pre-emphasis
	for (int i = endIndex; i >= startIndex + 1; i--)
		local_data[i - startIndex] = 1.0 * data[i] - data[i - 1] * preemph;
	local_data[0] = (currentFrameID == 1) 
					? 1.0 * data[startIndex] 
					: 1.0 * data[startIndex] - data[startIndex - 1] * preemph;

	// MFCC	
	mfcc(mfcc_result[currentFrameID - 1], local_data, 
		SAMPLE_RATE, nfilt, numcep, nfft, ceplifter, appendEnergy, window, endIndex - startIndex + 1);	

	free(local_data);
}

void main()
{
	PaStreamParameters  inputParameters,
						outputParameters;
	PaStream*			stream;
	PaError				err = paNoError;
	paTestData 			data;
	int					numSamples;
	int					numBytes;

	data.maxSampleIndex = numSamples = (int)ceil(NUM_SECONDS * (double)SAMPLE_RATE); 
	data.sampleIndex = 0;
	numBytes = numSamples * sizeof(SAMPLE);
	data.recordedSamples = (SAMPLE *) malloc(numBytes); 
	if (data.recordedSamples == NULL)
	{
		printf("Could not allocate record array.\n");
		goto done;
	}
	for (int i = 0; i < numSamples; i++) 
		data.recordedSamples[i] = 0;

	err = Pa_Initialize();
	if(err != paNoError) 
		goto done;

	inputParameters.device = Pa_GetDefaultInputDevice(); 
	if (inputParameters.device == paNoDevice) 
	{
		fprintf(stderr,"Error: No default input device.\n");
		goto done;
	}
	inputParameters.channelCount = 1;
	inputParameters.sampleFormat = PA_SAMPLE_TYPE;
	inputParameters.suggestedLatency = Pa_GetDeviceInfo(inputParameters.device)->defaultLowInputLatency;
	inputParameters.hostApiSpecificStreamInfo = NULL;
	
	if (numcep > nfilt || (winstep > winlen) || !((preemph >= 0) && (preemph <=1)) 
		|| !((appendEnergy == 0) || (appendEnergy == 1)) || !((window == 0) || (window == 1) || (window == 2))
		|| (nfilt <= 0)	|| (numcep <= 0) || (winlen <= 0) || (winstep <= 0) || (ceplifter < 0))
	{
		printf("Some error in MFCC parameters. Please check.\n");
		exit(0);
	}

	int frame_size = (int)round(winlen * SAMPLE_RATE);
	int frame_offset = (int)round(winstep * SAMPLE_RATE);
	int numframes;
	if (numSamples <= frame_size)
		numframes = 1;
	else
		numframes = 1 + (int)ceil((1.0 * numSamples - frame_size) / frame_offset);

	printf("\nFrame Size : %d, Frame offset : %d, Number of Frames : %d\n", 
												frame_size, frame_offset, numframes);
	if (nfft < frame_size)
	{
		printf("Some error in MFCC parameters. Please check.\n");
		exit(0);		
	}

	// Record some audio
	err = Pa_OpenStream(
						&stream,
						&inputParameters,
						NULL,
						SAMPLE_RATE,
						SAMPLES_PER_BUFFER,
						paClipOff,
						recordCallback,
						&data);

	if(err != paNoError) 
		goto done;

	err = Pa_StartStream(stream);
	if(err != paNoError) 
		goto done;

	printf("\n=== Now recording!! Please speak into the microphone. ===\n"); 

	double **mfcc_result = (double **)malloc(sizeof(double *) * numframes);
	for (int i = 0; i < numframes; i++)
		mfcc_result[i] = (double *)malloc(sizeof(double) * numcep);

	int currentFrameID = 1;
	while((err = Pa_IsStreamActive(stream)) == 1)
	{
		if (currentFrameID <= numframes && data.sampleIndex > (currentFrameID * frame_size))
		{
			calc_mfcc(data.recordedSamples, mfcc_result, currentFrameID, 
							numframes, frame_size, frame_offset, numSamples);
			currentFrameID += 1;
		}
	}

	if(err < 0) 
		goto done;

	err = Pa_CloseStream(stream);

	if(err != paNoError) 
		goto done;

	struct timeval tv1, tv2; 
	// gettimeofday(&tv1, 0);
	
	clock_t start, end;

	start = clock();

	// Continue with remaining mfccs
	while(currentFrameID <= numframes)
	{	
		calc_mfcc(data.recordedSamples, mfcc_result, currentFrameID, 
							numframes, frame_size, frame_offset, numSamples);
		currentFrameID += 1;
	}		

	// gettimeofday(&tv2, 0);
	end = clock();

	printf("Time taken to compute MFCCs after closing microphone stream - %.2f in number of clocks\n", (float)(end - start));
	// printf("Time taken to compute MFCCs after closing microphone stream - %.2f in milli seconds\n", (tv2.tv_usec - tv1.tv_usec) / 1000.0);
	printf("Clocks per sec %ld\n", CLOCKS_PER_SEC);

	// Write recorded data to a file. 
	#if WRITE_TO_FILE
	{
		FILE *fid;
		fid = fopen("recorded_data.dat", "w");
		for (int i = 0; i < numSamples; i++)
			fprintf(fid, "%d\n", data.recordedSamples[i]);
		fclose(fid);
		printf("Wrote data to 'recorded_data.dat'\n");

		fid = fopen("mfcc_result.dat", "w");
		for (int i = 0; i < numframes; i++)
			for (int j = 0; j < numcep; j++)
				fprintf(fid, "%.10e\n", mfcc_result[i][j]);
		fclose(fid);
		printf("Wrote MFCCs to 'mfcc_result.dat'\n");
	}
	#endif
	
	done:
	Pa_Terminate();
	if(data.recordedSamples)       
		free(data.recordedSamples);
	if(err != paNoError)
	{
		fprintf(stderr, "An error occured while using the portaudio stream\n");
		fprintf(stderr, "Error number: %d\n", err);
		fprintf(stderr, "Error message: %s\n", Pa_GetErrorText(err));
		err = 1;
	}

	for (int i = 0; i < numframes; i++)
		free(mfcc_result[i]);
	free(mfcc_result);
}