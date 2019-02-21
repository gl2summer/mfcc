# Library for computing MFCC

This is an implementation of the standard [MFCC algorithm](https://en.wikipedia.org/wiki/Mel-frequency_cepstrum) in C using [FFTW library](http://www.fftw.org/).<br/>
[real_time_mfcc.c](Example/real_time_mfcc.c) demonstrates real time usage of this library. 

Clone this repository and follow the instructions below to run it successfully on a Linux machine (tested on Ubuntu 16.04 LTS).
## Dependencies

### FFTW library

Download the latest stable release of FFTW [here](http://www.fftw.org/download.html).<br/>
Open terminal from the folder where the tar.gz file was downloaded into.
```
tar xvf fftw-3.3.6-pl2.tar.gz
cd fftw-3.3.6-pl2/
./configure --enable-shared --enable-threads && make -j
sudo make install
```
Depending on your machine you may find the library in /usr/local/lib.<br/>
Go to mfcc local repository's root directory.
```
mkdir ./lib
cp /usr/local/lib/libfftw3.a ./lib/libfftw3.a
``` 
### Portaudio - Optional

Install this to successfully run [real_time_mfcc.c](Example/real_time_mfcc.c).
```
sudo apt-get install libasound-dev
```
Download the latest stable release of Portaudio [here](http://www.portaudio.com/download.html).<br/>
Open terminal from the folder where the tar.gz file was downloaded into.
```
tar xvfz pa_stable_v190600_20161030.tgz
cd portaudio/
./configure && make -j
sudo make install
```
Depending on your machine you may find the library in /usr/local/lib.<br/>
Go to mfcc local repository's root directory.
```
cp /usr/local/lib/libportaudio.a ./lib/libportaudio.a
``` 

### numpy, scipy, python_speech_features - Optional

Install this to check the correctness of output of [real_time_mfcc.c](Example/real_time_mfcc.c).
```
sudo apt-get install python-pip
pip install numpy
pip install scipy
pip install python_speech_features
```

## Usage

### Compile source code
```
cd ./src/
gcc -c libmfcc.c -o ../lib/libmfcc.o
```

### Generating static library
```
cd ../lib/
ar x libportaudio.a
ar x libfftw3.a
ar rcs libmfcc.a *.o
rm *.o
```
### Test run
```
cd ../Example/
gcc real_time_mfcc.c -L../lib/ -lmfcc -lm -lrt -lasound -pthread -o real_time_mfcc
./real_time_mfcc
```
### Checking correctness of output
```
python mse_check.py
```

## Author
Akarsh Prabhakara
