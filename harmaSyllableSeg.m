function [syllables,FS,S,F,T,P] = harmaSyllableSeg(signal,FS,WINDOW,NOVERLAP,NFFT,MINDB)
% HARMASYLLABLESEG - Segments a signal stored in a WAV file into individual
% syllables.  Also graphs the spectrogram and signal with syllables
% highlighted in red to show what parts of the signal contain syllables.
%
% INPUT:
%
% - SIGNAL: mono signal in a one dimensional array
% - FS: Sampling frequency of mono signal
%   The following arguments are used by the spectrogram function type:
%   'help spectrogram' for more information on WINDOW,NOVERLAP, and NFFT
% - WINDOW: Either an integer value N or coefficients of a Window function
% stored in a length N matrix.  If an integer is passed a default hamming
% window of length N is used on each segment of the signal.
% - NOVERLAP: Number of samples each segment of the signal overlaps.
% - NFFT: Number of points used to calculate the DFT (discrete Fourier
% transform) of each segment.  This may be greater than the window length.
% In this case, each segment is zero padded to the NFFT length.
% - MINDB: Stopping criteria T (in dB) as defined in the original paper by Harma.
% A good default value for this parameter is ~20 dB.
%
% OUTPUT:
% - SYLLABLES: A struct array.  Each struct represents a single syllable and contains the following parameters:
%   - SIGNAL: An 1-dimensional array of doubles that represent the value of the signal over the range of this syllable. 
%   The following fields are in the order:
%   [Peak Peak-1 Peak-2...Peak+1 Peak+2...]
%   - SEGMENTS: The spectrogram index of each segment in this syllable.
%   - TIME: The time domain values of this syllable.
%   - FREQS: Peak frequency found in each segment.
%   - AMPS: Amplitude each peak frequency
% - FS: Sampling frequency of signal in the WAV file.
% - S: The spectrogram of the signal in the WAV file.
% - F: Frequency bins used in FFT.
% - T: Time domain values of each segment in the spectrogram.
% - P: Power spectral density of each segment in the spectrogram.
%
% Usage Example:
% [syllables,FS,S,F,T,P] = harmaSyllableSeg('[Path To WAV File]',kaiser(512),128,1024,20);
% 
% References:
% 1) Harma, A.; , "Automatic identification of bird species based on sinusoidal modeling of syllables,"
%    Acoustics, Speech, and Signal Processing, 2003. Proceedings. (ICASSP '03). 
%    2003 IEEE International Conference on , vol.5, no., pp. V- 545-8 vol.5, 6-10 April 2003
%    doi: 10.1109/ICASSP.2003.1200027
%    URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1200027&isnumber=26996
% 2) Lee, C. H., "Automatic Recognition of Bird Songs Using Cepstral Coefficients" 
%    Journal of Information Technology and Applications, 2006. Vol. 1 No. 1. p. 17-23
%    URL: http://140.126.5.184/Jita_web/publish/vol1_num1/05-20050044-text-sec.pdf
%
% Script By: Michael Lindemuth at the University of South Florida
%
% The code below may be licensed according to the BSD license.
% For a copy see: http://www.opensource.org/licenses/bsd-license.php

%% Extract Signal
%[signal,FS,~,~] = wavread(filename);
%signal = signal(:,1);

%% Calculate Spectrogram and Plot
[S,F,T,P] = spectrogram(signal,WINDOW,NOVERLAP,NFFT,FS);

% Note: plotting portion commented out because plotting syllables for every unsegmented file slows down the execution of the code. - Sutapa 

% figure;
% subplot(2,1,1);
% surf(T,F,10*log10(abs(P)),'EdgeColor','none');
mag = abs(S);
%surf(T,F,10*log10(mag),'EdgeColor','none');
% axis xy; axis tight; colormap(jet); view(0,90);
% xlabel('Time');
% ylabel('Frequency (Hz)');

%% Initialize Segmentation Parameters
N = 1;
syllables = [];
cutoff = [];
mag = abs(S);

%% Segment Signal into Syllables
while(1)
    % Find the maximum remaining magnitude in the spectrogram
    [freqMax,freqIndex] = max(mag,[],1);
    [argmax,segmentIndex] = max(freqMax,[],2);
    
    % Clear temp variables for this iteration
    times = [];
    segments = [];
    freqs = [];
    amps = [];
    
    % Setup temp variables with initial values
    segments(1) = segmentIndex;
    times(1) = T(segmentIndex);
    freqs(1) = F(freqIndex(segmentIndex));
    amps(1) = 20*log10(argmax);
    
    % Check if this is the first iteration, if so store the cutoff value
    % for the loop
    if(isempty(cutoff))
       cutoff = amps(1) - MINDB; 
    end
    
    % Is it time to stop looking for syllables?
    if(amps(1) < cutoff)
        break;
    end
    
    minAmp = amps(1) - MINDB;
    
    i = 1;
    %Look for all the values less than t with a high enough amplitude
    t = segmentIndex;
    while((t > 1) && (amps(i) >= minAmp))
        t = t-1; 
        i = i+1;
        segments(i) = t;
        times(i) = T(t);
        freqs(i) = F(freqIndex(t));
        amps(i) = 20*log10(freqMax(t));
    end
    
    % Remove the last index because it did not meet criteria
    if (i > 1)
        segments(i) = [];
        times(i) = [];
        freqs(i) = [];
        amps(i) = [];
        i = i-1;
    end
    
    % Look for all the values greater than t with a high enough amplitude
    t = segmentIndex;
    while((t < length(freqIndex)) && (amps(i) >= minAmp))
       t = t+1;
       i = i+1;
       segments(i) = t;
       times(i) = T(t);
       freqs(i) = F(freqIndex(t));
       amps(i) = 20*log10(freqMax(t));
    end
    
    % Remove the last index because it did not meet criteria
    if(i > 1)
        segments(i) = [];
        times(i) = [];
        freqs(i) = [];
        amps(i) = [];
        i = i-1;
    end
    
    % Store syllable parameters in struct
    syllables(N).signal = signal(round(min(times)*FS):round(max(times)*FS),1);
    syllables(N).powers = 10*log10(abs(P(:,sort(segments))));
    syllables(N).segments = segments;
    syllables(N).times = times;
    syllables(N).freqs = freqs;
    syllables(N).amps = amps;
    
    
    N = N+1;
    
    % Clear the magnitudes for this syllable so that it is not found again
    mag(:,segments) = 0;
end

%% Create plot of signal overlayed with segmented syllables in red markers
numSyllables = length(syllables);
% subplot(2,1,2)
% hold on;
% plot((1:length(signal))./FS,signal);
% for i=1:numSyllables
%    plot(sort(syllables(i).times),zeros(length(syllables(i).times),1),'Color','red','Marker','x');
% end
% hold off;
% axis tight;
% xlabel('Time (seconds)');
% ylabel('Signal Value');