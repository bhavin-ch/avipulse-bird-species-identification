%Code used for autosegmentation of syllables. Calls the harmaSyllableSeg
%function.

clc;clear all;close all;
[m2,m1,raw] = xlsread('unsegmented_ss.xlsx'); %Load names of unsegmented wavefiles

for i=1:length(m1(:,1))

   [y,Fs,bits]=wavread(char(m1(i,1)));%input: bird sound
   Freq(i)=Fs;
   hamming_length = ceil(Fs*0.003); % 3 ms Hamming window
   nooverlap = floor(hamming_length/2); % 1.5 ms overlap between consecutive windows
   mindb = 20;
      
   if length(y)<10000000 %to avoid MATLAB from hanging in case wavefile is too long - value 10000000 empirically chosen
   
       syllable = harmaSyllableSeg(y, Fs, hamming(hamming_length), nooverlap, 1024, mindb); %segmen the wavfile using Harma's algorithm
       l=[]; %empty array for storing length of each syllable
       stand=[]; %empty array for storing standard deviation of each syllable
       ampl=[]; %empty array for storing absolute mean ampltitude of each syllable
       [sylrow sylcol] = size(syllable);

       birdno=['Bird_' num2str(m2(i))]; %Class of bird
       fileprefix=[birdno '_' num2str(i)]; %index of file that is being segmented in current iteration
       for N = 1:sylcol
           x=syllable(N).signal;
           l(N)=length(x);%length of N-th syllable extracted from i-th wavefile
           stand(N)=std(x);%standard deviation of N-th syllable extracted from i-th wavefile
           ampl(N)=mean(abs(x));%mean of absolute amplitude of N-th syllable extracted from i-th wavefile
       end

       for N = 1:sylcol
            x=syllable(N).signal;
            if length(x)>=0.5*max(l) && std(x) >=mean(stand) && mean(abs(x)) >= mean(ampl)-std(ampl)%to rule out outliers
                filename=[fileprefix '_' num2str(N) '.wav'];     
                wavwrite(x, Fs, filename);%save wavefiles that are presumably not outliers
            end
       end
       
   else %For files with more than 10000000 samples : divide them into two halves and segment seperately.
        y1=y;
        y=y1(1:round(0.5*length(y1)));

        syllable = harmaSyllableSeg(y, Fs, hamming(hamming_length), nooverlap, 1024, mindb);
        l=[];
        stand=[];
        ampl=[];
        [sylrow sylcol] = size(syllable);

        birdno=['Bird_' num2str(m2(i)) '_1'];
        fileprefix=[birdno '_' num2str(i) ];
        for N = 1:sylcol
            x=syllable(N).signal;
            l(N)=length(x);
            stand(N)=std(x);
            ampl(N)=mean(abs(x));
        end

        for N = 1:sylcol
            x=syllable(N).signal;
            if length(x)>=0.5*max(l) && std(x) >=mean(stand) && mean(abs(x)) >= mean(ampl)-std(ampl)
                filename=[fileprefix '_' num2str(N) '.wav'];        
                wavwrite(x, Fs, filename);
            end
        end
        y=y1(round(0.5*length(y1)):length(y1));

        syllable = harmaSyllableSeg(y, Fs, hamming(hamming_length), nooverlap, 1024, mindb);
        l=[];
        stand=[];
        ampl=[];
        [sylrow sylcol] = size(syllable);

        birdno=['Bird_' num2str(m2(i)) '_2' ];
        fileprefix=[birdno '_' num2str(i)];
        
        for N = 1:sylcol
            x=syllable(N).signal;
            l(N)=length(x);
            stand(N)=std(x);
            ampl(N)=mean(abs(x));
        end

        for N = 1:sylcol
            x=syllable(N).signal;
            if length(x)>=0.5*max(l) && std(x) >=mean(stand) && mean(abs(x)) >= mean(ampl)-std(ampl)
                filename=[fileprefix '_' num2str(N) '.wav'];        
                wavwrite(x, Fs, filename);
            end
        end
    end
end
