%Feature extraction code prior to the training of classifier.

clc;clear all;

dlt=input('enter 0 to delete the feature excels. else enter any other number or character :'); %to prevent accidentally deleting feature files from previous run

if dlt==0
    delete('features.xlsx');
    delete('norm_features.xlsx');
end

[m2,m1,raw] = xlsread('autosegmented_syllables.xlsx');%load the filenames of segmented bird sound files

histbins=4;%number of bins for computing histogram

for i = 1:length(m1(:,1))
    [y,Fs]=wavread(char(m1(i,1)));%read i-th unsegmented wavefile

    birdclass(i)=m2(i);
    leny(i)=length(y);
    Freq(i)=Fs;
    Frame_size=ceil(0.003*Fs); %frame size of 3 ms 
    Frame_shift=floor(Frame_size/2); %overlap of 1.5 ms between two consecutive frames
    
    histvector=hist(y,histbins); %compute histogram
    
    c=melcepst(y, Fs, 'E0d'); %compute MFCC with Energy, 0th MFCC coefficient and delta coefficients - total length 28

    for j=1:28
        melvector(i,j)= mean(c(:,j));%mean MFCC for the whole syllable
    end

end

featmatrix = [melvector histvector]; % create the combined n*32 feature matrix where n=total number of wavefiles
[featrow, featcol]=size(featmatrix);
normfeatmatrix=zeros(featrow,featcol);

for a=1:featcol
  normfeatmatrix(:,a)=featmatrix(:,a)/max(featmatrix(:,a)); %normalize the feature matrix betwwen 0 and 1
end

features=[featmatrix birdclass']; %add the class vector as the rightmost clumn of the feature matrix, creating a n*33 matrix.
normfeatures=[normfeatmatrix birdclass'];

xlswrite('features.xlsx', features);
xlswrite('norm_features.xlsx', normfeatures);