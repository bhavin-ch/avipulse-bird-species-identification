%Extract features from test data and classify them using the trained model.

trainingset=xlsread('train.xlsx');
trainclass=trainingset(:,33); %Column vector corresponding to the class of birds present in training set
traindata=trainingset(:,1:32);%Features matrix of training data

O = NaiveBayes.fit(traindata, trainclass);

[m2,m1,raw] = xlsread('test.xlsx');%Load the list of syallable files for testing
histbins=4;
for i = 1:length(m1(:,1))
[y,Fs]=wavread(char(m1(i,1)));%input: Speech signal
Frame_size=ceil(0.003*Fs);      %Input: Frame size in millisecond

Frame_shift=floor(Frame_size/2);    %Input: Frame shift in millisecond
        [rowhist, colhist]=size(y);
        if colhist>1
            h=hist(y,histbins)';
        else h=hist(y,histbins);
        end
for k=1:histbins
            histvector(1,k)=mean(h(:,k));
        end
        c=melcepst(y, Fs, 'E0d');
        for j=1:28
            melvector(1,j)= mean(c(:,j));
        end
featmatrix = [melvector histvector];
[featrow, featcol]=size(featmatrix);
normfeatmatrix=zeros(featrow,featcol);

for a=1:featcol
  normfeatmatrix(:,a)=featmatrix(:,a)/max(featmatrix(:,a));
end
end

probability=posterior(O,normfeatmatrix);
for i=1:length(test)
for j=1:16 %number of classes
if max(probalility(i))==probalility(i,j)
predicted(i)=j; %the most likely class for the i-th test syllable
end
end
end

accuracy=0;
predicted=predicted';

for i=1:length(test)
if predicted(i)==testset(i,33)
accuracy = accuracy + 1; %increase accuracy count by 1, if the predicted class and actual class of a test syllable match
end
end
Accuracy_percentage=100*accuracy/length(test)
