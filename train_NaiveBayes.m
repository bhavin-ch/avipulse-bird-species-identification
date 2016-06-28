clc
clear all

%% Extract outlier-removed data

[num, txt, raw]=xlsread('norm_features.xlsx');
[numrow, numcol] = size(num);

X=num(:,1:(numcol-1));
Y=num(:,numcol);

%following loop introduces a very small difference in case for a particular
% class, one or more features have same value in all the cells. This allows
% the Naive Bayse classifier to function properly.
for a = 1:(numcol-1)
    for b=1:(numrow-1)
    if Y(b)==Y(b+1)
        if X(b,a)== X(b+1,a)
            if X(b,a)==0
                X(b+1,a)=0.00001;
            else
                X(b+1,a)=1.00001*X(b+1,a);
            end
        end
    end
    end
end

%% Fit Naive Bayes Model

O = NaiveBayes.fit(X,Y);
C = O.predict(X);
cmat = confusionmat(Y,C);
accuracy=trace(cmat)/sum(sum(cmat))%training accuracy
count=ones(size(cmat,1),1);
j=1;
number(1)=1;
for i=2:size(Y)
    if Y(i)==Y(i-1)
        number(j)=number(j)+1;
    else 
        j=j+1;
        number(j)=1;
    end
end

for i=1:size(cmat,1)
    if cmat(i,i)==number(i)
        error(i)=0;
    else
        error(i)=number(i)-cmat(i,i);
    end
end

datapoints=sum(number)