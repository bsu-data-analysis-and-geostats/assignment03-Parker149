%% QUESTION 1 

clear all: close all; clc
D=load("icevelocity.txt")
Z=D(:,1);
V=D(:,2);
degree=[0 1 2 3 4];
for n=1:length(degree)
    P=polyfit(Z,V,degree(n)); %fit model 
    vmod(:,n)=polyval(P,Z); %evaluate model
    RMSE(:,n)=sqrt((mean(vmod(:,n)-V).^2)); %model accuracy
end 
plot(Z,V);%,'ko','linewidth',2,'MarkerSize',10)
hold on;
plot(Z,vmod(:,1),'linewidth',2);% Zero degree
plot(Z,vmod(:,2),'linewidth',2);% First degree
plot(Z,vmod(:,3),'linewidth',2);% Second degree
plot(Z,vmod(:,4),'linewidth',2);% Third degree
plot(Z,vmod(:,5),'linewidth',2);% Fourth degree
display(RMSE)
legend({'Raw Data','Zero (0.1249)','First (0.0468)','Second (0.0648)','Third (0.0586)','Fourth (0.0047)'})

%% QUESTION 2 

clear all: close all; clc
D=load("icevelocity.txt");

percent=round(91*.9); % Number of rows to sample

degree=[0 1 2 3 4];

for k=1:1000
    for n=1:length(degree)
        r=randsample(91,percent,false); % randomly select .9 rows of 91
        Dsampled=D(r,:); % what are the depths and velocities of the 82 random rows we found
        Z=Dsampled(:,1);
        V=Dsampled(:,2);
        P=polyfit(Z,V,degree(n));
        vmod(:,n)=polyval(P,Z);
        RMSE(:,n)=sqrt((mean(vmod(:,n)-V).^2));
        report(k, n) = RMSE(1); % Store the coefficients      
    end
end

Zero = [nanmean(report(:,1));nanstd(report(:,1))];
First = [nanmean(report(:,2));nanstd(report(:,2))];
Second = [nanmean(report(:,3));nanstd(report(:,3))];
Third = [nanmean(report(:,4));nanstd(report(:,4))];
Fourth = [nanmean(report(:,5));nanstd(report(:,5))];

Parameters = table(Zero,First,Second,Third,Fourth)


%% QUESTION 3 

clear all; close all; clc;

D= load('icevelocity.txt')

pTrain=0.9

for k=1:5
    for n=1:1000
        [trainset,testset]=getTrainTest(D,pTrain);
        ztrain=trainset(:,1);
        vtrain=trainset(:,2);
        P=polyfit(ztrain,vtrain,k);%fit training data with a linear model
        ztest=testset(:,1);
        vtest=testset(:,2);
        vmodel=polyval(P,ztest); %evaluate on test data
        rmsCV(n, k) = sqrt(mean((vmodel - vtest).^2));
    end
end

subplot(1,5,1)
hist(rmsCV(:,1))
title('degree zero')
subplot(1,5,2)
hist(rmsCV(:,2))
title('degree one')
subplot(1,5,3)
hist(rmsCV(:,3))
title('degree two')
subplot(1,5,4)
hist(rmsCV(:,4))
title('degree three')
subplot(1,5,5)
hist(rmsCV(:,5))
title('degree four')

%% QUESTION 4 

clear all: close all; clc
D=load("icevelocity.txt");
z=D(:,1);
v=D(:,2);
zmod=0:0.5:180;
winsize=[3,10,50];

for i=1:length(winsize)
    j=winsize(i); 
    vmod=move_window_ave(z,v,zmod,j);
    result(:, i) = vmod;      
end
zmod=zmod.';

plot(z,v,'ko','MarkerSize',10);
hold on
plot(zmod,result(:,1),'linewidth',2);
hold on
plot(zmod,result(:,2),'b','linewidth',2);
hold on
plot(zmod,result(:,3),'g','linewidth',2);
legend({'Raw Data','Winsize 3','Winsize 10','Winsize 50'});

%% QUESTION 5 

clear all: close all; clc
D=load("icevelocity.txt")
z=D(:,1);
v=D(:,2);


zmod=0:0.5:180
winsize=[3,10,50]; %input our window size

for i=1:length(winsize)
    j=winsize(i)
    %pulls a function we created eslewhere 
    vmod=nonparametric_smooth(z,v,zmod,j)
    result(:, i) = vmod;      
end
zmod=zmod.';

plot(z,v,'ko','MarkerSize',10);
hold on
plot(zmod,result(:,1),'linewidth',2);
hold on
plot(zmod,result(:,2),'b','linewidth',2);
hold on
plot(zmod,result(:,3),'g','linewidth',2);
legend({'Raw Data','Winsize 3','Winsize 10','Winsize 50'});

%% QUESTION 6 

D = load('icevelocity.txt');
pTrain = 0.9; 
zmod = 0:0.5:180; 
winsize = 3:1:100;

for i = 1:length(winsize)
    [trainset, testset] = getTrainTest(D, pTrain);
    ztrain = trainset(:, 1);
    vtrain = trainset(:, 2);
    vmod_train = nonparametric_smooth(ztrain, vtrain, zmod, winsize(i));
    vmodel_test = nonparametric_smooth(ztrain, vtrain, testset(:, 1), winsize(i));
    rmsCV(i) = sqrt(mean((vmodel_test - testset(:, 2)).^2, 'omitnan'));
end


plot(winsize, rmsCV);
xlabel('Window Size')
ylabel('RMSE')
title('Optimal Window Size')

%% QUESTION 7 & 8

clear all; close all; clc;

D = load('icevelocity.txt');
z = D(:,1);
v = D(:,2);

Arange = 1e-18:.1e-18:10e-18; % range for A0 from class. A should be 1x10^-18 (1x10^-18:0.1x10^-18:10x10^-18)
Nrange = 2:.01:4; % range for N0, looked it up and it seems like its 3 so we will do 1-5

%runs through values of A and n in the equtation and takes the RMSE for
%each attempted A and n value 


for n= 1:length(Arange)
    for f=1:length(Nrange)
        vmod=(v(1))-(Arange(n)*(917*9.8*sin(10*pi/180)).^Nrange(f)).*(z.^(Nrange(f)+1)); %equation for ice flow
        %in matlab for degrees we need to convert to radians so we need to
        %multiple our sin(10) by pi/180
        RMSE=sqrt(mean((vmod-v).^2));
        RMSE_matrix(n, f) = RMSE;
    end
end

imagesc(Arange,Nrange,RMSE_matrix,[0 20]);colorbar 
%we look at this graph and see that most of it tells us nothing so we set an threshold of a 
% of a RMSE of only look at 0-20 meters
%% QUESTION 9
clear all; close all; clc;

D = load('icevelocity.txt');
z = D(:,1);
v = D(:,2);

% we need to give the vars starting points, chose the min of the ranges found above
A0= 1e-18;

%First we outline what function fminsearch will run through and what
%variable to epxlore, next we need to give it our fuction with what
%vairables are what 
fun = @(A) gradient_search_method(A,z,v);

[A_min, fval] = fminsearch(fun, A0);

display(A_min)

%% QUESTION 10

%we are only looking for the optimal A value so when using the function
%this time we will just replace n with its optimal value of 1 

clear; close all; clc;
D = load('icevelocity.txt');
A0= 1e-18;
for n = 1:1000   
    ninper=round(size(D,1)*0.9); %this tells us what 90% of the length of D is 
    rows=randperm(size(D,1),ninper); % this will select 82 random rows from the length of D
    samp=D(rows,:); %this takes those random rows and finds what values those are within D
    z=samp(:,1);
    v=samp(:,2);
    fun = @(A) gradient_search_method(A,z,v);
    [A_min, fval] = fminsearch(fun, A0);
    Avals(n,:)=A_min;
    RMSE(n,:)=fval;
end
subplot(1, 2, 1);
histogram(Avals)
title('OPT A Vals');
subplot(1, 2, 2);
histogram(RMSE)
title('RMSE');

%% QUESTION 11

%pull code from 7 and 10, in the future need to do this in a better way
%rather than clear all; close all; clc; every time 

% code from Q7
clear all; close all; clc;
D = load('icevelocity.txt');
z = D(:,1);
v = D(:,2);
Arange = 1e-18:.1e-18:10e-18; 
Nrange = 2:.01:4; 
for n= 1:length(Arange)
    for f=1:length(Nrange)
        vmod=(v(1))-(Arange(n)*(917*9.8*sin(10*pi/180)).^Nrange(f)).*(z.^(Nrange(f)+1));
        RMSE=sqrt(mean((vmod-v).^2));
        RMSE_matrix(n, f) = RMSE;
    end
end

%code from Q10
A0= 1e-18;
for n = 1:1000   
    ninper=round(size(D,1)*0.9); 
    rows=randperm(size(D,1),ninper); 
    samp=D(rows,:); 
    z=samp(:,1);
    v=samp(:,2);
    fun = @(A) gradient_search_method(A,z,v);
    [A_min, fval] = fminsearch(fun, A0);
    Avals(n,:)=A_min;
    RMSE(n,:)=fval;
end

% Plot RMSE matrix using imagesc
imagesc(Nrange, Arange, RMSE_matrix, [0 20]); colorbar;
hold on;

% to be able to plot our boxchart along our y-axis(Arange) on the optimal N
% value of three we need a x variable of only 3 
three = 3 * ones(size(Avals)); 
boxchart(three, Avals,'BoxEdgeColor','red');

%% QUESTION 12

% run code from Q11 first 

%get the means and standard devations 
Astat=[nanmean(Avals),nanstd(Avals)];
RMSEstat=[nanmean(RMSE),nanstd(RMSE)];

%this is how we can generate 1k values using the mu and sd
x = Astat(1,1) + Astat(1,2) *randn(1000,1); %add the ,1 so that it is just 1 column with out that it will make a 1k by 1k
y = RMSEstat(1,1) + RMSEstat(1,2) *randn(1000,1);
%produces 1k simulated values for A as x and for RMSE as y

subplot(1, 2, 1);
histogram(x)
title('Simulated A Values');
subplot(1, 2, 2);
histogram(y)
title('Simulated RMSE Values');

%% QUESTION 13
% use matlabs kstest to compare out monte carlo values of A to our
% simulated values of A

h = kstest2(Avals,x)

%h=0 Fails to reject the null hypothesis (the samples are from the same distribution).
%h=1 Rejects the null hypothesis (the samples are from different distributions).



