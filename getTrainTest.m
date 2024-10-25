function[trainset,testset]=getTrainTest(data,pTrain)
% partions data into training and test sets
%INPUT  data = [ns,nd] dataset with ns samples and nd dimensions
% pTrain= percent of data to use for training
%output= trainst and testset

[ns,nd]=size(data);
nsTrain=round(pTrain*ns); %number of data samples for training
Ix=1:ns;%vector of indeces
IxTrain=randsample(Ix,nsTrain,true); %indicies for training 
trainset=data(IxTrain,:); % data for training
IxTest=ones(size(Ix)); %vector for binary ones
IxTest(IxTrain)=0; %remove the indicies to the training data
IxTest=logical(IxTest); %make IxTest a logical variable
testset=data(IxTest,:); %data for testing 



