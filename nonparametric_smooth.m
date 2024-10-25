function ymod=nonparametric_smooth(x,y,xmod,winsize)
%this function smooths a data set of 1 variable using a weighted moving
%average
%inputs x=independent variable (n,1)
%       y= depednet varaible (n,1)
%   xmod = location for estimates (*,1)
%  winsize= size of half the window (1,1)
%OUTPUT = ymod = non-parameteric smooth esimates

x=x(:);y=y(:);xmod=xmod(:); %force to be columns
ymod=zeros(size(xmod))*NaN; % initialize
for i=1:length(xmod)
    dist=sqrt((x-xmod(i)).^2); %distance function 
    ival=find(dist<winsize);
    if length(ival)>0
        weights=15/16*(1-(dist(ival)/winsize).^2).^2; %bisquare kernel
        ymod(i)=sum(weights.*y(ival))./sum(weights);%modeled ests
    end
end