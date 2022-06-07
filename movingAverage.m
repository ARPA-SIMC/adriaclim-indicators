function newData = movingAverage(data,k,varargin)
% data1 = movingAverage(data,k);
% 
% Computes the moving average of vector 'data' forming averages of length 'k'
% 
% Input : 
% 	data        : data matrix 
% 	k           : size of the window of averaged blocks 
% 	varargin	: optional parameters
% 	'missing'  - Maximum rate of missing data within each window to be
% 	             considered not a NaN
 	
% Output :
% 	data1       : averaged data matrix
% 
% Example:
% 
%   dat=movingAverage(data,30,'missing',0.1);
%   mothly values for those months with less than 10 percent of missing data 

missing=0.1; %Maximum rate of missing data inside each box
for i=1:2:length(varargin)
    switch lower(varargin{i}),
        case 'missing', missing=varargin{i+1};
    end
end

% data=data(:);
[ndata,Nest]=size(data);
ini=1; fin=ndata;
ndata=floor((fin-ini+1)/k);
newData=zeros(ndata,Nest)*NaN; 
for j=1:ndata
    d=data((ini+(j-1)*k):(ini+j*k-1),:);
    aux=sum(isnan(d),1);
    id=find(aux/size(d,1)<=missing);
%     id=find(isnan(d)==0);
    newData(j,id)=nanmean(data((ini+(j-1)*k):(ini+j*k-1),id));
end
