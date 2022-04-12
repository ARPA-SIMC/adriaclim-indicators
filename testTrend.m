function [pValue,trend,T]=testTrend(data,varargin)
% [pValue,trend]=testTrend(data,varargin);
% 
% The function computes the pvalue for the data trend
% 
% Input data:
% 	data        : data to be used to compute the trend
% 	varargin	: poptional parameters
% 		'test'	- 'Spearman' or 'MannKendall' (MannKendall as default)
% 		'period'-'day','month' y 'year'. How to group the data to compute the trend. With default id uses 'day'.
%       'missing' - missing data parameter for the 'movingAverage' function
% 	
% Output data :
% 	pValue : p values
%            
% 	trend  : trend values
% 	T        : vector containing the value of the statistics
% 
% Example:
% 
% 		[pValue,trend]=testTrend(data,'test','Spearman','period','year','missing',0.1);

test='MannKendall';
period='day';
missing=0.1;
autocorrelation=0;
for i=1:2:length(varargin)
    switch lower(varargin{i}),
        case 'test', test=varargin{i+1};
        case 'period', period=varargin{i+1};
        case 'missing', missing=varargin{i+1};
        case 'autocorrelation', autocorrelation=varargin{i+1};
    end
end

Nest=size(data,2);
trend=zeros(Nest,1)+NaN;
pValue=zeros(Nest,1)+NaN;
switch lower(test) 
    case('spearman')
        [pValue,trend,T]=SP(data,'period',period,'missing',missing);
    case('mannkendall')
        [pValue,trend,T]=MK(data,'period',period,'missing',missing,'autocorrelation',autocorrelation);
end
