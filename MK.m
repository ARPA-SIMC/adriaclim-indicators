function [pValue,trend,u,newData]=MK(data,varargin)
% [pVal,u,newData]=MK(data,varargin)   Mann-Kendall test for serial trend
% 
% Input
% 	data(:,k)   : Daily data series for each of the 'k' stations.
% 	varargin	: optional parameters
% 	  'missing' - missing parameter for the 'movingAverage' function
%     'period'  - 'day', 'month' or 'year'
% 	
% Output    
% 	pValue  : p-value of the test for each of the stations.
% 	trend    : 1xNest matrix with serial trend.
%   	u           : test's statistic
%   newData   : Transformed data used in the test

% Examples
% 
% 		pVal=MK(data,'period','year');
%       [pVal,u,dataY]=MK(data,'period','year','missing',0.1);

missing=0; 
period='';
autocorrelation=0;
for i=1:2:length(varargin)
    switch lower(varargin{i}),
        case 'period', period=varargin{i+1};
        case 'missing', missing=varargin{i+1};
        case 'autocorrelation', autocorrelation=varargin{i+1};
    end
end
[ndata,Nest]=size(data);
pValue=NaN*zeros(1,Nest);
trend=NaN*zeros(1,Nest);
u=NaN*zeros(1,Nest);
newData=[];
k=1;
switch(period)
    case('day')
        k=1;
    case('month')
        k=30;
    case('year')
        k=365;
end
newData=movingAverage(data,k,'missing',missing);
stations=find(~isnan(nansum(data,1)));

for i=stations,%1:Nest
    dato=newData(find(isnan(newData(:,i))==0),i);
	b=regress(dato,[ones(length(dato),1) [1:length(dato)]']);
	trend(i)=b(2);
    N=length(dato);
    t=0;
    for l=1:N
        for j=l+1:N
            t=t+sign(dato(j)-dato(l));
        end
    end
    var=(1/18)*N*(N-1)*(2*N+5);
    if length(unique(dato))<N
        aux1=unique(dato);
        for l=1:length(aux1)
            a=length(find(dato==aux1(l)));
            if a>1
                var=var-(1/18)*a*(a-1)*(2*a+5);
            end
        end
    end
	switch autocorrelation
		case {1,'true','yes'}
			% Consideramos la autocorrelacion de la serie:
			acf=repmat(NaN,N-1,2);rankeData=crank(dato);
			for l=1:N-1,
				[a1,a2]=corrcoef(rankeData(1:N-l),rankeData(l+1:N));
				acf(l,:)=[a1(1,size(a1,2)),a2(1,size(a2,2))];
			end
			acf(acf(:,2)>0.05,1)=0;
			ess=1+(2/(N*(N-1)*(N-2)))*nansum(((N-[1:N-1]').*(N-1-[1:N-1]').*(N-2-[1:N-1]')).*abs(acf(:,1)));
			var=var*ess;
	end
    if t==0
        if var~=0
            u(i)=0;
        else
            pValue(i)=0;trend(i)=0;
        end
    elseif t>0
        u(i)=(t-1)/sqrt(var);
    elseif t<0
        u(i)=(t+1)/sqrt(var);
    end
    if u(i)>0
        pValue(i)=2*(1-normcdf(u(i)));
    end
    if u(i)<=0
        pValue(i)=2*normcdf(u(i));
    end
end

function [r,ra]=crank(x)

u=unique(x);
[xs,z1]=sort(x);
[z1,z2]=sort(z1);
r=(1:length(x))';
r=r(z2);
ra=0;
for i=1:length(u)
	s=find(u(i)==x);
	r(s)=nanmean(r(s));
	ra=ra+length(s)*(length(s)-1)*(length(s)+1)/2;
end
