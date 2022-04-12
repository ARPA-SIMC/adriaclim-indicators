function [dataM,dateM]=aggregateData(data,datDate,grpBy,varargin)
%[dataM,dateM]=aggregateData(data,datDate,grpBy,varargin)
%   DATA Data to be aggregate by rows. "datDate" are the serial dates (obtained from datenum) 
%   associated to each row of "data". GROUPBY is the agrupation criteria, for annual write
%   Y, for monthly M and for seasonally S (S1,S2,S3,S4 correspond to
%   DJF,MAM,JJA and SON). 
%   
%   OPTIONS 
%        'aggFun': Aggragating funtion, by default 'nanmean'. Alternatives are 
%                  nanmin, nansum, sum, min, max,nanstd, and user defined
%        'aggDate': Optional aggregating dates
%        'missing': maximun percentage of missing data (i.e. NaN) per group   
%        'step'   : aggregation block ['Y','S','M',{'D'}]
%   Dates are can be also introduced as a cell with 2 dates
%   {beginDate,endDate}. In this case daily step is considered, and the
%   serial date is build
aggFun='nanmean';
missing=0.1;
aggDate=datDate;
step='D';
for i=1:2:length(varargin)
    switch lower(varargin{i}),
		case 'aggfun',  aggFun = varargin{i+1};
		case 'missing', missing = varargin{i+1};
		case 'aggdate', aggDate = varargin{i+1};
		case 'step', step = varargin{i+1};
		otherwise, error('Unknown optional argument: %s',varargin{i});
    end
end
dataM=[];
dateM=[];
if(~isnumeric(datDate) && ~iscell(datDate) && length(datDate)~=2)
    error('Incorrect data date format');
end
if(~isempty(aggDate) && ~isnumeric(datDate) && ~iscell(datDate) && length(datDate)~=2)
    error('Incorrect aggregation date format');
end

if(iscell(datDate))
	switch step
		case {'Y','S','M','D','24:00'}
			datDate=datenum(datDate{1}):datenum(datDate{2});
		otherwise
			datDate=datenum(datDate{1}):datenum(stepvec(step)):datenum(datDate{2});
	end
end
if(iscell(aggDate))
	switch step
		case {'Y','S','M','D','24:00'}
			aggDate=datenum(aggDate{1}):datenum(aggDate{2});
		otherwise
			aggDate=datenum(aggDate{1}):datenum(stepvec(step)):datenum(aggDate{2});
	end
end
fechaIni=datevec(aggDate(1));
fechaFin=datevec(aggDate(end));
switch grpBy       
	case 'Y'
		format1='yyyy';
		startDate=datenum([fechaIni(1) 1 1 0 0 0]);
		endDate=datenum([fechaFin(1) 12 31 23 0 0]);
	case 'S'
		format1='YYYYSS';
		monthIni=[12 3 6 9];
		monthFin=[2 5 8 11];
		suma=0;
		if ismember(fechaIni(2),[12 1 2])
			season=1;if fechaIni(2)~=12,suma=-1;end
		elseif ismember(fechaIni(2),[3 4 5])
			season=2;
		elseif ismember(fechaIni(2),[6 7 8])
			season=3;
		else
			season=4;
		end
		startDate=datenum([fechaIni(1)+suma monthIni(season) 1 0 0 0]);
		suma=0;
		if ismember(fechaFin(2),[12 1 2])
			season=1;if fechaFin(2)==12,suma=1;end
		elseif ismember(fechaFin(2),[3 4 5])
			season=2;
		elseif ismember(fechaFin(2),[6 7 8])
			season=3;
		else
			season=4;
		end
		endDate=datenum([fechaFin(1)+suma monthFin(season) 1 0 0 0]);
	case 'M'
		format1='yyyymm';
		startDate=datenum([fechaIni(1:2) 1 0 0 0]);
		endDate=datenum([fechaFin(1) fechaFin(2)+1 1 0 0 0])-datenum([0 0 0 1 0 0]);
	case {'D','24:00'}
		format1='yyyymmdd';
		startDate=datenum([fechaIni(1:3) 0 0 0]);
		endDate=datenum([fechaFin(1:3) 23 0 0]);
	case {'12h','12H','12:00'}
		format1='yyyymmdd12';
		startDate=datenum(fechaIni);
		endDate=datenum(fechaFin);
	case {'6h','6H','06:00'}
		format1='yyyymmdd66';
		startDate=datenum(fechaIni);
		endDate=datenum(fechaFin);
	case {'3h','3H','03:00'}
		format1='yyyymmdd33';
		startDate=datenum(fechaIni);
		endDate=datenum(fechaFin);
	case {'1h','1H','H','01:00'}
		format1='yyyymmddhh';				
		startDate=datenum(fechaIni);
		endDate=datenum(fechaFin);
end 
switch step
	case 'Y'
		format='yyyy';
	case 'S'
		format='YYYYSS';
	case 'M'
		format='yyyymm';
	case {'D','24:00'}
		format='yyyymmdd';
	case {'12h','12H','12:00'}
		format='yyyymmdd12';
	case {'6h','6H','06:00'}
		format='yyyymmdd66';
	case {'3h','3H','03:00'}
		format='yyyymmdd33';
	case {'1h','1H','H','01:00'}
		format='yyyymmddhh';
end

fechaIni=datesym(max(datenum(fechaIni),datDate(1)),format);
fechaFin=datesym(min(datenum(fechaFin),datDate(end)),format);
datDateSym=unique(datesym(datDate,format),'rows');
aggDateSym=unique(datesym(aggDate,format),'rows');
[dateM,J1,J2]=datesym2datesym(aggDateSym,step,grpBy);
nq=histc(J1,unique(J1));ndata=length(nq);
fechaIni1=strmatch(fechaIni,aggDateSym);
fechaFin1=strmatch(fechaFin,aggDateSym);
J1=J1(fechaIni1:fechaFin1);I1=unique(J1);
[a,b,c]=intersect(aggDateSym(fechaIni1:fechaFin1,:),datDateSym,'rows');
dataM=NaN*zeros(ndata,size(data,2));
if ~isempty(b)
	for i=1:length(I1),
		ind=intersect(find(J1==I1(i)),b);
        if(~isempty(ind))
            if all(isnan(data(c(ind),:)))
                dataM(I1(i),:)=NaN;
            else
                dataM(I1(i),:)=mifeval(aggFun,data(c(ind),:));
                %we consider NaN if missing is more than TH*N missing
                dataM(I1(i),sum(~isnan(data(c(ind),:)),1)<=(nq(I1(i))*(1-missing)))=NaN;
            end            
        end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=mifeval(fun,x,dim)

if any(strcmp(fun,{'min','max'}))
    y=feval(fun,x,[],1);
else                
    y=feval(fun,x);
end
