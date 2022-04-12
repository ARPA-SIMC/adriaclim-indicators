function S = datesym_dev(D,dateform,pivotyear)
%DATESYM String representation of date.
%   DATESYM number   DATEFORM string         Example
%      0             'dd-mmm-yyyy HH:MM:SS'   01-Mar-2000 15:45:17 
%      1             'dd-mmm-yyyy'            01-Mar-2000  
%      2             'mm/dd/yy'               03/01/00     
%      3             'mmm'                    Mar          
%      4             'm'                      M            
%      5             'mm'                     03            
%      6             'mm/dd'                  03/01        
%      7             'dd'                     01            
%      8             'ddd'                    Wed          
%      9             'd'                      W            
%     10             'yyyy'                   2000         
%     11             'yy'                     00           
%     12             'mmmyy'                  Mar00        
%     13             'HH:MM:SS'               15:45:17     
%     14             'HH:MM:SS PM'             3:45:17 PM  
%     15             'HH:MM'                  15:45        
%     16             'HH:MM PM'                3:45 PM     
%     17             'QQ-YY'                  Q1-96        
%     18             'QQ'                     Q1           
%     19             'dd/mm'                  01/03        
%     20             'dd/mm/yy'               01/03/00     
%     21             'mmm.dd,yyyy HH:MM:SS'   Mar.01,2000 15:45:17 
%     22             'mmm.dd,yyyy'            Mar.01,2000  
%     23             'mm/dd/yyyy'             03/01/2000 
%     24             'dd/mm/yyyy'             01/03/2000 
%     25             'yy/mm/dd'               00/03/01 
%     26             'yyyy/mm/dd'             2000/03/01 
%     27             'QQ-YYYY'                Q1-1996        
%     28             'mmmyyyy'                Mar2000        
%     29 (ISO 8601)  'yyyy-mm-dd'             2000-03-01
%     30 (ISO 8601)  'yyyymmddTHHMMSS'        20000301T154517 
%     31             'yyyy-mm-dd HH:MM:SS'    2000-03-01 15:45:17 
%    
%   See also DATE, DATENUM, DATEVEC, DATETICK.

% handle input arguments
error(nargchk(1,3,nargin)) % handle number of input arguments.

% Convert strings and clock vectors to date numbers.
if isstr(D) | iscell(D) | (size(D,2)==6 & ...
	all(all(D(:,1:5)==fix(D(:,1:5)))) & all(abs(sum(D,2)-2000)<500))
	if nargin < 3
		D = datenum(D); 
	else
		D = datenum(D,pivotyear); 
	end
end

% Determine format if none specified.  If all the times are zero,
% display using date only.  If all dates are all zero display time only.
% Otherwise display both time and date.

D = D(:);
if (nargin < 2) | (dateform == -1)
	if all(floor(D)==D),
		dateform = 1;
	elseif all(floor(D)==0)
		dateform = 16;
	else
		dateform = 0;
	end
end

% Handle the empty case properly.  Return an empty which is the same
% length of the string that is normally returned for each dateform.
if isempty(D)
	try
		S= zeros(0,length(datestr(now,dateform)));
	catch
		error(lasterr);
	end
	return;
end

% Determine from string.
if isstr(dateform)
	switch dateform
		case 'dd-mmm-yyyy HH:MM:SS', dateform = 0;
		case 'dd-mmm-yyyy', dateform = 1;
		case 'mm/dd/yy', dateform = 2;
		case 'mmm', dateform = 3;
		case 'm', dateform = 4;
		case 'mm', dateform = 5;
		case 'mm/dd', dateform = 6;
		case 'dd', dateform = 7;
		case 'ddd', dateform = 8;
		case 'd', dateform = 9;
		case 'yyyy', dateform = 10;
		case 'yy', dateform = 11;
		case 'mmmyy', dateform = 12;
		case 'HH:MM:SS', dateform = 13;
		case 'HH:MM:SS PM', dateform = 14;
		case 'HH:MM', dateform = 15;
		case 'HH:MM PM', dateform = 16;
		case 'QQ-YY', dateform = 17;
		case 'QQ', dateform = 18;
		case 'dd/mm', dateform = 19;
		case 'dd/mm/yy', dateform = 20;
		case 'mmm.dd,yyyy HH:MM:SS', dateform = 21;
		case 'mmm.dd,yyyy', dateform = 22;
		case 'mm/dd/yyyy', dateform = 23;
		case 'dd/mm/yyyy', dateform = 24;
		case 'yy/mm/dd', dateform = 25;
		case 'yyyy/mm/dd', dateform = 26;
		case 'QQ-YYYY', dateform = 27;
		case 'mmmyyyy', dateform = 28; 
		case 'yyyy-mm-dd', dateform = 29;
		case 'yyyymmddTHHMMSS', dateform = 30;
		case 'yyyy-mm-dd HH:MM:SS', dateform = 31;
		case 'yyyymm', dateform = 32;
		case 'yyyydd', dateform = 33;
		case 'mmdd', dateform = 34;
		case 'yyyymmdd', dateform = 35;
		case 'YYYYQQ', dateform = 36;
		case 'YYYYSS', dateform = 37;
		case 'yyyymmdd12', dateform = 38;
		case 'yyyymmdd66', dateform = 39;
		case 'yyyymmdd33', dateform = 40;
		case 'yyyymmddhh', dateform = 41;
		otherwise
			error(['Unknown date format: ' dateform])
   end
end

mths = real(['Jan';'Feb';'Mar';'Apr';'May';'Jun';
   'Jul';'Aug';'Sep';'Oct';'Nov';'Dec']);

% Obtain components using mex file, rounding to integer number of seconds.

[y,mo,d,h,min,s] = datevec(D);  mo(mo==0) = 1;

% Vectorization is done within sprintf.  
switch dateform
case 0   % 'dd-mmm-yyyy HH:MM:SS'
	f = '%.2d-%c%c%c-%c%.3d %.2d:%.2d:%.2d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[d,mths(mo,:),c,mod(abs(y),1000),h,min,round(s)]');
case 1   % 'dd-mmm-yyyy'
	f = '%.2d-%c%c%c-%c%.3d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[d,mths(mo,:),c,mod(abs(y),1000)]');
case 2   % 'mm/dd/yy'
	f = '%.2d/%.2d/%.2d';
	S = sprintf(f,[mo,d,mod(abs(y),100)]');
case 3   % 'mmm'
	f = '%c%c%c';
	S = sprintf(f,[mths(mo,:)]');
case 4   % 'm'
	f = '%c';
	S = sprintf(f,[mths(mo,1)]');
case 5   % 'mm'
	f = '%.2d';
	S = sprintf(f,[mo]');
case 6   % 'mm/dd'
	f = '%.2d/%.2d';
	S = sprintf(f,[mo,d]');
case 7   % 'dd'
	f = '%.2d';
	S = sprintf(f,[d]');
case 8   % 'ddd'
	f = '%c%c%c';
	[ignore,w] = weekday(D);
	S = sprintf(f,[w]');
case 9   % 'd'
	f = '%c';
	[ignore,w] = weekday(D);
	S = sprintf(f,[w(:,1)]');
case 10   % 'yyyy'
	f = '%c%.3d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[c,mod(abs(y),1000)]');
case 11   % 'yy'
	f = '%.2d';
	S = sprintf(f,[mod(abs(y),100)]');
case 12   % 'mmmyy'
	f = '%c%c%c%.2d';
	S = sprintf(f,[mths(mo,:),mod(abs(y),100)]');
case 13   % 'HH:MM:SS'
	f = '%.2d:%.2d:%.2d';
	S = sprintf(f,[h,min,round(s)]');
case 14   % 'HH:MM:SS PM'
	f = '%2d:%.2d:%.2d %cM';
	c = h; c(h<12) = 'A'; c(h>=12) = 'P';
	h = mod(h-1,12) + 1;
	S = sprintf(f,[h,min,round(s),c]');
case 15   % 'HH:MM'
	f = '%.2d:%.2d';
	S = sprintf(f,[h,min]');
case 16   % 'HH:MM PM'
	f = '%2d:%.2d %cM';
	c = h; c(h<12) = 'A'; c(h>=12) = 'P';
	h = mod(h-1,12) + 1;
	S = sprintf(f,[h,min,c]');
case 17   % 'QQ-YY'
	q = floor((mo-1)/3)+1;
	f = 'Q%.1d-%.2d';
	S = sprintf(f,[q,mod(abs(y),100)]');
case 18   % 'QQ'
	q = floor((mo-1)/3)+1;
	f = 'Q%.1d';
	S = sprintf(f,[q]');
case 19   % 'dd/mm'
	f = '%.2d/%.2d';
	S = sprintf(f,[d,mo]');
case 20   % 'dd/mm/yy'
	f = '%.2d/%.2d/%.2d';
	S = sprintf(f,[d,mo,mod(abs(y),100)]');
case 21   % 'mmm.dd,yyyy HH:MM:SS'
	f = '%c%c%c.%.2d,%c%.3d %.2d:%.2d:%.2d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[mths(mo,:),d,c,mod(abs(y),1000),h,min,round(s)]');
case 22   % 'mmm.dd,yyyy'
	f = '%c%c%c.%.2d,%c%.3d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[mths(mo,:),d,c,mod(abs(y),1000)]');
case 23   % 'mm/dd/yyyy'
	f = '%.2d/%.2d/%c%.3d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[mo,d,c,mod(abs(y),1000)]');       
case 24   % 'dd/mm/yyyy'
	f = '%.2d/%.2d/%c%.3d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[d,mo,c,mod(abs(y),1000)]'); 
case 25   % 'yy/mm/dd'
	f = '%.2d/%.2d/%.2d';
	S = sprintf(f,[mod(abs(y),100),mo,d]');
case 26   % 'yyyy/mm/dd'
	f = '%c%.3d/%.2d/%.2d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[c,mod(abs(y),1000),mo,d]');  
case 27   % 'QQ-YYYY'
	q = floor((mo-1)/3)+1;
	f = 'Q%.1d-%c%.3d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[q,c,mod(abs(y),1000)]');
case 28   % 'mmmyyyy'
	f = '%c%c%c%c%.3d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[mths(mo,:),c,mod(abs(y),1000)]');      
case 29   % 'yyyy-mm-dd'
	f = '%c%.3d-%.2d-%.2d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[c,mod(abs(y),1000),mo,d]');  
case 30   % 'yyyymmddTHHMMSS'
	f = '%c%.3d%.2d%.2dT%.2d%.2d%.2d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[c,mod(abs(y),1000),mo,d,h,min,round(s)]');
case 31   % 'yyyy-mm-dd HH:MM:SS'
	f = '%c%.3d-%.2d-%.2d %.2d:%.2d:%.2d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[c,mod(abs(y),1000),mo,d,h,min,round(s)]');
case 32   % 'yyyymm'
	f = '%c%.3d%.2d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[c,mod(abs(y),1000),mo]');
case 33   % 'yyyydd'
	f = '%c%.3d%.2d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[c,mod(abs(y),1000),d]');
case 34   %'mmdd'
	f = '%.2d%.2d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[mo,d]');
case 35   %'yyyymmdd'
	f = '%c%.3d%.2d%.2d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[c,mod(abs(y),1000),mo,d]');
case 36   % 'YYYYQQ'
	q = floor((mo-1)/3)+1;
	f = '%c%.3dQ%.1d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[c,mod(abs(y),1000),q]');
case 37   % 'YYYYSS'
	q = floor(rem(mo,12)/3)+1;
	y(mo==12)=y(mo==12)+1;
	%q = floor(q/3)+1;
	f = '%c%.3dS%.1d';
	c = '0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S = sprintf(f,[c,mod(abs(y),1000),q]');
case 38   % 'yyyymmdd12'
	h=div(h,12)*12;
	f='%c%.3d%.2d%.2d%.2d';
	c='0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S=sprintf(f,[c,mod(abs(y),1000),mo,d,h]');
case 39   % 'yyyymmdd66'
	h=div(h,6)*6;
	f='%c%.3d%.2d%.2d%.2d';
	c='0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S=sprintf(f,[c,mod(abs(y),1000),mo,d,h]');
case 40   % 'yyyymmdd33'
	h=div(h,3)*3;
	f='%c%.3d%.2d%.2d%.2d';
	c='0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S=sprintf(f,[c,mod(abs(y),1000),mo,d,h]');
case 41   % 'yyyymmddhh'
	f='%c%.3d%.2d%.2d%.2d';
	c='0'+fix(mod(y,10000)/1000); c(y<0) = '-';
	S=sprintf(f,[c,mod(abs(y),1000),mo,d,h]');
otherwise
   error(['Unknown date format number: ' int2str(dateform)])
end

S=reshape(S',length(S)/length(D),length(D))';
