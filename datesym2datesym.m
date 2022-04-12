function [dateList,I1,I2]=datesym2datesym(dateList,format1,format2)

if strcmp(format1,format2)
	I1=[1:size(dateList,1)]';
	I2=I1;
	dateList=dateList;
else
	switch format1
		case 'Y'
			error('The format2 must be greater than the format1')
		case 'S'
			if strcmp(format2,'Y')
				[dateList,I2,I1]=unique(dateList(:,1:4),'rows');
			else
				error('The format2 must be greater than the format1')
			end
		case 'M'
			switch format2
				case 'Y'
					[dateList,I2,I1]=unique(dateList(:,1:4),'rows');
				case 'S'
					meses=str2num(dateList(:,5:6));
					ind=find(meses==12);
					if ~isempty(ind)
						aux=str2num(dateList(ind,1:4));
						aux=aux+1;dateList(ind,1:4)=num2str(aux);
					end
					months=[12 1 2;3 4 5;6 7 8;9 10 11];
					meses=str2num(dateList(:,5:6));
					for i=1:4
						ind=find(ismember(meses,months(i,:)));
						if ~isempty(ind)
							dateList(ind,5:6)=repmat(['S' num2str(i)],length(ind),1);
						end
					end
					[dateList,I2,I1]=unique(dateList,'rows');
				otherwise
					error('The format2 must be greater than the format1')
			end
		case {'D','24:00'}
			switch format2
				case 'Y'
					[dateList,I2,I1]=unique(dateList(:,1:4),'rows');
				case 'S'
					meses=str2num(dateList(:,5:6));
					ind=find(meses==12);
					if ~isempty(ind)
						aux=str2num(dateList(ind,1:4));
						aux=aux+1;dateList(ind,1:4)=num2str(aux);
					end
					months=[12 1 2;3 4 5;6 7 8;9 10 11];
					meses=str2num(dateList(:,5:6));
					for i=1:4
						ind=find(ismember(meses,months(i,:)));
						if ~isempty(ind)
							dateList(ind,5:6)=repmat(['S' num2str(i)],length(ind),1);
						end
					end
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case 'M'
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
                otherwise
                    if strmatch(format2(end),'D') & size(format2,2)>1
                        I2=[1:str2num(format2(1:end-1)):size(dateList,1)]';
                        I1=kron([1:length(I2)]',ones(str2num(format2(1:end-1)),1));
                        dateList=[dateList(I2,1:6) repmat('W',length(I2),1) num2str([1:length(I2)]')];
                    else
                        error('The format2 must be greater than the format1')
                    end
			end
		case {'12h','12H','12:00'}
			switch format2
				case 'Y'
					[dateList,I2,I1]=unique(dateList(:,1:4),'rows');
				case 'S'
					meses=str2num(dateList(:,5:6));
					ind=find(meses==12);
					if ~isempty(ind)
						aux=str2num(dateList(ind,1:4));
						aux=aux+1;dateList(ind,1:4)=num2str(aux);
					end
					months=[12 1 2;3 4 5;6 7 8;9 10 11];
					meses=str2num(dateList(:,5:6));
					for i=1:4
						ind=find(ismember(meses,months(i,:)));
						if ~isempty(ind)
							dateList(ind,5:6)=repmat(['S' num2str(i)],length(ind),1);
						end
					end
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case 'M'
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case {'D','24:00'}
					[dateList,I2,I1]=unique(dateList(:,1:8),'rows');
				otherwise
					error('The format2 must be greater than the format1')
			end
		case {'6h','6H','06:00'}
			switch format2
				case 'Y'
					[dateList,I2,I1]=unique(dateList(:,1:4),'rows');
				case 'S'
					meses=str2num(dateList(:,5:6));
					ind=find(meses==12);
					if ~isempty(ind)
						aux=str2num(dateList(ind,1:4));
						aux=aux+1;dateList(ind,1:4)=num2str(aux);
					end
					months=[12 1 2;3 4 5;6 7 8;9 10 11];
					meses=str2num(dateList(:,5:6));
					for i=1:4
						ind=find(ismember(meses,months(i,:)));
						if ~isempty(ind)
							dateList(ind,5:6)=repmat(['S' num2str(i)],length(ind),1);
						end
					end
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case 'M'
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case {'D','24:00'}
					[dateList,I2,I1]=unique(dateList(:,1:8),'rows');
				case {'12h','12H','12:00'}
					horas=str2num(dateList(:,9:10));
					i00=find(horas<12);i12=find(horas>=12);
					dateList(i00,9:10)=repmat('00',length(i00),1);
					dateList(i12,9:10)=repmat('12',length(i00),1);
					[dateList,I2,I1]=unique(dateList(:,1:10),'rows');
				otherwise
					error('The format2 must be greater than the format1')
			end
		case {'3h','3H','03:00'}
			switch format2
				case 'Y'
					[dateList,I2,I1]=unique(dateList(:,1:4),'rows');
				case 'S'
					meses=str2num(dateList(:,5:6));
					ind=find(meses==12);
					if ~isempty(ind)
						aux=str2num(dateList(ind,1:4));
						aux=aux+1;dateList(ind,1:4)=num2str(aux);
					end
					months=[12 1 2;3 4 5;6 7 8;9 10 11];
					meses=str2num(dateList(:,5:6));
					for i=1:4
						ind=find(ismember(meses,months(i,:)));
						if ~isempty(ind)
							dateList(ind,5:6)=repmat(['S' num2str(i)],length(ind),1);
						end
					end
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case 'M'
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case {'D','24:00'}
					[dateList,I2,I1]=unique(dateList(:,1:8),'rows');
				case {'12h','12H','12:00'}
					horas=str2num(dateList(:,9:10));
					i00=find(horas<12);i12=find(horas>=12);
					dateList(i00,9:10)=repmat('00',length(i00),1);
					dateList(i12,9:10)=repmat('12',length(i00),1);
					[dateList,I2,I1]=unique(dateList(:,1:10),'rows');
				case {'6h','6H','06:00'}
					horas=str2num(dateList(:,9:10));
					ind1=[];
					for i=1:4
						ind=find(horas<i*6);
						ind=setdiff(ind,ind1);ind1=union(ind,ind1);
						sd=num2str(100+(i-1)*6);
						if ~isempty(ind)
							dateList(ind,9:10)=repmat(sd(2:3),length(ind),1);
						end
					end
					[dateList,I2,I1]=unique(dateList(:,1:10),'rows');
				otherwise
					error('The format2 must be greater than the format1')
			end
		case {'1h','1H','H','01:00'}
			switch format2
				case 'Y'
					[dateList,I2,I1]=unique(dateList(:,1:4),'rows');
				case 'S'
					meses=str2num(dateList(:,5:6));
					ind=find(meses==12);
					if ~isempty(ind)
						aux=str2num(dateList(ind,1:4));
						aux=aux+1;dateList(ind,1:4)=num2str(aux);
					end
					months=[12 1 2;3 4 5;6 7 8;9 10 11];
					meses=str2num(dateList(:,5:6));
					for i=1:4
						ind=find(ismember(meses,months(i,:)));
						if ~isempty(ind)
							dateList(ind,5:6)=repmat(['S' num2str(i)],length(ind),1);
						end
					end
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case 'M'
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case {'D','24:00'}
					[dateList,I2,I1]=unique(dateList(:,1:8),'rows');
				case {'12h','12H','12:00'}
					horas=str2num(dateList(:,9:10));
					i00=find(horas<12);i12=find(horas>=12)
					dateList(i00,9:10)=repmat('00',length(i00),1);
					dateList(i12,9:10)=repmat('12',length(i00),1);
					[dateList,I2,I1]=unique(dateList(:,1:10),'rows');
				case {'6h','6H','06:00'}
					horas=str2num(dateList(:,9:10));
					ind1=[];
					for i=1:4
						ind=find(horas<i*6);
						ind=setdiff(ind,ind1);ind1=union(ind,ind1);
						sd=num2str(100+(i-1)*6);
						if ~isempty(ind)
							dateList(ind,9:10)=repmat(sd(2:3),length(ind),1);
						end
					end
					[dateList,I2,I1]=unique(dateList(:,1:10),'rows');
				case {'3h','3H','03:00'}
					horas=str2num(dateList(:,9:10));
					ind1=[];
					for i=1:8
						ind=find(horas<i*3);
						ind=setdiff(ind,ind1);ind1=union(ind,ind1);
						sd=num2str(100+(i-1)*3);
						if ~isempty(ind)
							dateList(ind,9:10)=repmat(sd(2:3),length(ind),1);
						end
					end
					[dateList,I2,I1]=unique(dateList(:,1:10),'rows');
				otherwise
					error('The format2 must be greater than the format1')
			end
		case {'standard'}
			switch format2
				case 'Y'
					[dateList,I2,I1]=unique(dateList(:,1:4),'rows');
				case 'S'
					meses=str2num(dateList(:,5:6));
					ind=find(meses==12);
					if ~isempty(ind)
						aux=str2num(dateList(ind,1:4));
						aux=aux+1;dateList(ind,1:4)=num2str(aux);
					end
					months=[12 1 2;3 4 5;6 7 8;9 10 11];
					meses=str2num(dateList(:,5:6));
					for i=1:4
						ind=find(ismember(meses,months(i,:)));
						if ~isempty(ind)
							dateList(ind,5:6)=repmat(['S' num2str(i)],length(ind),1);
						end
					end
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case 'M'
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case '360_day'
					x=dateList;
					dayIni=str2num(x(1,7:8));
					dayFin=str2num(x(end,7:8));
					month=str2num(unique(x(:,1:6),'rows'));
                    indices=1:min(30*length(month),size(dateList,1));
					dateList=dateList(indices,:);
                    aux=num2str(kron(month,ones(30,1)));
					dateList(indices,1:6)=aux(indices,:);
					aux=num2str([1:30]','%02d');
                    aux=repmat(aux,length(month),1);
                    dateList(:,7:8)=aux(1:size(dateList,1),:);
					if dayIni<=30
						dateList=dateList(dayIni:end,:);
					else
						dateList=dateList(31:end,:);
					end
					if dayFin<=30
						dateList=dateList(1:end-(30-dayFin),:);
					end
					[a1,I1,I2]=intersect(x,dateList,'rows');
				case 'no_leap'
					x=dateList;
					yearIni=str2num(x(1,1:4));
					yearFin=str2num(x(end,1:4));
					fechas=zeros((yearFin-yearIni+1)*365,6);
					fechas(:,1)=kron([yearIni:yearFin]',ones(365,1));						
					fechas(:,2)=repmat([ones(31,1);2*ones(28,1);3*ones(31,1);4*ones(30,1);5*ones(31,1);6*ones(30,1);7*ones(31,1);8*ones(31,1);9*ones(30,1);10*ones(31,1);11*ones(30,1);12*ones(31,1)],yearFin-yearIni+1,1);
					fechas(:,3)=repmat([[1:31]';[1:28]';[1:31]';[1:30]';[1:31]';[1:30]';[1:31]';[1:31]';[1:30]';[1:31]';[1:30]';[1:31]'],yearFin-yearIni+1,1);
					dayIni=[str2num(x(1,1:4)) str2num(x(1,5:6)) str2num(x(1,7:8)) 0 0 0];
					dayFin=[str2num(x(end,1:4)) str2num(x(end,5:6)) str2num(x(end,7:8)) 0 0 0];
					if isequal(dayIni,[dayIni(1) 2 29 0 0 0])
						dayIni=[dayIni(1) 3 1 0 0 0];
					end
					if isequal(dayFin,[dayFin(1) 2 29 0 0 0])
						dayFin=[dayFin(1) 3 1 0 0 0];
					end
					[a1,indIni,a3]=intersect(fechas,dayIni,'rows');
					[a1,indFin,a3]=intersect(fechas,dayFin,'rows');
					fechas=datenum(fechas(indIni:indFin,:));
					dateList=unique(datesym(fechas,'yyyymmdd'),'rows');
					[a1,I1,I2]=intersect(x(:,1:8),dateList(:,1:8),'rows');
				otherwise
					error('The format2 must be greater than the format1')
			end
		case {'no_leap'}
			switch format2
				case 'Y'
					[dateList,I2,I1]=unique(dateList(:,1:4),'rows');
				case 'S'
					meses=str2num(dateList(:,5:6));
					ind=find(meses==12);
					if ~isempty(ind)
						aux=str2num(dateList(ind,1:4));
						aux=aux+1;dateList(ind,1:4)=num2str(aux);
					end
					months=[12 1 2;3 4 5;6 7 8;9 10 11];
					meses=str2num(dateList(:,5:6));
					for i=1:4
						ind=find(ismember(meses,months(i,:)));
						if ~isempty(ind)
							dateList(ind,5:6)=repmat(['S' num2str(i)],length(ind),1);
						end
					end
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case 'M'
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case '360_day'
					x=dateList(:,1:8);
					dayIni=str2num(x(1,7:8));
					dayFin=str2num(x(end,7:8));
					month=str2num(unique(x(:,1:6),'rows'));
                    indices=1:min(30*length(month),size(dateList,1));
					dateList=dateList(indices,:);
                    aux=num2str(kron(month,ones(30,1)));
					dateList(indices,1:6)=aux(indices,:);
					aux=num2str([1:30]','%02d');
                    aux=repmat(aux,length(month),1);
                    dateList(:,7:8)=aux(1:size(dateList,1),:);
					if dayIni<=30
						dateList=dateList(dayIni:end,:);
					else
						dateList=dateList(31:end,:);
					end
					if dayFin<=30
						dateList=dateList(1:end-(30-dayFin),:);
					end
					[a1,I1,I2]=intersect(x(:,1:8),dateList(:,1:8),'rows');
				case 'standard'
					x=dateList;
                    if size(x,2)>8,
                        h1=str2num(x(1,9:10));h2=str2num(x(end,9:10));hh='yyyymmddhh';h=datenum([0 0 0 1 0 0]);
                        dayIni=datenum([str2num(x(1,1:4)) str2num(x(1,5:6)) str2num(x(1,7:8)) h1 0 0]);
                        dayFin=datenum([str2num(x(end,1:4)) str2num(x(end,5:6)) str2num(x(end,7:8)) h2 0 0]);
                    else,
                        h1=0;h2=0;hh='yyyymmdd';
                        dayIni=datenum([str2num(x(1,1:4)) str2num(x(1,5:6)) str2num(x(1,7:8)) 0 0 0]);h=datenum([0 0 1 0 0 0]);
                        dayFin=datenum([str2num(x(end,1:4)) str2num(x(end,5:6)) str2num(x(end,7:8)) 0 0 0]);
                    end
                    dateList=datesym([dayIni:h:dayFin]',hh);
                    [a1,I1,I2]=intersect(x,dateList,'rows');
% 				case 'standard'
%                     x=dateList;
% 					if size(x,2)>8,h=str2num(x(1,9:10));hh='yyyymmddhh';else,h=0;hh='yyyymmdd';end
% 					dayIni=datenum([str2num(x(1,1:4)) str2num(x(1,5:6)) str2num(x(1,7:8)) h 0 0]);
% 					dayFin=datenum([str2num(x(end,1:4)) str2num(x(end,5:6)) str2num(x(end,7:8)) h 0 0]);
% 					dateList=datesym([dayIni:dayFin]',hh);
% 					[a1,I1,I2]=intersect(x(:,1:8),dateList(:,1:8),'rows');
				otherwise
					error('The format2 must be greater than the format1')
			end
		case {'360_day'}
			switch format2
				case 'Y'
					[dateList,I2,I1]=unique(dateList(:,1:4),'rows');
				case 'S'
					meses=str2num(dateList(:,5:6));
					ind=find(meses==12);
					if ~isempty(ind)
						aux=str2num(dateList(ind,1:4));
						aux=aux+1;dateList(ind,1:4)=num2str(aux);
					end
					months=[12 1 2;3 4 5;6 7 8;9 10 11];
					meses=str2num(dateList(:,5:6));
					for i=1:4
						ind=find(ismember(meses,months(i,:)));
						if ~isempty(ind)
							dateList(ind,5:6)=repmat(['S' num2str(i)],length(ind),1);
						end
					end
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case 'M'
					[dateList,I2,I1]=unique(dateList(:,1:6),'rows');
				case 'no_leap'
					x=dateList(:,1:8);
					yearIni=str2num(x(1,1:4));
					yearFin=str2num(x(end,1:4));
					fechas=zeros((yearFin-yearIni+1)*365,6);
					fechas(:,1)=kron([yearIni:yearFin]',ones(365,1));						
					fechas(:,2)=repmat([ones(31,1);2*ones(28,1);3*ones(31,1);4*ones(30,1);5*ones(31,1);6*ones(30,1);7*ones(31,1);8*ones(31,1);9*ones(30,1);10*ones(31,1);11*ones(30,1);12*ones(31,1)],yearFin-yearIni+1,1);
					fechas(:,3)=repmat([[1:31]';[1:28]';[1:31]';[1:30]';[1:31]';[1:30]';[1:31]';[1:31]';[1:30]';[1:31]';[1:30]';[1:31]'],yearFin-yearIni+1,1);
					dayIni=[str2num(x(1,1:4)) str2num(x(1,5:6)) str2num(x(1,7:8)) 0 0 0];
					dayFin=[str2num(x(end,1:4)) str2num(x(end,5:6)) str2num(x(end,7:8)) 0 0 0];
					if ismember(dayIni,[dayIni(1) 2 29 0 0 0;dayIni(1) 2 30 0 0 0])
						dayIni=[dayIni(1) 3 1 0 0 0];
					end
					if ismember(dayFin,[dayFin(1) 2 29 0 0 0;dayFin(1) 2 30 0 0 0])
						dayFin=[dayFin(1) 3 1 0 0 0];
					end
					[a1,indIni,a3]=intersect(fechas,dayIni,'rows');
					[a1,indFin,a3]=intersect(fechas,dayFin,'rows');
					fechas=datenum(fechas(indIni:indFin,:));
					dateList=unique(datesym(fechas,'yyyymmdd'),'rows');
					[a1,I1,I2]=intersect(x(:,1:8),dateList(:,1:8),'rows');
				case 'standard'
					x=dateList;
                    if size(x,2)>8,
                        h1=str2num(x(1,9:10));h2=str2num(x(end,9:10));hh='yyyymmddhh';h=datenum([0 0 0 1 0 0]);
                        dayIni=datenum([str2num(x(1,1:4)) str2num(x(1,5:6)) str2num(x(1,7:8)) h1 0 0]);
                        dayFin=datenum([str2num(x(end,1:4)) str2num(x(end,5:6)) str2num(x(end,7:8)) h2 0 0]);
                    else,
                        h1=0;h2=0;hh='yyyymmdd';
                        dayIni=datenum([str2num(x(1,1:4)) str2num(x(1,5:6)) str2num(x(1,7:8)) 0 0 0]);h=datenum([0 0 1 0 0 0]);
                        dayFin=datenum([str2num(x(end,1:4)) str2num(x(end,5:6)) str2num(x(end,7:8)) 0 0 0]);
                    end
                    if isequal(dayIni,datenum([str2num(x(1,1:4)) 2 30 h1 0 0]))
                        dayIni=dayIni+1;
                    end
                    if isequal(dayFin,datenum([str2num(x(end,1:4)) 2 30 h2 0 0]))
                        dayFin=dayFin-1;
                    end
                    dateList=datesym([dayIni:h:dayFin]',hh);
                    [a1,I1,I2]=intersect(x,dateList,'rows');
% 					if size(x,2)>8,h=str2num(x(1,9:10));hh='yyyymmddhh';else,h=0;hh='yyyymmdd';end
% 					dayIni=datenum([str2num(x(1,1:4)) str2num(x(1,5:6)) str2num(x(1,7:8)) h 0 0]);
% 					dayFin=datenum([str2num(x(end,1:4)) str2num(x(end,5:6)) str2num(x(end,7:8)) h 0 0]);
% 					if isequal(dayIni,datenum([str2num(x(1,1:4)) 2 30 h 0 0]))
% 						dayIni=dayIni+1;
% 					end
% 					if isequal(dayFin,datenum([str2num(x(end,1:4)) 2 30 h 0 0]))
% 						dayFin=dayFin-1;
% 					end
% 					dateList=datesym([dayIni:dayFin]',hh);
% 					[a1,I1,I2]=intersect(x,dateList,'rows');
				otherwise
					error('The format2 must be greater than the format1')
			end
	end
end
