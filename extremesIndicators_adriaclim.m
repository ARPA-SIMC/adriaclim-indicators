function indicator=extremesIndicators_adriaclim(missing,aggregation,Tg,Tn,Tx,Pr,Hum,Ssa,Nv,Ws,Wsmax,Wd,Tg_th,Tn_th,Tx_th,Pr_th,period,period_th,nombres)

% This function estimates the extreme index defined by the ETCCDI.
% The input are:
% - outputpath: folder path for output
% - modelName: name of the model
% - varagin: optional inputs.
% - variables: {'Tx';'Tn';'Tg';'Pr';'Hum';'Ssa';'Nv';'Ws';'Wd';'Wsmax'}, ndata x Nest dimensions matrix with the daily data. Each row represents an observed day and
% each column a station or grid point. The data units must be: ºC for temperature (Tx,Tn,Tg), mm for precipitation (Pr), % for relative humidity (Hum), mm for surface snow amount (Ssa), mm for snowfall flux (Nv), m/s for wind speed (Ws) and degrees between (0º-360º) for wind direction (Wd).
% - aggregation: is the aggrupation criteria, for annual write Y, for monthly M and for seasonally S (S1,S2,S3,S4 correspond to DJF,MAM,JJA and SON).
% - missing: maximun percentage of missing data (i.e. NaN) per group
% - names: cell with the index names. 
% %%%%%%%%%% Indicators defined in: https://www.ecad.eu/indicesextremes/indicesdictionary.php %%%%%%%%%%%
% Note: The symbol * indicates "calculated for a 5 day window centred on each calendar day in the 1971-2000 period or the period on which the threshold is calculated".   !!!!!!!!!! PERIOD !!!!!
% Note: The symbol ** indicates "calculated for a 31 day window centred on each calendar day in the 1971-2000 period or the period on which the threshold is calculated".   !!!!!!!!!! PERIOD !!!!!

%Indicators:
% 'wsdi': Warm-spell days (days) - Number of days per period where, in intervals of at least 6
% 								 consecutive days Tx>90th percentile* of max temp.
% 'hwn': number of heat waves - Number of times per period where, in intervals of at least 3
% 								 consecutive days Tx>90th percentile** of max temp.
% 'hwtxdx': (°C) - maximum value between the averages of the maximum temperatures averaged for each
%             heat wave. The heat wave is identified by the exceeding, for at least 3 consecutive days,
%             of the 90th percentile** of the maximum temperatures.
% 'txx': Maximun value of daily maximum temperature (C).
% 'txn': Minimun value of daily maximum temperature (°C).
% 'tg':  Mean of daily mean temperature (°C).
% 'tr': Tropical nights (days) - Number of days with Tn>20°C.
% 'tx90p': Warm day-times (days) - Days with Tx>90th percentile* of daily max temp.
% 'fg': Mean of daily wind speed (m/s).
% 'prcptot': Precipitation sum with Pr>1mm (mm).
% 'rx1day': Max 1-day precipitation amount (mm) - Maximun daily value.
% 'rx5day': Max 5-day precipitation amount (mm) - Maximun 5-days value.
% 'r95ptot': Precipitation fraction due to very wet days (%).
% 'cdd': Consecutive dry days (days) - Largest number of consecutive days with Pr<1mm.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO LOAD OCTAVE PACKAGES NEEDED
pkg load statistics
pkg load io
pkg load netcdf
% This function calls the following scripts: 
% aggregateData.m , aggregateDataDay.m and aggregateDataDay_M.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Identify the reference period in order to compute the threshold values
if ~isempty(period_th)
    % vettore con tutte le date del periodo considerato in formato matlab
    period_txt=datestr(period_th,30); %tutte le date del periodo considerato in formato testo
    period_days=[1:length(period_th)];
    [days5Dw,I1,I2]=unique(period_txt(period_days,5:8),'rows'); %period_txt(period_days,5:8) prende la parte della stringa della data contenete mese e giorno !!!!!!!
    ndays5Dw=size(days5Dw,1);
else
    period_th=[];
end

%Number of indicators to be computed
Nindex=length(nombres);
if Nindex<1
    error('At least an indicator name is necessary')
end

%Aggregation Type: Annual, Monthly, Seasonal 
datesAgg=datevec(period(:));
switch aggregation
    case 'Y'
        years=cellstr(num2str(unique(datesAgg(:,1))));
        [a1,a2,a3]=unique(datesAgg(:,1),'rows');
        datesAgg=a3;Nagg=length(a1);clear a1 a2 a3
    case 'S'
        years_cell=cellstr(num2str(unique(datesAgg(:,1))));
        nseas=length(years_cell)*4+1;
        years=repmat({'NaN'},nseas,1);
        seas_cell={'DJF '; 'MAM '; 'JJA '; 'SON '};
        for iyear=1:length(years_cell)
            years((4*(iyear-1)+1):(4*iyear),1)=strcat(seas_cell,repmat(years_cell(iyear),4,1));
        end
        estaciones=[12 1 2;3 4 5;6 7 8;9 10 11];
        aux=find(datesAgg(:,2)<=2);datesAgg(aux,1)=datesAgg(aux,1)-1;datesAgg(aux,2)=1;
        aux=find(datesAgg(:,2)==12);datesAgg(aux,2)=1;
        for i=2:4,datesAgg(find(ismember(datesAgg(:,2),estaciones(i,:))),2)=i;end
        [a1,a2,a3]=unique(datesAgg(:,1:2),'rows');
        datesAgg=a3;Nagg=length(a1);clear a1 a2 a3
    case 'M'
        [a1,a2,a3]=unique(datesAgg(:,1:2),'rows');
        datesAgg=a3;Nagg=length(a1);clear a1 a2 a3
end


%The Indicators are computed below:
for i=1:Nindex

    indicator(i).Name=nombres{i};

    switch lower(nombres{i})

        case {'tg' 'bio1'}
            indicator(i).Index=aggregateData(Tg,period,aggregation,'aggfun','nanmean','missing',missing);
        case 'txn'
            indicator(i).Index=aggregateData(Tx,period,aggregation,'aggfun','nanmin','missing',missing);
        case {'txx' 'bio5'}
            indicator(i).Index=aggregateData(Tx,period,aggregation,'aggfun','nanmax','missing',missing);
        case 'tx90p'
            [ndata,Nest]=size(Tx);
            
            % The threshold values are computed
            threshold.tx90p=NaN(ndays5Dw,Nest);
            indices=[ndays5Dw-1 ndays5Dw 1:ndays5Dw 1 2];
            for indays5Dw=1:ndays5Dw
                ind=find(ismember(I2,indices(indays5Dw:indays5Dw+4)));
                threshold.tx90p(indays5Dw,:)=prctile(Tx_th(ind,:),90);
            end
            clear indices indays5Dw ind
            %The codes extract the epoch where to calculate the thresholds
            indicator(i).Index=NaN(Nagg,Nest);
            fechas=datesym(period,'yyyymmdd');
            diasRef=datesym(period_th(1):period_th(end),'yyyymmdd');
            diasRef=unique(diasRef(:,5:end),'rows');
            fullaux=NaN(ndata,Nest);
            %Days and Number of Days that satisfy the threshold condition
            for k=1:Nest
                aux=NaN(ndata,1);
                for j=1:size(diasRef,1)
                    ind=strmatch(diasRef(j,:),fechas(:,5:end));
                    ind1=find(~isnan(Tx(ind,k)) & Tx(ind,k)>threshold.tx90p(j,k));
                    ind2=find(~isnan(Tx(ind,k)) & Tx(ind,k)<=threshold.tx90p(j,k));
                    aux(ind(ind1))=1;
                    aux(ind(ind2))=0;
                end
                fullaux(:,k)=aux;
            end
            indicator(i).Index=aggregateData(fullaux,period,aggregation,'aggfun','nansum','missing',missing); % Number of days per aggregation (yearly/monthly/seasonal) exceeding the threshold
        case 'wsdi'
            [ndata,Nest]=size(Tx);

            threshold.tx90p=NaN(ndays5Dw,Nest);
            indices=[ndays5Dw-1 ndays5Dw 1:ndays5Dw 1 2];
            for indays5Dw=1:ndays5Dw
                ind=find(ismember(I2,indices(indays5Dw:indays5Dw+4)));
                threshold.tx90p(indays5Dw,:)=prctile(Tx_th(ind,:),90);
            end
            clear indices indays5Dw ind
            %The codes extract the epoch where to calculate the thresholds
            indicator(i).Index=NaN(Nagg,Nest);
            fechas=datesym(period,'yyyymmdd');
            diasRef=datesym(period_th(1):period_th(end),'yyyymmdd');
            diasRef=unique(diasRef(:,5:end),'rows');
            %Days and Number of Days that satisfy the threshold condition   
            for k=1:Nest
                aux=NaN(ndata,1);
                for j=1:size(diasRef,1)
                    ind=strmatch(diasRef(j,:),fechas(:,5:end));
                    ind1=find(~isnan(Tx(ind,k)) & Tx(ind,k)>threshold.tx90p(j,k));
                    ind2=find(~isnan(Tx(ind,k)) & Tx(ind,k)<=threshold.tx90p(j,k));
                    aux(ind(ind1))=1;
                    aux(ind(ind2))=0;
                end
                for j=1:Nagg
                    ind=find(datesAgg==j)';
                    % At least 6 consecutive days must exceed the threshold
                    indicator(i).Index(j,k)=0;ndays=0;
                    for l=ind
                        if aux(l)==1
                            ndays=ndays+1;
                        else
                            if ndays>5
                                indicator(i).Index(j,k)=indicator(i).Index(j,k)+ndays;
                            end
                            ndays=0;
                        end
                    end
                    if all(isnan(Tx(ind))) 
                        indicator(i).Index(j,k)=NaN;
                    end
                end
            end
        case {'hwn'}
            [ndata,Nest]=size(Tx);

            threshold.tx90p=NaN(ndays5Dw,Nest);
            indices=[ndays5Dw-14:ndays5Dw 1:ndays5Dw 1:15];
            for indays5Dw=1:ndays5Dw
                ind=find(ismember(I2,indices(indays5Dw:indays5Dw+30))); %the threshold has been computed on a 31-days moving window 
                threshold.tx90p(indays5Dw,:)=prctile(Tx_th(ind,:),90);
            end
            clear indices indays5Dw ind
            
            indicator(i).Index=NaN(Nagg,Nest);
            fechas=datesym(period,'yyyymmdd');
            diasRef=datesym(period_th(1):period_th(end),'yyyymmdd');
            diasRef=unique(diasRef(:,5:end),'rows');
            for k=1:Nest
                aux=NaN(ndata,1);
                for j=1:size(diasRef,1)
                    ind=strmatch(diasRef(j,:),fechas(:,5:end));
                    ind1=find(~isnan(Tx(ind,k)) & Tx(ind,k)>threshold.tx90p(j,k));
                    ind2=find(~isnan(Tx(ind,k)) & Tx(ind,k)<=threshold.tx90p(j,k));
                    aux(ind(ind1))=1;
                    aux(ind(ind2))=0;
                end
                for j=1:Nagg
                    ind=find(datesAgg==j)';
                    indicator(i).Index(j,k)=0;ndays=0;
                    for l=ind
                        if aux(l)==1
                            ndays=ndays+1;
                        else
                            if ndays>2    % At least 3 consecutive days must exceed the threshold

                                indicator(i).Index(j,k)=indicator(i).Index(j,k)+1;
                            end
                            ndays=0;
                        end
                    end
                    if all(isnan(Tx(ind))) 
                        indicator(i).Index(j,k)=NaN;
                    end
                end
            end
        case {'hwtxdx'}
            [ndata,Nest]=size(Tx);

            threshold.tx90p=NaN(ndays5Dw,Nest);
            indices=[ndays5Dw-14:ndays5Dw 1:ndays5Dw 1:15];
            for indays5Dw=1:ndays5Dw
                ind=find(ismember(I2,indices(indays5Dw:indays5Dw+30))); % The threshold has been computed on a 31-days moving window
                threshold.tx90p(indays5Dw,:)=prctile(Tx_th(ind,:),90);
            end
            clear indices indays5Dw ind

            indicator(i).Index=NaN(Nagg,Nest);
            fechas=datesym(period,'yyyymmdd');
            diasRef=datesym(period_th(1):period_th(end),'yyyymmdd');
            diasRef=unique(diasRef(:,5:end),'rows');
            for k=1:Nest
                aux=NaN(ndata,1);
                auxtx=NaN(ndata,1);
                for j=1:size(diasRef,1)
                    ind=strmatch(diasRef(j,:),fechas(:,5:end));
                    ind1=find(~isnan(Tx(ind,k)) & Tx(ind,k)>threshold.tx90p(j,k));
                    ind2=find(~isnan(Tx(ind,k)) & Tx(ind,k)<=threshold.tx90p(j,k));
                    aux(ind(ind1))=1;
                    aux(ind(ind2))=0;
                    auxtx(ind(ind1))=Tx(ind(ind1));
                    auxtx(ind(ind2))=0;
                end
                for j=1:Nagg
                    ind=find(datesAgg==j)';
                    indicator(i).Index(j,k)=NaN;ndays=0;
                    count=1;
                    for l=ind
                        if aux(l)==1
                            ndays=ndays+1;
                            auxtx_tomediate(ndays)=auxtx(l);
                        else
                            if ndays>2 %At least 3 consecutive days must exceed the threshold
                                hwtx_mean(count)=mean(auxtx_tomediate);
                                count=count+1;
                            end
                            ndays=0;
                            clear auxtx_tomediate
                        end
                    end
                    if exist('hwtx_mean','var')
                        indicator(i).Index(j,k)=max(hwtx_mean);
                        clear hwtx_mean
                    end
                    if all(isnan(Tx(ind))) 
                        indicator(i).Index(j,k)=NaN;
                    end
                end
            end
        case 'tr'
            [ndata,Nest]=size(Tn);
            indicator(i).Index=NaN(Nagg,Nest);
            fullaux=NaN(ndata,Nest);
            for k=1:Nest
                aux=NaN(ndata,1);
                ind=find(~isnan(Tn(:,k)) & Tn(:,k)>20);aux(ind)=1;
                ind=find(~isnan(Tn(:,k)) & Tn(:,k)<=20);aux(ind)=0;
                fullaux(:,k)=aux;
            end
            indicator(i).Index=aggregateData(fullaux,period,aggregation,'aggfun','nansum','missing',missing);
            
        case 'cdd'
            [ndata,Nest]=size(Pr);
            indicator(i).Index=NaN(Nagg,Nest);
            fechas=datesym(period,'yyyymmdd');
            for k=1:Nest
                aux=NaN(ndata,1);
                ind=find(~isnan(Pr(:,k)) & Pr(:,k)<1);
                aux(ind)=1;
                clear ind
                for j=1:Nagg
                    ind=find(datesAgg==j)';
                    indicator(i).Index(j,k)=0;
                    maximo=0;
                    for l=ind
                        if aux(l)==1
                            indicator(i).Index(j,k)=indicator(i).Index(j,k)+1;
                        else
                            if indicator(i).Index(j,k)>maximo
                                maximo=indicator(i).Index(j,k);
                                indicator(i).Index(j,k)=0;
                            else
                                indicator(i).Index(j,k)=0;
                            end
                        end
                    end

                    if all(isnan(Pr(ind,k)))
                        indicator(i).Index(j,k)=NaN;
                    else
                        indicator(i).Index(j,k)=max(maximo,indicator(i).Index(j,k));
                    end
                    clear ind
                end
            end
        case {'prcptot' 'bio12'}
            Pr_wetday=Pr;Pr_wetday(Pr_wetday<1)=NaN;
            indicator(i).Index=aggregateData(Pr_wetday,period,aggregation,'aggfun','nansum','missing',missing);
        case 'r95ptot'
            Pr_th(Pr_th<1)=NaN;% Only days with precipitation higher than 1 mm are considered (Pr>=1mm)

            switch aggregation
                case 'Y'
                    threshold.r95p=prctile(Pr_th,95);
                    [ndata,Nest]=size(Pr);
                    indicator(i).Index=NaN(Nagg,Nest);
                    Pr_wetday=Pr;Pr_wetday(Pr_wetday<1)=NaN;
                    RR=aggregateData(Pr_wetday,period,aggregation,'aggfun','nansum','missing',missing);
                    for k=1:Nest
                        aux=Pr(:,k);
                        ind=find(~isnan(Pr(:,k)) & Pr(:,k)<=threshold.r95p(k));
                        aux(ind)=0;
                        aux=aggregateData(aux,period,aggregation,'aggfun','nansum','missing',missing);
                        ind=find(~isnan(RR(:,k)) & RR(:,k)~=0);
                        indicator(i).Index(ind,k)=100*aux(ind)./RR(ind,k);
                    end
                case 'S'
                    [ndata,Nest]=size(Pr);
                    %SEASONAL INDICATORS%
                    % Computing the threshold for each season
                    day_season=aggregateDataDay(Pr_th,period_th,aggregation,'aggfun','nanmean','missing',missing);
                    aux=NaN(size(day_season,2), size(Pr_th,2));
                    for seas=1:size(day_season,2)
                        aux1=day_season(:,seas);
                        aux1(isnan(aux1))=[];
                        aux(seas,:)=prctile(Pr_th(aux1,:),95);
                        clear aux1
                    end
                    threshold.r95p=aux; % For each point four percentile thresholds are defined (one per season) 
                    indicator(i).Index=NaN(Nagg,Nest);
                    Pr_wetday=Pr;Pr_wetday(Pr_wetday<1)=NaN;
                    RR=aggregateData(Pr_wetday,period,aggregation,'aggfun','nansum','missing',missing);
                    day_season_pr=aggregateDataDay(Pr,period,aggregation,'aggfun','nanmean','missing',missing);
                    for k=1:Nest
                        aux2=NaN(length(Pr),1);
                        for seas=1:4
                            aux1=day_season_pr(:,seas);
                            aux1(isnan(aux1))=[];

                            ind=find(~isnan(Pr(aux1,k)) & Pr(aux1,k)<=threshold.r95p(seas,k));
                            ind2=find(~isnan(Pr(aux1,k)) & Pr(aux1,k)>threshold.r95p(seas,k));
                            aux2(aux1(ind),1)=0;
                            aux2(aux1(ind2),1)=Pr(aux1(ind2),k);
                            clear ind ind2 aux1
                        end
                        aux3=aggregateData(aux2,period,aggregation,'aggfun','nansum','missing',missing);
                        ind=find(~isnan(RR(:,k)) & RR(:,k)~=0);
                        indicator(i).Index(ind,k)=100*aux3(ind)./RR(ind,k); 
                        clear aux2 aux3
                    end
                case 'M'
                    [ndata,Nest]=size(Pr);
                    %The percentile is computed for each season distribution 
                    day_month=aggregateDataDay_M(Pr_th,period_th,aggregation,'aggfun','nanmean','missing',missing);
                    aux=NaN(size(day_month,2), size(Pr_th,2));
                    for imonth=1:size(day_month,2)
                        aux1=day_month(:,imonth);
                        aux1(isnan(aux1))=[];
                        aux(imonth,:)=prctile(Pr_th(aux1,:),95);
                        clear aux1
                    end
                    %MONTHLY INDICATORS%
                    threshold.r95p=aux; %The percentiles are computed for each month: 12 thresholds defined.
                    indicator(i).Index=NaN(Nagg,Nest);
                    Pr_wetday=Pr;Pr_wetday(Pr_wetday<1)=NaN;
                    RR=aggregateData(Pr_wetday,period,aggregation,'aggfun','nansum','missing',missing);
                    day_month_pr=aggregateDataDay_M(Pr,period,aggregation,'aggfun','nanmean','missing',missing);
                    for k=1:Nest
                        aux2=NaN(length(Pr),1);
                        for imonth=1:12
                            aux1=day_month_pr(:,imonth);
                            aux1(isnan(aux1))=[];

                            ind=find(~isnan(Pr(aux1,k)) & Pr(aux1,k)<=threshold.r95p(imonth,k));
                            ind2=find(~isnan(Pr(aux1,k)) & Pr(aux1,k)>threshold.r95p(imonth,k));
                            aux2(aux1(ind),1)=0;
                            aux2(aux1(ind2),1)=Pr(aux1(ind2),k);
                            clear ind ind2 aux1
                        end
                        [aux3,dateM]=aggregateData(aux2,period,aggregation,'aggfun','nansum','missing',missing);
                        ind=find(~isnan(RR(:,k)) & RR(:,k)~=0);
                        indicator(i).Index(ind,k)=100*aux3(ind)./RR(ind,k); %100 PER precipitazione solo dei giorni che superano il percentile/totale
                        clear aux2 aux3
                    end
            end
        case 'rx1day'
            indicator(i).Index=aggregateData(Pr,period,aggregation,'aggfun','nanmax','missing',missing);
        case 'rx5day'
            [ndata,Nest]=size(Pr);
            indicator(i).Index=NaN(Nagg,Nest);
            fechas=datesym(period,'yyyymmdd');
            for j=1:Nagg
                ind=find(datesAgg==j)';
                aux=NaN(length(ind)-4,Nest);
                for k=1:Nest
                    for l=1:length(ind)-4
                        if any(all(isnan(Pr(ind(l):ind(l)+4,k))))
                            aux(l,k)=NaN;
                        else
                            aux(l,k)=nansum(Pr(ind(l):ind(l)+4,k));
                        end
                    end
                end
                indicator(i).Index(j,:)=nanmax(aux);
                clear ind
            end
            %Hum
        case {'rh' 'hurs'}
            indicator(i).Index=aggregateData(Hum,period,aggregation,'aggfun','nanmean','missing',missing);
        case 'humidexv'
            exphum=6.11*10.^(7.5*Tx./(237.7+Tx)).*Hum./100;
            humidex=Tx+0.5555*(exphum-10);
            indicator(i).Index=aggregateData(humidex,period,aggregation,'aggfun','nanmean','missing',missing);
        case 'huxwf'
            % The perceived temperature (humidex) is here computed
            exphum=6.11*10.^(7.5*Tx./(237.7+Tx)).*Hum./100;
            humidex=Tx+0.5555*(exphum-10);
            [ndata,Nest]=size(humidex);
            indicator(i).Index=NaN(Nagg,Nest);
            for k=1:Nest
                aux=NaN(ndata,1);
                ind=find(~isnan(humidex(:,k)) & humidex(:,k)>35);aux(ind)=1;
                ind=find(~isnan(humidex(:,k)) & humidex(:,k)<=35);aux(ind)=0;
                for j=1:Nagg
                    ind=find(datesAgg==j)';
                    indicator(i).Index(j,k)=0;ndays=0;
                    for l=ind
                        if aux(l)==1
                            ndays=ndays+1;
                        else
                            if ndays>3
                                indicator(i).Index(j,k)=indicator(i).Index(j,k)+ndays;
                            end
                            ndays=0;
                        end
                    end
                    if all(isnan(Tx(ind))) || all(isnan(Hum(ind))) 
                        indicator(i).Index(j,k)=NaN;
                    end
                end
            end
        case 'fg'
            indicator(i).Index=aggregateData(Ws,period,aggregation,'aggfun','nanmean','missing',missing);
            % Number of days with wind speed >  5 (moderate), 10 (strong), 15 (very strong) or 25 (storm) m/s 
        case 'fgmax'
            indicator(i).Index=aggregateData(Wsmax,period,aggregation,'aggfun','nanmean','missing',missing);

        otherwise, warning(sprintf('Unknown indicator'));

    end

end

end