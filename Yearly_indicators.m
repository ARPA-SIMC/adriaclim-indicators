% TO CLEAN THE CACHE MEMORY 
clear all; close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Authors: Giusy Fedele and Veronica Villani
%Affiliation: Regional Models and Geo-Hydrological Impacts Division (REMHI), Centro Euro-Mediterraneo sui Cambiamenti Climatici, 81100 Caserta, Italy

%Disclaimer for codes provided by REMHI division:
%The users must be expert on their specific application. Therefore, REMHI division is not responsible for inappropriate use of codes and analysis provided. To use the codes provided, it is requested to mention REMHI Division as developer and provider.
%The MeteoLab open-source Matlab toolbox (https://meteo.unican.es/trac/MLToolbox/wiki) has been used to estimate the indicators above-mentioned, adapting the available codes to the AdriaClim Project needs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %%%   OCTAVE  %%%
% TO LOAD OCTAVE PACKAGES NEEDED
%Put "%" symbol if working in Octave in the following three lines:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load statistics
pkg load io
pkg load netcdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off', 'all'); %Disable warnings
% DEFINING PATHS
dir_input='/Users/giusyfedele/Dropbox (CMCC)/ADRIACLIM_xGiusy/periodi/' ; %To define the path where the input data are collected 
dir_output='/Users/giusyfedele/Dropbox (CMCC)/ADRIACLIM_xGiusy/periodi/indicatori_annuali/prova/' ; %To define  the path where to locate the outputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE INDICATORS:

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

regions={'Marche';'Puglia'}; %To define Pilot Area

%To define vector containing the indicators to be computed (vettindex)
%To define vector containing the respective unit of measurement (udm)

vettindex= {'wsdi' 'hwn' 'hwtxdx' 'txx' 'txn' 'tg' 'tr' 'tx90p' 'fg' ...
    'prcptot' 'rx1day' 'rx5day' 'r95ptot' 'cdd'};
udm= {'days' 'number of heat waves' '°C' '°C' '°C' '°C' 'days' 'days' 'm/s' ...
    'mm' 'mm' 'mm' '%' 'days'};

%For each Pilot Area the indicators are computed for two 30-years epochs here defined respectively: 1981-2010 and 1991-2020.    
    
for iregion=1:length(regions)
    periods={'1981_2010' '1991_2020'};
    startDates_CHAR={'01-Jan-1981' '01-Jan-1991'};
    endDates_CHAR={'31-Dec-2010' '31-Dec-2020'};
    
    % A for cycle loop is here defined in order to calculate every indicator for each epoch already defined. 

    for iperiod=1:length(periods)

        nome_file=[dir_input,regions{iregion} '_ERA5_Land_1d_daily_' char(periods(iperiod)) '.mat']; 
        load(nome_file)

        % UNIT OF MEASURE CONVERSION       
        tas=t2mean-273.15;      %Mean Temperature at 2 meters - Conversion from degrees Kelvin (°K) to degrees Celsius (°C). N.B. To comment this line if temperatures are already in degrees Celsius.
        tasmin=t2min-273.15;    %Min Temperature at 2 meters - Conversion from degrees Kelvin (°K) to degrees Celsius (°C). N.B. To comment this line if temperatures are already in degrees Celsius.
        tasmax=t2max-273.15;    %Max Temperature at 2 meters - Conversion from degrees Kelvin (°K) to degrees Celsius (°C). N.B. To comment this line if temperatures are already in degrees Celsius.
        pr=prec_cum*1000;       %Cumulated Precipitation at surface - Conversion from meters (m) to millimiters (mm). N.B. To comment this line if the precipitation rate is already in millimiters.

        % TO DEFINE THE VARIABLES WHERE A THRESHOLD WILL BE IMPOSED
        tas_th=t2mean-273.15;   %Mean Temperature at 2 meters - Conversion from degrees Kelvin (°K) to degrees Celsius (°C). N.B. To comment this line if temperatures are already in degrees Celsius.
        tasmin_th=t2min-273.15; %Min Temperature at 2 meters - Conversion from degrees Kelvin (°K) to degrees Celsius (°C). N.B. To comment this line if temperatures are already in degrees Celsius.
        tasmax_th=t2max-273.15; %Max Temperature at 2 meters - Conversion from degrees Kelvin (°K) to degrees Celsius (°C). N.B. To comment this line if temperatures are already in degrees Celsius.
        pr_th=prec_cum*1000;    %Cumulated Precipitation at surface - Conversion from meters (m) to millimiters (mm). N.B. To comment this line if the precipitation rate is already in millimiters.
        date_th=date;           %To import datetime
        
        % TO CLEAR VARIABLES NOT NEEDED ANYMORE TO FREE UP MEMORY
        clear date prec_cum t2mean t2max t2min u10 v10 

        %FROM DATETIME TO DATENUM
        startDates=startDates_CHAR{iperiod};
        endDates=endDates_CHAR{iperiod};
        period=[datenum(startDates):1:datenum(endDates)];

        % TO APPLY THE FUNCTION: extremesIndicators_adriaclim
        indicators=extremesIndicators_adriaclim(1,'Y',tas,tasmin,tasmax,pr,[],[],[],wind10,[],[],tas_th,tasmin_th,tasmax_th,pr_th,period,date_th,vettindex(6));

        clear t2mean t2max t2min pr u10 v10 wind10 tas_th tasmin_th tasmax_th pr_th date_th startDates endDates

        % TO SAVE INDICATORS AS .mat
        filename=strcat('Indicators_',num2str(regions{iregion}),'_',num2str(periods{iperiod}),'.mat');
        % TO SAVE INDICATORS AS .nc
        filename_nc=strcat('Indicators_',num2str(regions{iregion}),'_',num2str(periods{iperiod}),'.nc');

        % TO CALCULATE TREND
        for iindex=1:length(indicators)
            trend_tot=nanmean(indicators(iindex).Index,2);
            [pvalue, trendp, ~]=testTrend(trend_tot,'test','mannkendall');
            if pvalue<0.05
                fid=fopen([dir_output,'TestTrend_' num2str(regions{iregion}),'_',num2str(periods{iperiod}) '.txt'],'a');
                fprintf(fid,'\t The trend %2.2f of the %s indicator is statistically significant\n',trendp,indicators(iindex).Name);
                fclose(fid);
                clear fid
            else
                fid=fopen([dir_output,'TestTrend_' num2str(regions{iregion}),'_',num2str(periods{iperiod}) '.txt'],'a');
                fprintf(fid,'\t The trend %2.2f of the %s indicator is not statistically significant\n',trendp,indicators(iindex).Name);
                fclose(fid);
                clear fid
            end
        end      

        % To define dimensions %
        nt=length(period); npoints=size(tas,2);
       
       % TO SAVE INDICATORS AS .netcdf
        
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
        % Note: The symbol * indicates "calculated for a 5 day window centred on each calendar day in the 1971-2000 period or the period on which the threshold is calculated".   
        % Note: The symbol ** indicates "calculated for a 31 day window centred on each calendar day in the 1971-2000 period or the period on which the threshold is calculated".   

        % To call each indicator from the structured file "indicators"
        i_wsdi=struct2cell(indicators(1));      i_wsdi=cell2mat(i_wsdi(2));
        i_hwn=struct2cell(indicators(2));       i_hwn=cell2mat(i_hwn(2));
        i_hwtxdx=struct2cell(indicators(3));    i_hwtxdx=cell2mat(i_hwtxdx(2));
        i_txx=struct2cell(indicators(4));       i_txx=cell2mat(i_txx(2));
        i_txn=struct2cell(indicators(5));       i_txn=cell2mat(i_txn(2));
        i_tg=struct2cell(indicators(6));        i_tg=cell2mat(i_tg(2));
        i_tr=struct2cell(indicators(7));        i_tr=cell2mat(i_tr(2));
        i_tx90p=struct2cell(indicators(8));     i_tx90p=cell2mat(i_tx90p(2));
        i_fg=struct2cell(indicators(9));        i_fg=cell2mat(i_fg(2));
        i_prcptot=struct2cell(indicators(10));  i_prcptot=cell2mat(i_prcptot(2));
        i_rx1day=struct2cell(indicators(11));   i_rx1day=cell2mat(i_rx1day(2));
        i_rx5day=struct2cell(indicators(12));   i_rx5day=cell2mat(i_rx5day(2));
        i_r95ptot=struct2cell(indicators(13));  i_r95ptot=cell2mat(i_r95ptot(2));
        i_cdd=struct2cell(indicators(14));      i_cdd=cell2mat(i_cdd(2));
        
        
        fillvalue=1.0e+36; % To substitute NaN with fillvalue where present
        nccreate(filename_nc,'wsdi','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'wsdi',i_wsdi);
        nccreate(filename_nc,'hwn','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'hwn',i_hwn);
        nccreate(filename_nc,'hwtxdx','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'hwtxdx',i_hwtxdx);
        nccreate(filename_nc,'txx','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'txx',i_txx);
        nccreate(filename_nc,'txn','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'txn',i_txn);
        nccreate(filename_nc,'tg','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'tg',i_tg);
        nccreate(filename_nc,'tr','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'tr',i_tr);
        nccreate(filename_nc,'tx90p','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'tx90p',i_tx90p);
        nccreate(filename_nc,'fg','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'fg',i_fg);
        nccreate(filename_nc,'prcptot','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'prcptot',i_prcptot);
        nccreate(filename_nc,'rx1day','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'rx1day',i_rx1day);
        nccreate(filename_nc,'rx5day','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'rx5day',i_rx5day);
        nccreate(filename_nc,'r95ptot','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'r95ptot',i_r95ptot);
        nccreate(filename_nc,'cdd','Dimensions', {'x',nt,'y',npoints},'FillValue',fillvalue); ncwrite(filename_nc,'cdd',i_cdd);
        nccreate(filename_nc,'period','cdd','Dimensions', {'x',nt},'FillValue',fillvalue); ncwrite(filename_nc,'period',period);
        ncdisp(filename_nc)
        
        save(filename,'indicators'); % save indicators as .mat
        
end
        % TO CLEAR VARIABLES NOT NEEDED ANYMORE TO FREE UP MEMORY
        clear period indicators

    end

    % TO CLEAR VARIABLES NOT NEEDED ANYMORE TO FREE UP MEMORY
    clear periods startDates_CHAR endDates_CHAR

