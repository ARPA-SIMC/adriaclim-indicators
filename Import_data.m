clear all; close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Authors: Giusy Fedele and Veronica Villani
%Affiliation: Regional Models and Geo-Hydrological Impacts Division (REMHI), Centro Euro-Mediterraneo sui Cambiamenti Climatici, 81100 Caserta, Italy

%Disclaimer for codes provided by REMHI division:
%The users must be expert on their specific application. Therefore, REMHI division is not responsible for inappropriate use of codes and analysis provided. To use the codes provided, it is requested to mention REMHI Division as developer and provider.
%The MeteoLab open-source Matlab toolbox (https://meteo.unican.es/trac/MLToolbox/wiki) has been used to estimate the indicators above-mentioned, adapting the available codes to the AdriaClim Project needs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TO IMPORT NETCDF FILES AND CONVERT THEM IN .MAT 
% LOADING THE PACKAGES NEEDED TO RUN THE SCRIPTS %
pkg load statistics
pkg load io
pkg load netcdf
% Please note that the input and output files are saved as:
%[dir_input, regions{iregion} '/' variables{ivar} '_1d_daily_ERA5_1950_2021.nc']
%[dir_input, regions{iregion} '/' variables{ivar} '_1d_daily_ERA5_1950_2021.mat']
%where dir_input is the input directory path, 
%regions{iregion} indicate the name of the Pilot Area; 
%variables{ivar} indicate the name of the variable (Example:t2max,t2min, ...) 

dir_input='/Users/giusyfedele/Dropbox (CMCC)/ADRIACLIM_xGiusy/ERA5_Land/'; %To define the path where the input data are collected 

regions={'Marche' 'Puglia'}; % To define the name of the Pilot Area


variables={'prec_cum' 't2mean' 't2max' 't2min' 'u10' 'v10' 'wind10'};%To define vector containing the indicators to be computed
variables_infilenc={'tp' 't2m' 't2m' 't2m' 'u10' 'v10' 'u10'};%To define vector containing the name of each variable in the netcdf files
udm={'m' 'k' 'k' 'k' 'm/s' 'm/s' 'm/s' };%To define vector containing the respective unit of measurement 

for ivar=1:length(variables)

    for iregion=1:length(regions)

        nome_file=[dir_input, regions{iregion} '/' variables{ivar} '_1d_daily_ERA5_1950_2021.nc']; % To define the name of the input file 
        variab = ncread(nome_file,variables_infilenc{ivar}); % It reads the input file
        
        % THE MATRIX with the number of points (rows) x number of times (columns) is defined
        variab=rot90(variab);
        variab = reshape(variab,[],size(squeeze(variab),3));
        variab = variab'; % The matrix is trasposed so that we get: number of times (rows) x number of points (columns) 


        % DELETE THE MODEL GRID POINTS WHICH HAVE ALL NAN VALUES
        Nest=size(variab,2);
        Ndays=size(variab,1);
        j=1;
        for i=1:Nest
            ctrl=sum(isnan(variab(:,i)));
            if(ctrl==Ndays)
                indici(j)=i;
                j=j+1;
            end
        end        
        variab(:,indici)=[];

        if(ivar==1)
            %TO READ COORDINATES
            lon= ncread(nome_file,'longitude');
            lat= ncread(nome_file,'latitude');
            [LON, LAT]=meshgrid(lon, flipud(lat));
            coord=[double(LON(:)) double(LAT(:))];
            
            coord(indici,:)=[]; % The coordinates without values are ignored

            save(['ERA5_Land_coord_' regions{iregion} '.mat'],'coord'); %save coordinate as.mat
            
            %save coordinates as netcdf
            fillvalue=1.0e+36; % To substitute NaN with fillvalue where present
            filename_nc=strcat('coordinates_',num2str(regions{iregion}),'.nc');
            nlat=size(LON,1); nlon=size(LON,2);
            nccreate(filename_nc,'longitude','Dimensions', {'x',nlat,'y',nlon},'FillValue',fillvalue); ncwrite(filename_nc,'longitude',LON);
            nccreate(filename_nc,'latitude','Dimensions', {'x',nlat,'y',nlon},'FillValue',fillvalue); ncwrite(filename_nc,'latitude',LAT);           
        end
        
        clear indici %restore indici

        eval([variables{ivar} '=variab;']);
        clear variab

        %SAVE DATA IN .MAT
        nome_file_mat=[dir_input, regions{iregion} '/' variables{ivar} '_1d_daily_ERA5_1950_2021.mat']; % To define the name of the output file
        save(nome_file_mat,variables{ivar});

    end

end

