function [pValue,trend,T]=testTrend(data,varargin)
% [pValue,trend]=testTrend(data,varargin);
% 
% Funcion que calcula el p-valor y la tendencia de los datos.
% 
% En la entrada : 
% 	data        : son las series de datos de las cuales queremos estudiar la
%                 tendencia.
% 	varargin	: paramétros opcionales
% 		'test'	-	toma los valores 'Spearman' o 'MannKendall' (MannKendall por defecto)
% 		'period'-	toma los valores 'day','month' y 'year'. Indica
% 		            como agrupamos los datos para la aplicacion de los tests. Por
%            		defecto se le asigna 'day'.
%       'missing' - parametro de missing data para la funcion 'movingAverage'
% 	
% En la salida :
% 	pValue : vector de longitud Nest (numero de estaciones) con los
%            p-valores del test en cada una de las estaciones.
% 	trend  : vector de longitud Nest que contiene la tendencia de la serie
% 	T        : vector de longitud Nest que contiene el valor del estadistico
% 
% Ejemplo de llamada:
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
% Nest=size(data,2);
% trend=zeros(Nest,1)+NaN;
% pValue=zeros(Nest,1)+NaN;
% for k=1:Nest
    % switch lower(test) 
        % case('spearman')
            % [pValue(k),trend(k)]=SP(data(:,k),'period',period,'missing',missing);
        % case('mannkendall')
            % [pValue(k),trend(k)]=MK(data(:,k),'period',period,'missing',missing);
    % end
% end

Nest=size(data,2);
trend=zeros(Nest,1)+NaN;
pValue=zeros(Nest,1)+NaN;
switch lower(test) 
    case('spearman')
        [pValue,trend,T]=SP(data,'period',period,'missing',missing);
    case('mannkendall')
        [pValue,trend,T]=MK(data,'period',period,'missing',missing,'autocorrelation',autocorrelation);
end
