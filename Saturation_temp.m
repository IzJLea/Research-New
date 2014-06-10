%% system input

Pin=11.380; % MPa 

hin=1500; %enthalpy entering channel. kJ/kg

%% Retrieving H2O physical data (for interpolation)
%Saturation Pressure

[~, ~, raw] = xlsread('C:\Users\Izaak\Documents\Research\H2O_TempSat.xls','Sheet1','B32:B52');

WaterSaturationPropertiesTemperatureTable8 = reshape([raw{:}],size(raw));

clearvars raw;

Ph2o=WaterSaturationPropertiesTemperatureTable8;

% Temperature Data
[~, ~, raw] = xlsread('C:\Users\Izaak\Documents\Research\H2O_TempSat.xls','Sheet1','A32:A52');

WaterSaturationPropertiesTemperatureTable2 = reshape([raw{:}],size(raw));

clearvars raw;

Th2o=WaterSaturationPropertiesTemperatureTable2;


Tsys=Lagint(Ph2o,Th2o,Pin);

display(Tsys);
