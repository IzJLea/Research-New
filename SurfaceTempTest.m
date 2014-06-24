%% Temperature calculation

% channel as one uniform property channel

Qchannel=5.5; %MW

Tin=267.2; %C

M=25.8; % kg/s

Lchannel=6; %m


Tsingle=Tsurfaces(Tin,Qchannel,M,Lchannel);

plot(Tsingle);

