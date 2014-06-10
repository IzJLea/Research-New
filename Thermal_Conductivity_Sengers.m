%% system information input

Minput=25.8;    %kg/s coolant flowrate

Rhoinl=1000;   %kg/m^3 liquid water density at saturation conditions

Rhoinv=52.6;    %kg/m^3 water vapour density at saturation conditions

Hlin=1450;  %kJ/kg enthalpy of liquid water at saturation conditions

Hvin=2706;  %kJ/kg enthalpy of gaseous water at saturation conditions

Hin=1550;   %kJ/kg enthalpy of fluid entering channel

Dh=0.0074;  %m equivalent hydraulic diameter of channel

Pin=11.380e6; %Pa pressure at channel inlet 

Tin=273.15; %K saturation temperature at inlet pressure

%% Two Phase Pressure drop
% water viscosity calculation

Ts=647.27; % K  reference constant temperature

rhos=317.763; % kg/m^3 refernce constant critical density

Ps=22.115e6; % Pa reference constant pressure

mus=1e-6;    %Pa*s reference constant viscosity

Tr=Tin/Ts; % Dimensionless reference temperature

rhor=Rhoinl/rhos; % Dimensionless reference density

Pr=Pin/Ps; % Dimensionless reference pressure

lams=1; %W/m.K  reference conductivity

% lam=lam0xlam1+lam2
%lam0 calculation

k=linspace(0,3,4);

ak=[2.02223, 14.11166, 5.25597, -2.01870];

Trk=Tr.^(-k);

Dem0=ak.*Trk;

lam0=mus*sqrt(Tr)/(sum(Dem0));

display(lam0);

%lam1 calculation

i=linspace(0,4,5);

j=linspace(0,5,6);

j=j';

bij=[1.3293046, 1.7018363, 5.2246158, 8.7127675, -1.8525999;
    -0.40452437, -2.2156845, -10.124111, -9.5000611, 0.93404690;
    0.24409490, 1.6511057, 4.9874687, 4.3786606, 0.0;
    0.018660751, -0.76736002, -0.27297694, -0.91783782, 0.0;
    -0.12961068, 0.37283344, -0.43083393, 0.0, 0.0;
    0.044809953, -0.11203160, 0.13333849, 0.0, 0.0];


Tinv=((1/Tr)-1).^i;

rhorm=((rhor-1).^j);

lam10=rhorm*Tinv;

lam11=lam10.*bij;

lam12=sum(lam11);

lam13=sum(lam12);

lam1=exp(rhor*lam13);

display(lam1);

%lam2 calculation

% water viscosity calculation

Ts=647.27; % K  reference constant temperature

rhos=317.763; % kg/m^3 refernce constant critical density

Ps=22.115e6; % Pa reference constant pressure

mus=1e-6;    %Pa*s reference constant viscosity

Tr=Tin/Ts; % Dimensionless reference temperature

rhor=Rhoinl/rhos; % Dimensionless reference density

Pr=Pin/Ps; % Dimensionless reference pressure

% mu=mu0xmu1 
%mu0 calculation

k=linspace(0,3,4);

ak=[0.0181583, 0.0177624, 0.0105287, -0.0036744];

Trk=Tr.^(-k);

Dem0=ak.*Trk;

mu0=mus*sqrt(Tr)/(sum(Dem0));

display(mu0);

%mu1 calculation

i=linspace(0,5,6);

j=linspace(0,4,5);

j=j';

bij=[0.501938,0.162888,-0.130356,0.907919,-0.551119,0.146543; 
    0.235622,0.789393,0.673665,1.207552,0.0670665,-0.0843370;
    -0.274637,-0.743539,-0.959456,-0.687343,-0.497089,0.195286;
    0.145831,0.263129,0.347247,0.213486,0.100754,-0.032932;
    -0.0270448,-0.0253093,-0.0267758,-0.0822904,0.0602253,-0.0202595];


Tinv=((1/Tr)-1).^i;

rhorm=((rhor-1).^j);

mu10=rhorm*Tinv;

mu11=mu10.*bij;

mu12=sum(mu11);

mu13=sum(mu12);

mu1=exp(rhor*mu13);

display(mu1);


mu=mu0*mu1;

display(mu);



lam=(lam0*lam1)+lam2;

display(lam);
% Calculation of Ybp from Levy model (1967)

