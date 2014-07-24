%
% Water properties for given pressure and temperature
function [rho,beTa,my,nu,k,alpha,cp,h,hf_p,hfg_p,rhoL,rhoV,hg_p] = fluid_props(P,T)
%Fluid Properties at temperature T[C] and pressure P[atm]

if nargin==1
disp('Only pressure given, therefore all values taken as values at saturation at given pressure')
T=XSteam('Tsat_p',P);
end

rho=XSteam('rho_pt',P,T); %[kg/m^3] Water density
if isnan(rho)==1
disp('The given temperature gives saturated conditions at given pressure, saturated liquid density given by default');
rho=XSteam('rhoL_p',P);
end

beTa=XSteam('v_pt',P,T)^-1*(XSteam('v_pt',P,T+1e-3)-XSteam('v_pt',P,T))/1e-3; %[1/K] Thermal/Volumetric expansion coefficient
%for now beta is calculated from my approximation to the derivative that defines it. It's within less than 1% error between 10<T<88 Celcius
if ((beTa>0.00075048) || (isnan(beTa)==1)) %0.00075048 according to [miniREFPROP]
beTa=0.00075048; %beTa at saturation for atm pressure less than 1% difference for pressures between 1 and 45 atmospheres
end

k=XSteam('tc_pT',P,T); %[W/m-K] Thermal conductivity
if isnan(k)==1
k=XSteam('tcL_p',P);
end

my=XSteam('my_pT',P,T); %[Pa-s]=[kg/(s*m)] Viscosity
if isnan(my)==1
k=XSteam('tcL_p',P);
end

nu=my/rho; %[m^2/s] Kinematic viscosity

cp=XSteam('Cp_pT',P,T)*1e3; %[J/(kg-K)] Water specific heat capacity
if isnan(cp)==1
cp=XSteam('CpL_p',P)*1e3;
end

alpha=k/(rho*cp); %[m^2/s] Thermal diffusivity

h=XSteam('h_pT',P,T)*1e3; %[J/kg] Specific enthalpy
if isnan(my)==1
disp('The given temperature gives saturated conditions at given pressure, please choose either saturated liquid enthalpy or saturated vapor enthalpy');
h=0;   
end

hf_p=XSteam('hL_p',P)*1e3; %[J/kg] Liquid saturation enthalpy vaporization @ given pressure
hg_p=XSteam('hV_p',P)*1e3;
hfg_p=hg_p-hf_p; %[J/kg] Latent heat of vaporization @ given pressure
rhoL=XSteam('rhoL_p',P); %[kg/m^3] Saturated liquid density
rhoV=XSteam('rhoV_p',P); %[kg/m^3] Saturated vapor density

end
