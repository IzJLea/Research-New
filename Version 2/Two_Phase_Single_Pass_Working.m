%%Inputs and property calculation 

Qchannel=5:0.1:12; %MW

Qchannelimp=Qchannel*1000*3412.14; %Btu/h

Tin=267.2; %C

Achannel=0.0034;%m^2

Achannelimp=Achannel*3.2808^2; %ft^2

doutclad=0.0138; %m

roc=doutclad/2;

dfuel=0.0122; %m

tclad=0.00038; %m

Lchannel=5.94; %m

%% Initial single pass k calculation- 
% pressure differences in each pass section
%



DP = [165;
262.5000;
520.5000;
258;
174];


%%calculate values for k for each section

% reference mass flowrate

M=25.8; % kg/s

Rho=780.6; %kg/m^3

keff=DP./1000./(M^2);

reff=keff.*Rho;

reffT=sum(reff);

Pout=10; %MPa

Pin=((Pout*1000)+sum(DP))/1000;  %MPa


%% Physical data import
% Temperature for water properties (C)


Tsat=[100
110
120
130
140
150
160
170
180
190
200
210
220
230
240
250
260
270
280
290
300
310
320
330
340
350
360
370
373.9500];

Tsat=reshape(Tsat,29,1);

% Saturation Pressure (MPa)

Psat=[0.1014
0.1434
0.1987
0.2703
0.3615
0.4762
0.6182
0.7922
1.0028
1.2552
1.5549
1.9077
2.3196
2.7971
3.3469
3.9762
4.6923
5.5030
6.4166
7.4418
8.5879
9.8651
11.2840
12.8580
14.6010
16.5290
18.6660
21.0440
22.0640];

Psat=reshape(Psat,29,1);

% Fluid specific volume

svfluid=[0.0010
0.0011
0.0011
0.0011
0.0011
0.0011
0.0011
0.0011
0.0011
0.0011
0.0012
0.0012
0.0012
0.0012
0.0012
0.0013
0.0013
0.0013
0.0013
0.0014
0.0014
0.0014
0.0015
0.0016
0.0016
0.0017
0.0019
0.0022
0.0031];

svfluid=reshape(svfluid,29,1);

% density conversion (kg/m^3)

rhof=1./svfluid;

% vapour specific volume

svvapour=[1.6718
1.2093
0.8912
0.6680
0.5085
0.3925
0.3068
0.2426
0.1938
0.1564
0.1272
0.1043
0.0861
0.0715
0.0597
0.0501
0.0422
0.0356
0.0302
0.0256
0.0217
0.0183
0.0155
0.0130
0.0108
0.0088
0.0069
0.0050
0.0031];

svvapour=reshape(svvapour,29,1);

% density conversion (kg/m^3)

rhov=1./svvapour;

% fluid enthalpy (kJ/kg)

hf=[419.1700
461.4200
503.8100
546.3800
589.1600
632.1800
675.4700
719.0800
763.0500
807.4300
852.2700
897.6300
943.5800
990.1900
1.0376e+03
1.0858e+03
1135
1.1853e+03
1.2369e+03
1290
1345
1.4022e+03
1.4622e+03
1.5259e+03
1.5945e+03
1.6709e+03
1.7617e+03
1.8907e+03
2.0843e+03];

hf=reshape(hf,29,1);

% vapour enthalpy (kJ/kg)

hv=[2.6756e+03
2.6911e+03
2.7059e+03
2.7201e+03
2.7334e+03
2.7459e+03
2.7574e+03
2.7679e+03
2.7772e+03
2.7853e+03
2792
2.7973e+03
2.8009e+03
2.8029e+03
2803
2.8009e+03
2.7966e+03
2.7897e+03
2.7799e+03
2.7667e+03
2.7496e+03
2.7279e+03
2.7006e+03
2666
2.6218e+03
2.5636e+03
2.4815e+03
2.3345e+03
2.0843e+03];

% latent heat (kJ/kg)

hfv=[2.2564e+03
2.2297e+03
2.2021e+03
2.1737e+03
2.1442e+03
2.1137e+03
2.0819e+03
2.0488e+03
2.0141e+03
1.9779e+03
1.9397e+03
1.8997e+03
1.8573e+03
1.8127e+03
1.7654e+03
1.7151e+03
1.6616e+03
1.6044e+03
1543
1.4767e+03
1.4046e+03
1.3257e+03
1.2384e+03
1.1401e+03
1.0273e+03
892.7000
719.8000
443.8000
0];

% Compressed water Pressures

Pcomp=[5;10;15]; % MPa

% Compressed water enthalpy (10 MPa, 260 C)

hcomp=[1.1349e+03;1.1343e+03;1.5924e3]; %kJ/kg

% Temp for surface tension

Tsigma=[190
195
200
205
210
215
220
225
230
235
240
245
250
255
260
265
270
275
280
285
290
295
300
305
310
315
320
325
330
335
340
345
350
355
360
365
370
374];

% Surface tension mN/m

sigma=[39.9500
38.8300
37.6900
36.5500
35.4100
34.2500
33.1000
31.9300
30.7700
29.6000
28.4200
27.2400
26.0600
24.8700
23.6700
22.4800
21.3000
20.1100
18.9400
17.7700
16.6100
15.4500
14.3000
13.1700
12.0400
10.9200
9.8100
8.7300
7.6600
6.6100
5.5900
4.6000
3.6500
2.7500
1.9000
1.1300
0.4500
0];

% Steam heat capacity Pressures (MPa)

PCp=[617.8;791.7;1002.1;1254.4;1553.8;2318;3344;4688;6412;8581;11274;14586;18651]./1000;

% liquid Heat capacity (J/kg.K)

Cpl=[4.340;4.370;4.410;4.460;4.500;4.610;4.760;4.970;5.280;5.750;6.540;8.240;14.96];

% vapour Heat capacity (J/kg.K)

Cpv=[2420; 2490; 2590;2710;2840;3110;3520;4070;4835;5980;7900;11870;25800]; 

% liquid thermal conductivity (W/m.K)

kl=[0.680;0.677;0.673;0.669;0.663;0.650;0.632;0.609;0.581;0.548;0.509;0.469;0.427];

% vapour thermal conductivity (W/m.K)

kv=[0.0331;0.0347;0.0364;000382;0.0401;0.0442;0.0487;0.0540;0.0605;0.0695;0.0836;0.110;0.178];

% liquid dynamic viscosity (kg/m.s)

mul=[0.170e-3;0.160e-3;0.150e-3;0.142e-3;0.134e-3;0.122e-3;0.111e-3;0.102e-3;0.094e-3;0.086e-3;0.078e-3;0.070e-3;0.060e-3];

% vapor dynamic viscosity (kg/m.s)

muv=[1.434e-5;1.468e-5;1.502e-5;1.537e-5;1.571e-5;1.641e-5;1.712e-5;1.788e-5;1.870e-5;1.965e-5;2.084e-5;2.255e-5;2.571e-5];

% liquid Prandtl Number

Prl=[1.09;1.03;0.983;0.947;0.910;0.865;0.836;0.832;0.854;0.902;1.00;1.23;2.06];

% vapor Prandtl Number

Prv=[1.05;1.05;1.07;1.09;1.11;1.15;1.24;1.35;1.49;1.69;1.97;2.43;3.73];

%% System Properties
% Compressed water enthalpy 

hin=interp1(Pcomp,hcomp,Pout); %kJ/kg

% saturation temperature 

Tsys=interp1(Psat,Tsat,Pout); %C

% saturation liquid enthalpy

hfsys=interp1(Psat,hf,Pout); %kJ/kg

% saturation vapour enthalpy

hvsys=interp1(Psat,hv,Pout); %kJ/kg

% saturation latent heat

hfvsys=interp1(Psat,hfv,Pout); %kJ/kg

% saturation liquid density

rhofsys=interp1(Psat,rhof,Pout);  %kg/m^3

% saturation gas density

rhovsys=interp1(Psat,rhov,Pout); %kg/m^3

% Cp Liquid  

Cplsys=interp1(PCp,Cpl,Pout); %kJ/kg.K

% Cp vapour

Cpvsys=interp1(PCp,Cpv,Pout); %kJ/kg.K

% liquid Thermal Conductivity

klsys=interp1(PCp,kl,Pout); %W/m.K

% vapour Termal Conductivity 

kvsys=interp1(PCp,kv,Pout); %W/m.K

% Liquid Dynamic Viscosity

mulsys=interp1(PCp,mul,Pout); %kg/m.s

% Vapour Dynamic Viscosity

muvsys=interp1(PCp,muv,Pout); %kg/m.s

% Liquid Prandtl Number

Prlsys=interp1(PCp,Prl,Pout); 

% Vapour Prandtl Number

Prvsys=interp1(PCp,Prv,Pout);

% hout calculation

hout=hin+(Qchannel.*1000./M); %kJ/kg

% surface tension

sigmasys=interp1(Tsigma,sigma,Tsys);

% Hydraulic Diameter

Dh=0.0074; %m

% Heat transfer area

Aht=9.1224; %m^2

%% Mass flow calculation

x=(hout-hfsys)./(hvsys-hfsys);

n=length(Qchannel);

Mchannel=zeros(1,n);

alphas=zeros(1,n);

for in=1:n
    if x(in)<=0
    
    Mchannel(in)=sqrt((Pin-Pout)*rhofsys/reffT);
    else
    check=1;
    eta=1;
    alpha=x(in);
    
    while eta>=0.001
        
        xLevy=((alpha*(1-(2*alpha)))+(alpha*sqrt((1-(2*(x(in))))+(alpha*((2*rhovsys/rhofsys*((1-alpha)^2))))+(alpha*(1-(2*alpha))))))/((2*rhovsys/rhofsys*(1-alpha))+(alpha*(1+(2*alpha))));
    
        eta=x(in)-xLevy;
        
        alpha=alpha-eta;
        
    end
    clear eta
    LF=((1-x(in))^1.75)/((1-alpha)^2);
    
    Mchannel(in)=sqrt((Pin-Pout)*rhofsys/reffT/LF);
    
    alphas(in)=alpha;
    
    clear LF
    
    clear alpha
    end
    
    
end



plot(Qchannel,Mchannel);

%% Determination of system reynolds number

% dynamic viscosity
musys=zeros(1,length(x));

for lm=1:length(x)
    
    if x(lm)<=0
        
        musys(lm)=mulsys;
    else
        musys(lm)=((x(lm)/muvsys)+((1-x(lm))/mulsys))^(-1);
    
    end
end
clear lm

% density 

rhosys=zeros(1,length(alphas));

for lr=1:length(alphas)
    
    if alphas(lr)<=0
        
        rhosys(lr)=rhofsys;
        
    else
        
        rhosys(lr)=(alphas(lr)*rhovsys)+((1-alphas(lr))*rhofsys);
        
    end
end

clear lr

% Reynold's number calculation

Reynolds=Mchannel.*rhosys.*Dh./musys;


%liquid only Reynold's number

ReynoldsLO=zeros(1,length(x));

for Rind=1:length(x)
           
    ReynoldsLO(Rind)=Mchannel(Rind)*rhofsys*Dh/mulsys;
end

%% Volumetric heat generation in fuel

Qvol=Qchannel./Lchannel/pi()/dfuel^2/4*1000;


%% bulk temperature 
% This is actually the highest fluid temperature reached in the channel,
% but as there is negligible axial heat transfer this is where the clad and
% centerline temperatures will be calculated so that the highest possible
% values are found. Tbulk is typically the average temperature but this
% will not give the highest possible value. heat transfer properties will
% be calculated at the maximum possible values for the channel (node)



Tbulk=zeros(1,length(x));

for lt=1:length(x)
    
    if x(lt)<=0
        
        Tbulk(lt)=(Qchannel(lt)*1000/Mchannel(lt)/Cplsys)+Tin;
        
    else Tbulk(lt)=Tsys;
    end
    
    if Tbulk(lt)>Tsys
        Tbulk(lt)=Tsys;
    end
end


hl=zeros(1,length(x));

for lh=1:length(x)
    
    hl(lh)=0.023*ReynoldsLO(lh)^0.8*Prlsys^0.4*klsys/Dh;
    
end

hsys=zeros(1,length(x));

for hsysind=1:length(x)
    
    if x(hsysind)<=0
        
        hsys(hsysind)=hl(hsysind);
        
    else
        
        Co=((1-x(hsysind))/x(hsysind))^0.8*(rhovsys/rhofsys)^0.5;
        
        Fr=(Mchannel(hsysind)/Achannel)^2/(rhofsys^2*9.81*Dh);
        
        if Fr>0.04
            
            C=0;
        else
            C=0.3;
        end
       
        Bo=Qchannel(hsysind)/(pi()*37*doutclad*Lchannel)/(Mchannel(hsysind)/Achannel*1000*hfvsys);
        
        htpconv=hl(hsysind)*((1.1360*Co^-0.9*(25*Fr)^C)+(667.2*Bo^0.7));
        
        htpnuc=hl(hsysind)*((0.6683*Co^-0.2*(25*Fr)^C)+(1058.0*Bo^0.7));
        
        if htpnuc>htpconv
            
            hsys(hsysind)=htpnuc;
            
        else
            
            hsys(hsysind)=htpconv;
        end
    end
end

%% Outer Clad Temperature

Tclado=zeros(1,length(x));

Q=Qchannel.*1000./37;

for Tcladoind=1:length(x)
    
    Tclado(Tcladoind)=(Q(Tcladoind)/(hsys(Tcladoind)*doutclad*Lchannel*pi()))+Tbulk(Tcladoind);
        
end

%% Inner Clad Temperature


Tcladi=zeros(1,length(x));

ric=roc-tclad;

for Tcladind=1:length(x)

    kzirc=(7.51+(0.362e-3*Tclado(Tcladind))-(0.618e-7*Tclado(Tcladind)^2)+(0.718e-11*Tclado(Tcladind)^3))*10^-3;

    Tcladi(Tcladind)=(Q(Tcladind)*log(roc/ric)/(2*pi()*Lchannel*kzirc))+Tclado(Tcladind);
    
end

%% outer fuel meat temperature

Tfuelo=zeros(1,length(x));

dman=ric-(dfuel/2); %m

djump=10e-6; %m



for Tfueloind=1:length(x)
    
    kgap=0.0476+(0.362e-3*Tcladi(Tfueloind))-(0.618e-7*Tcladi(Tfueloind)^2)+(0.718e-11*Tcladi(Tfueloind)^3)*10^-3; %kW/m.C

    hgap=kgap/(dman+djump);
    
    Tfuelo(Tfueloind)=(Q(Tfueloind)/(dfuel/2*Lchannel)/hgap)+Tcladi(Tfueloind);
end

%% Fuel centerline Temperature

Tc=zeros(1,length(x));

rfuel=dfuel/2;

for Tcind=1:length(x)
    
    divf=1000;

    reval=linspace(rfuel,0,divf);

    Tfuel=zeros(1,divf);

    Tfuel(1)=Tfuelo(Tcind);

    for pev=2:divf
    
        Tfuel(pev)=(Qvol(Tcind)/4/kUO2(Tfuel(pev-1))*(reval(pev-1)^2-reval(pev)^2))+Tfuel(pev-1);
    
    end

    Tc(Tcind)=Tfuel(divf);
end


