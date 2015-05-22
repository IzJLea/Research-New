%% Steam generator module for full system analysis
function res=Steam_Generator(Hsys,PROH,PRIH,SGfrac,Tvap,div)

%SGfrac is the fraction which the SG is filled. Design is simplified to
%have heat transfer surfaces throughout the entire generator so that
%Atransfer at any fraction fill is SGfrac*Atransfer

% pressure drop information from generic CANDU-9 flow information
% assume 50% of the steam generator secondary side is filled (50% removal
% capacity)

%% maximum heat removal from CANDU-9

Qmax=3347; %MW

Qmaxremoved=Qmax*SGfrac*1000;


%% Height reference

Hheader=11;

HSG=15;

Volume=73.252+33.196; %Combination of the SG volume and the piping volume

%% SG resistance calculation

dPSG=551; % calculated from reference data pressure drops: 200 from ROH to 
          % SGoutlet, 121 from outlet to pump suction, 30 from discharge to 
          % RIH and a guessed 200 from the stopped pump friction losses
rhoref=780.6; % reference density 

mflowref=2667.5; % reference mass flow rate

RSG=dPSG*rhoref/mflowref^2;

%% Main calculations

dP=PROH-PRIH;

Peval=(PRIH+PROH)/2;

rhoV=XSteam('rho_pT',Peval/100,Tvap);

if isnan(rhoV)
    rhoV=XSteam('rhoV_p',PROH/100);
end

rhoL=XSteam('rhoL_p',PRIH/100);

rhoave=(rhoV+rhoL)/2;

if dP>=0
    
    mflow=sqrt(dP*rhoave/RSG);
    
else
    
    mflow=-sqrt(-dP*rhoave/RSG);
end

for i=1:1000
    
    Hsysnew=Hsys-(div*mflow/XSteam('rhoL_p',PRIH/100)*Volume/(HSG-Hheader));

    PRIH=Peval+(9.81*Hsysnew*XSteam('rhoL_p',PRIH/100));

    hin=XSteam('h_pT',PROH/100,Tvap);

    hout=XSteam('hL_p',PRIH/100);

    Qremoved=abs(mflow*(hout-hin));

    e=(Qremoved-Qmaxremoved)/2;
    
    if e<0.1
        
        break
        
    end

    mflow=mflow+(e/(hout-hin));
    
    if i==1000
        converged=0;
    else
        converged=1;
    end
    
    dP=PROH-PRIH;

    Peval=(PRIH+PROH)/2;

    rhoV=XSteam('rho_pT',Peval/100,Tvap);

    if isnan(rhoV)
        rhoV=XSteam('rhoV_p',PROH/100);
    end

    rhoL=XSteam('rhoL_p',PRIH/100);

    rhoave=(rhoV+rhoL)/2;

    if dP>=0
    
        mflow=sqrt(dP*rhoave/RSG);
    
    else
    
        mflow=-sqrt(-dP*rhoave/RSG);
    end
    
end

res=[PRIH;mflow;Hsysnew;converged];
    
    







