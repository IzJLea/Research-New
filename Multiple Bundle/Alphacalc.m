function res=Alphacalc(PSH,PVH,Tenter,Tout,Qchannel,RCH,RF,Hchannel)


Peval=(((PSH+PVH)/2)+(XSteam('rho_pT',PSH/100,Tenter)*9.81*Hchannel/1000))/100; % system evaluation pressure

hin=XSteam('hL_T',Tenter);

rhovap=XSteam('rho_pT',Peval,Tout);

rhoin=XSteam('rho_ph',Peval,hin);

deltaP=(PSH-PVH)+((rhoin-rhovap)*9.81*Hchannel/1000);

rhoave=(rhovap+rhoin)/2;

a=rhoave/rhovap;

x=(XSteam('hL_p',Peval)-hin)/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

wsmx=Qchannel*1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

for i=1:50
    
    alpha=1/(1+((1+x)/wsmx*sqrt(deltaP*rhoave/(RCH+(a*RF)))));
    
    a=alpha^2*rhoave/rhovap;
end

mflow=wsmx*(1+x)/(1-alpha);

res=[alpha;mflow];
