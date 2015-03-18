function res=ALPHACALC(PSH,PVH,Tenter,Tout,Qchannel)
sub=(XSteam('hL_p',Peval)-XSteam('h_pT',Peval,Tenter))/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

wsmx=Qchannel*1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval));

rhoave=(XSteam('rhoV_p',Peval)+XSteam('rho_pT',Peval,Tout))/2;

alpha=0.1;

for i=1:200

    a=alpha^2*rhoave/XSteam('rhoV_p',Peval);
    
    alpha=1/(1+((1+sub)/wsmx*sqrt((PSH-PVH)*rhoave/(RCH+(a*RF)))));
       
   
end

Mflow=alpha*sqrt((PSH-PVH)*rhoave/(RCH+(a*RF)));

res=[alpha;Mflow];


