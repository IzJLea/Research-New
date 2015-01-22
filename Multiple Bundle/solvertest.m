Peval=110;

Tvap=XSteam('Tsat_p',Peval);

Qbundle=Bundpower(1,5);

hinp=hin(1,5);

Mfluxvapin=0;

if Tvap<=XSteam('Tsat_p',Peval)+0.1
    Tvap=Tvap+0.1;
end
syms alph hout Mfluxvap
A=Mfluxvapin+(1-alph)*Qbundle/(alph*Aflow*(XSteam('hV_p',Peval)-XSteam('hL_p',Peval))*1000)-Mfluxvap;

B=(hout-hinp)*Mflux*Aflow-(1-alph)*Qbundle;

C=(XSteam('h_pT',Peval,Tvap)*alph+(1-alph)*XSteam('hL_p',Peval))*1000-hout;



[alph,hout,Mfluxvap]=solve(A,B,C,alph,hout,Mfluxvap);



