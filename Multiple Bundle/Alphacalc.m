function res=Alphacalc(Peval,Tvap,Qbundle,hin,Mfluxvapin,Mflux)

if Tvap<=XSteam('Tsat_p',Peval)+0.1
    Tvap=Tvap+0.1;
end

A=Mfluxvapin+(1-alpha)*Qbundle/(alpha*Aflow*(XSteam('hV_p',Peval)-XSteam('hL_p',Peval))*1000)-Mfluxvap;

B=(hout-hin)*Mflux*Aflow-(1-alpha)*Qbundle;

C=(XSteam('h_pT',Peval,Tvap)*alpha+(1-alpha)*XSteam('hL_p',Peval))*1000-hout;



res=solve(A,B,C,alpha,hout,Mfluxvap);


