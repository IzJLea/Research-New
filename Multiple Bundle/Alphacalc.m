function res=Alphacalc(power,Aflow,Peval,Mflux,gammain,alphain,hvapin,hvap)

syms x y
A=(2*x*y)-x-y;

B=(gammain*alphain*Mflux*Aflow*hvapin)-(y*x*Mflux*Aflow*hvap)+((1-x)*power/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval))*XSteam('hV_p',Peval));

funct=[A;B];

res=solve(funct,0.5,0.5);







