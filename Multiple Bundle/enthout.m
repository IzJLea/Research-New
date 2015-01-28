function res=enthout(alpha,hvap,Peval)

gamma=0.5+1/alpha;

res=(gamma*alpha*hvap)+((1-gamma)*(1-alpha)*1000*XSteam('hL_p',Peval));