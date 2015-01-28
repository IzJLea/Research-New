function res=Alphafunct(x,y)

% requires inputs of Peval,Bundpower,Mflux,Aflow,hvap
M=power/Mflux/Aflow;

hl=XSteam('hL_p',Peval)*1000;


A=x*y+(1-x)*(1-y)-1;

B=hinit-(x*y*hvap)-((1-x)*(1-y)*hl)+((1-x)*M);

res=[A;B];