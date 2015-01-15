function res=R3(Tcool,TPT,Aipt,Dh,ript,ropt,length,Mflux,Peval)

res=Rcool(Tcool,Aipt,Dh,Mflux,Peval)+0.5*Rzirc(TPT,ript,ropt,length);

end