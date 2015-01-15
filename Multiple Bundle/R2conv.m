function res=R2conv(Tclad,Tcool,ric,roc,Afuel,Dh,Mflux,Peval,length)

res=0.5*Rzirc(Tclad,ric,roc,length)+Rcool(Tcool,Afuel,Dh,Mflux,Peval);

end