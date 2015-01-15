function res=R2(Tclad,Tcool,TPT,Afuel,Arad,Aipt,Dh,Mflux,Peval,ric,roc,length)


res=0.5*Rzirc(Tclad,ric,roc,length)+inv(inv(Rcool(Tcool,Afuel,Dh,Mflux,Peval))+inv(Rrad(Tclad,TPT,Arad,Aipt)));

end
