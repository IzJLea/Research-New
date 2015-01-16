function res=R2(Tclad,Tcool,TPT,Afuel,Arad,Aipt,Dh,Mflux,Peval,ric,roc,length)

rinv1=inv(Rcool(Tcool,Afuel,Dh,Mflux,Peval));

rinv2=inv(Rrad(Tclad,TPT,Arad,Aipt));


res=0.5*Rzirc(Tclad,ric,roc,length)+inv(rinv1+rinv2);

end
