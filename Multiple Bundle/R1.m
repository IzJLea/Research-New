function res=R1(Tfuel, Tclad,ric,roc,length)

res=0.5*Rfuel(Tfuel,length)+Rgap(Tclad,Tfuel)+0.5*Rzirc(Tclad,ric,roc,length);

end
