function res=R2rad(Tclad,TPT,ric,roc,Arad,Aipt,length)

res=0.5*Rzirc(Tclad,ric,roc,length)+Rrad(Tclad,TPT,Arad,Aipt);

end