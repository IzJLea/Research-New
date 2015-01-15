function res=R5(TCT,Tmod,rict,roct,Aopt,hmod,length)

res=0.5*Rzirc(TCT,rict,roct,length)+Rmod(Tmod,Aopt,hmod);

end