function res=R5(TCT,rict,roct,hmod,length)

Aoct=2*pi()*roct*length;

res=0.5*Rzirc(TCT,rict,roct,length)+Rmod(hmod,Aoct);

end