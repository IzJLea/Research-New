function q=Qloss(Tcool,TPT,TCT,Tmod,ript,ropt,rict,roct,Dh,Lbund,Mflux,Peval)

Aipt=2*pi()*ript*Lbund;

Aoct=2*pi()*roct*Lbund;

Rout=R3(Tcool,TPT,Aipt,Dh,ript,ropt,Lbund,Mflux,Peval)+R4(TPT,TCT,ript,ropt,rict,roct,Lbund)+R5(TCT,Tmod,rict,roct,hmod,Lbund);

q=(Tcool-Tmod)/Rout;