function Tinit=Tinit_SS(Tvap,Tmod,hvap,hmod)

Lbund=0.5; %m
Qloss=Qloss_single_channel2(Tvap,Tmod,hvap,hmod,Lbund);

Tpti=Tvap-(Qloss(1)*Qloss(2));

Tpto=Tvap-(Qloss(1)*(Qloss(2)+Qloss(3)));

Tcti=Tvap-(Qloss(1)*(Qloss(2)+Qloss(3)+Qloss(4)));

Tcto=Tvap-(Qloss(1)*(Qloss(2)+Qloss(3)+Qloss(4)+Qloss(5)));

Tinit=[Tpti,Tpto,Tcti,Tcto, Qloss(1)];


