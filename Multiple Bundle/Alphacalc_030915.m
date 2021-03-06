Qchannel=0.15;
qlin=25;
bundles=12;
Lbund=0.5;
Tenter=300;
Peval=110;
Tsat=318.08;
Lchannel=6;
alpha=0.1155;
Mflow=0.0972;
Aflow=0.0035; %m^2 Flow area in channel
hin=XSteam('hL_T',Tenter);
hlsat=XSteam('hL_p',Peval);
hvsat=XSteam('hV_p',Peval);
hfg=hvsat-hlsat;
Lsat=Mflow*(hlsat-hin)/qlin;
Kmax=(Lchannel-Lsat)/Lbund;
Qboil=qlin*(Lchannel-Lsat);
K=zeros(1,ceil(Kmax));
K(1,1)=rem(Kmax,11);
DPT=0.10338;


L2p=Lchannel-Lsat;

FullBunds=floor(L2p/Lbund);

L1=L2p-(FullBunds*Lbund);

L=zeros(FullBunds+1,1);

L(1,1)=L1;

for i=2:length(L)
    
    L(i,1)=L(i-1,1)+Lbund;
end


