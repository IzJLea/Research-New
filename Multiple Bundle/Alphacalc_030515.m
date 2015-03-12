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

Ms=Mflow/Aflow; %Mass flux

l=ceil(Kmax);

for k=2:l
    
    K(1,k)=K(1,k-1)+1;
end


%% Functions and derivatives
% functions

F1=@(x) Ms*(x(1)-1)-x(2)*(1-2*x(1));

F2=@(x) x(2)*Ms/Qboil-x(1);

F3=@(x) Ms*(x(3)-1)-(x(2)+x(4))*(1-2*x(3));

F4=@(x) x(4)*Ms/Qboil-x(3);

F5=@(x) Ms*(x(5)-1)-(x(2)+x(4)+x(6))*(1-2*x(5));

F6=@(x) x(6)*Ms/Qboil-x(5);

F7=@(x) Ms*(x(7)-1)-(x(2)+x(4)+x(6)+x(8))*(1-2*x(7));

F8=@(x) x(8)*Ms/Qboil-x(7);

F9=@(x) Ms*(x(9)-1)-(x(2)+x(4)+x(6)+x(8)+x(10))*(1-2*x(9));

F10=@(x) x(10)*Ms/Qboil-x(9);

F11=@(x) Ms*(x(11)-1)-(x(2)+x(4)+x(6)+x(8)+x(10)+x(12))*(1-2*x(11));

F12=@(x) x(12)*Ms/Qboil-x(11);

F13=@(x) Ms*(x(13)-1)-(x(2)+x(4)+x(6)+x(8)+x(10)+x(12)+x(14))*(1-2*x(13));

F14=@(x) x(14)*Ms/Qboil-x(13);

F15=@(x) Ms*(x(15)-1)-(x(2)+x(4)+x(6)+x(8)+x(10)+x(12)+x(14)+x(16))*(1-2*x(15));

F16=@(x) x(16)*Ms/Qboil-x(15);

F17=@(x) Ms*(x(17)-1)-(x(2)+x(4)+x(6)+x(8)+x(10)+x(12)+x(14)+x(16)+x(18))*(1-2*x(17));

F18=@(x) x(18)*Ms/Qboil-x(17);

F19=@(x) Ms*(x(19)-1)-(x(2)+x(4)+x(6)+x(8)+x(10)+x(12)+x(14)+x(16)+x(18)+x(20))*(1-2*x(19));

F20=@(x) x(20)*Ms/Qboil-x(19);

F21=@(x) Ms*(x(21)-1)-(x(2)+x(4)+x(6)+x(8)+x(10)+x(12)+x(14)+x(16)+x(18)+x(20)+x(22))*(1-2*x(21));

F22=@(x) x(22)*Ms/Qboil-x(21);

F23=@(x) Ms*(x(23)-1)-(x(2)+x(4)+x(6)+x(8)+x(10)+x(12)+x(14)+x(16)+x(18)+x(20)+x(22)+x(24))*(1-2*x(23));

F24=@(x) x(24)*Ms/Qboil-x(23);

dF1dx1=@(x) Ms+2*x(2);

dF1dx2=@(x) 2*x(1)-1;

dF2dx1=@(x) -1;

dF2dx2=@(x) Ms/Qboil;

dF3dx3=@(x) Ms+(2*(x(2)+x(4)));

dF3dx4=@(x) 2*x(3)-1;

dF4dx3=@(x) -1;

dF4dx4=@(x) Ms/Qboil;

dF5dx5=@(x) Ms+(2*(x(2)+x(4)+x(6)));

dF5dx6=@(x) 2*x(5)-1;

dF6dx5=@(x) -1;

dF6dx6=@(x) Ms/Qboil;

dF7dx7=@(x) Ms+(2*(x(2)+x(4)+x(6)+x(8)));

dF7dx8=@(x) 2*x(7)-1;

dF8dx7=@(x) -1;

dF8dx8=@(x) Ms/Qboil;

dF9dx9=@(x) Ms+(2*(x(2)+x(4)+x(6)+x(8)+x(10)));

dF9dx10=@(x) 2*x(9)-1;

dF10dx9=@(x) -1;

dF10dx10=@(x) Ms/Qboil;

dF11dx11=@(x) Ms+(2*(x(2)+x(4)+x(6)+x(8)+x(10)+x(12)));

dF11dx12=@(x) 2*x(11)-1;

dF12dx11=@(x) -1;

dF12dx12=@(x) Ms/Qboil;

dF13dx13=@(x) Ms+(2*(x(2)+x(4)+x(6)+x(8)+x(10)+x(12)+x(14)));

dF13dx14=@(x) 2*x(13)-1;

dF14dx13=@(x) -1;

dF14dx14=@(x) Ms/Qboil;

dF15dx15=@(x) Ms+(2*(x(2)+x(4)+x(6)+x(8)+x(10)+x(12)+x(14)+x(16)));

dF15dx16=@(x) 2*x(15)-1;

dF16dx15=@(x) -1;

dF16dx16=@(x) Ms/Qboil;

dF17dx17=@(x) Ms+(2*(x(2)+x(4)+x(6)+x(8)+x(10)+x(12)+x(14)+x(16)+x(18)));

dF17dx18=@(x) 2*x(17)-1;

dF18dx17=@(x) -1;

dF18dx18=@(x) Ms/Qboil;

dF19dx19=@(x) Ms+(2*(x(2)+x(4)+x(6)+x(8)+x(10)+x(12)+x(14)+x(16)+x(18)+x(20)));

dF19dx20=@(x) 2*x(19)-1;

dF20dx19=@(x) -1;

dF20dx20=@(x) Ms/Qboil;

dF21dx21=@(x) Ms+(2*(x(2)+x(4)+x(6)+x(8)+x(10)+x(12)+x(14)+x(16)+x(18)+x(20)+x(22)));

dF21dx22=@(x) 2*x(21)-1;

dF22dx21=@(x) -1;

dF22dx22=@(x) Ms/Qboil;

dF23dx23=@(x) Ms+(2*(x(2)+x(4)+x(6)+x(8)+x(10)+x(12)+x(14)+x(16)+x(18)+x(20)+x(22)+x(24)));

dF23dx24=@(x) 2*x(23)-1;

dF24dx23=@(x) -1;

dF24dx24=@(x) Ms/Qboil;

Jac=zeros(2*ceil(Kmax),2*ceil(Kmax));

F=zeros(2*ceil(Kmax),1);

x=zeros(2*cel(Kmax),1);

for i=1:2*ceil(Kmax)
    
    if mod(i,2)
        
        x(i)=alpha;
    else
        
        x(i)=a

for i=1:1
    
    F(1,1)=F1(x)
