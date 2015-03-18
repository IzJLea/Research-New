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

B=qlin/hfg/Mflow;

C=-2*DPT*(alpha-0.5)/L2p;

for i=2:length(L)
    
    L(i,1)=L(i-1,1)+Lbund;
end



F1=@(x) B*L1/2/DPT*(DPT+x(13))-x(1);

F2=@(x) B*Lbund/2/DPT*(x(13)+x(14))-x(2);

F3=@(x) B*Lbund/2/DPT*(x(14)+x(15))-x(3);

F4=@(x) B*Lbund/2/DPT*(x(15)+x(16))-x(4);

F5=@(x) B*Lbund/2/DPT*(x(16)+x(17))-x(5);

F6=@(x) B*Lbund/2/DPT*(x(17)+x(18))-x(6);

F7=@(x) B*Lbund/2/DPT*(x(18)+x(19))-x(7);

F8=@(x) B*Lbund/2/DPT*(x(19)+x(20))-x(8);

F9=@(x) B*Lbund/2/DPT*(x(20)+x(21))-x(9);

F10=@(x) B*Lbund/2/DPT*(x(21)+x(22))-x(10);

F11=@(x) B*Lbund/2/DPT*(x(22)+x(23))-x(11);

F12=@(x) sum(x(1:12))-1;

F13=@(x) C*L(1,1)+DPT-x(13);

F14=@(x) C*L(2,1)+DPT-x(14);

F15=@(x) C*L(3,1)+DPT-x(15);

F16=@(x) C*L(4,1)+DPT-x(16);

F17=@(x) C*L(5,1)+DPT-x(17);

F18=@(x) C*L(6,1)+DPT-x(18);

F19=@(x) C*L(7,1)+DPT-x(19);

F20=@(x) C*L(8,1)+DPT-x(20);

F21=@(x) C*L(9,1)+DPT-x(21);

F22=@(x) C*L(10,1)+DPT-x(22);

F23=@(x) C*L(11,1)+DPT-x(23);

F24=@(x) sum(x(13:24,1))/Kmax-(alpha*Lchannel*L2p);

dF1dx1=-1;

dF1dx13=B*L1/2;

dF2dx2=-1;

dF2dx13=B*L1/2;

dF2dx14=B*L1/2;

dF3dx3=-1;

dF3dx14=B*L1/2;

dF3dx15=B*L1/2;

dF4dx4=-1;

dF4dx15=B*L1/2;

dF4dx16=B*L1/2;

dF5dx5=-1;

dF5dx16=B*L1/2;

dF5dx17=B*L1/2;

dF6dx6=-1;

dF6dx17=B*L1/2;

dF6dx18=B*L1/2;

dF7dx7=-1;

dF7dx18=B*L1/2;

dF7dx19=B*L1/2;

dF8dx8=-1;

dF8dx19=B*L1/2;

dF8dx20=B*L1/2;

dF9dx9=-1;

dF9dx20=B*L1/2;

dF9dx21=B*L1/2;

dF10dx10=-1;

dF10dx21=B*L1/2;

dF10dx22=B*L1/2;

dF11dx11=-1;

dF11dx22=B*L1/2;

dF11dx23=B*L1/2;

dF12dx12=-1;

dF12dx23=B*L1/2;

dF12dx24=B*L1/2;

dF13dx13=-1;

dF14dx14=-1;

dF15dx15=-1;

dF16dx16=-1;

dF17dx17=-1;

dF18dx18=-1;

dF19dx19=-1;

dF20dx20=-1;

dF21dx21=-1;

dF22dx22=-1;

dF23dx23=-1;

dF24dx23=1;

dF24dx22=1;

dF24dx21=1;

dF24dx20=1;

dF24dx19=1;

dF24dx18=1;

dF24dx17=1;

dF24dx16=1;

dF24dx15=1;

dF24dx14=1;

dF24dx13=1;

dF24dx24=1;
%% Creation of Jacobian matrix and Function Matrix

Jac=zeros(2*i,2*i);

F=zeros(2*i,1);

x=ones(2*i,1);

iterations=10000;

iters=zeros(1,iterations);

errs=zeros(1,iterations);

for j=1:iterations
    
    iters(1,j)=j;
    
    Jac(1,1)=dF1dx1;
    
    Jac(1,13)=dF1dx13;
    
    Jac(2,2)=dF2dx2;
    
    Jac(2,13)=dF2dx13;
    
    Jac(2,14)=dF2dx14;
    
    Jac(3,3)=dF3dx3;
    
    Jac(3,14)=dF3dx14;
    
    Jac(3,15)=dF3dx15;
    
    Jac(4,4)=dF4dx4;
    
    Jac(4,15)=dF4dx15;
    
    Jac(4,16)=dF4dx16;
    
    Jac(5,5)=dF5dx5;
    
    Jac(5,16)=dF5dx16;
    
    Jac(5,17)=dF5dx17;
    
    Jac(6,6)=dF6dx6;
    
    Jac(6,17)=dF6dx17;
    
    Jac(6,18)=dF6dx18;
    
    Jac(7,7)=dF7dx7;
    
    Jac(7,18)=dF7dx18;
    
    Jac(7,19)=dF7dx19;
    
    Jac(8,8)=dF8dx8;
    
    Jac(8,19)=dF8dx19;
    
    Jac(8,20)=dF8dx20;
    
    Jac(9,9)=dF9dx9;
    
    Jac(9,20)=dF9dx20;
    
    Jac(9,21)=dF9dx21;
    
    Jac(10,10)=dF10dx10;
    
    Jac(10,21)=dF10dx21;
    
    Jac(10,22)=dF10dx22;
    
    Jac(11,11)=dF11dx11;
    
    Jac(11,22)=dF11dx22;
    
    Jac(11,23)=dF11dx23;
    
    Jac(12,12)=dF12dx12;
    
    Jac(12,23)=dF12dx23;
    
    Jac(12,24)=dF12dx24;
    
    Jac(13,1)=dF13dx1;
    
    Jac(13,13)=dF13dx13;
    
    Jac(14,1)=dF14dx1;
    
    Jac(14,2)=dF14dx2;
    
    Jac(14,14)=dF14dx14;
    
    Jac(15,1)=dF15dx1;
    
    Jac(15,2)=dF15dx2;
    
    Jac(15,3)=dF15dx3;
    
    Jac(15,15)=dF15dx15;
    
    Jac(16,1)=dF16dx1;
    
    Jac(16,2)=dF16dx2;
    
    Jac(16,3)=dF16dx3;
    
    Jac(16,4)=dF16dx4;
    
    Jac(16,16)=dF16dx16;
    
    Jac(17,1)=dF17dx1;
    
    Jac(17,2)=dF17dx2;
    
    Jac(17,3)=dF17dx3;
    
    Jac(17,4)=dF17dx4;
    
    Jac(17,5)=dF17dx5;
    
    Jac(17,17)=dF17dx17;
    
    Jac(18,1)=dF18dx1;
    
    Jac(18,2)=dF18dx2;
    
    Jac(18,3)=dF18dx3;
    
    Jac(18,4)=dF18dx4;
    
    Jac(18,5)=dF18dx5;
    
    Jac(18,6)=dF18dx6;
    
    Jac(18,18)=dF18dx18;
    
    Jac(19,1)=dF19dx1;
    
    Jac(19,2)=dF19dx2;
    
    Jac(19,3)=dF19dx3;
    
    Jac(19,4)=dF19dx4;
    
    Jac(19,5)=dF19dx5;
    
    Jac(19,6)=dF19dx6;
    
    Jac(19,7)=dF19dx7;
    
    Jac(19,19)=dF19dx19;
    
    Jac(20,1)=dF20dx1;
    
    Jac(20,2)=dF20dx2;
    
    Jac(20,3)=dF20dx3;
    
    Jac(20,4)=dF20dx4;
    
    Jac(20,5)=dF20dx5;
    
    Jac(20,6)=dF20dx6;
    
    Jac(20,7)=dF20dx7;
    
    Jac(20,8)=dF20dx8;
    
    Jac(20,20)=dF20dx20;
    
    Jac(21,1)=dF21dx1;
    
    Jac(21,2)=dF21dx2;
    
    Jac(21,3)=dF21dx3;
    
    Jac(21,4)=dF21dx4;
    
    Jac(21,5)=dF21dx5;
    
    Jac(21,6)=dF21dx6;
    
    Jac(21,7)=dF21dx7;
    
    Jac(21,8)=dF21dx8;
    
    Jac(21,9)=dF21dx9;
    
    Jac(21,21)=dF21dx21;
    
    Jac(22,1)=dF22dx1;
    
    Jac(22,2)=dF22dx2;
    
    Jac(22,3)=dF22dx3;
    
    Jac(22,4)=dF22dx4;
    
    Jac(22,5)=dF22dx5;
    
    Jac(22,6)=dF22dx6;
    
    Jac(22,7)=dF22dx7;
    
    Jac(22,8)=dF22dx8;
    
    Jac(22,9)=dF22dx9;
    
    Jac(22,10)=dF22dx10;
    
    Jac(22,22)=dF22dx22;
    
    Jac(23,1)=dF23dx1;
    
    Jac(23,2)=dF23dx2;
    
    Jac(23,3)=dF23dx3;
    
    Jac(23,4)=dF23dx4;
    
    Jac(23,5)=dF23dx5;
    
    Jac(23,6)=dF23dx6;
    
    Jac(23,7)=dF23dx7;
    
    Jac(23,8)=dF23dx8;
    
    Jac(23,9)=dF23dx9;
    
    Jac(23,10)=dF23dx10;
    
    Jac(23,11)=dF23dx11;
    
    Jac(23,23)=dF23dx23;
    
    Jac(24,13)=dF24dx13;
    
    Jac(24,14)=dF24dx14;
    
    Jac(24,15)=dF24dx15;
    
    Jac(24,16)=dF24dx16;
    
    Jac(24,17)=dF24dx17;
    
    Jac(24,18)=dF24dx18;
    
    Jac(24,19)=dF24dx19;
    
    Jac(24,20)=dF24dx20;
    
    Jac(24,21)=dF24dx21;
    
    Jac(24,22)=dF24dx22;
    
    Jac(24,23)=dF24dx23;
    
    Jac(24,24)=dF24dx24;
    
    F(1,1)=F1(x);
    
    F(2,1)=F2(x);
    
    F(3,1)=F3(x);
    
    F(4,1)=F4(x);
    
    F(5,1)=F5(x);
    
    F(6,1)=F6(x);
    
    F(7,1)=F7(x);
    
    F(8,1)=F8(x);
    
    F(9,1)=F9(x);
    
    F(10,1)=F10(x);
    
    F(11,1)=F11(x);
    
    F(12,1)=F12(x);
    
    F(13,1)=F13(x);
    
    F(14,1)=F14(x);
    
    F(15,1)=F15(x);
    
    F(16,1)=F16(x);
    
    F(17,1)=F17(x);
    
    F(18,1)=F18(x);
    
    F(19,1)=F19(x);
    
    F(20,1)=F20(x);
    
    F(21,1)=F21(x);
    
    F(22,1)=F22(x);
    
    F(23,1)=F23(x);
    
    F(24,1)=F24(x);
    
    deltas=-Jac\F;
    
    x=x+(2*deltas);
    
    err=max(abs(deltas));
    
    errs(1,j)=err;
    
    if err<0.0001
        
         display('Converged!');
                 
         %send_mail_message('izaak.lea','CONVERGED!','CONVERGED','R1.m') 
            
         break
    end
    
    plot(iters,errs)
end

%send_mail_message('izaak.lea','NOT CONVERGED!','NOT CONVERGED','R1.m') 
    
display('Not Converged!');
    
    
    