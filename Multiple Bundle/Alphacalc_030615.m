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


Ms=Mflow/Aflow; %Mass flux

l=ceil(Kmax);

for k=2:l
    
    K(1,k)=K(1,k-1)+1;
end

L2p=K*0.5;

L1=L2p(1,1);

L2=L2p(1,11)-L1; % Length of 'horizontal' section

L3=L2p(1,12)-L1-L2; %single bundle length

Ltotal=L2p(1,12);

alpha1=0.5*alpha;

alpha2=alpha;

alpha3=((alpha*Ltotal)-(L1*alpha1)-(L2*alpha))/L3;

%% Functions and derivatives

%initial Functions

F1=@(x) (qlin*L1*x(13)/hfg)-x(1);

F2=@(x) (qlin*L3*x(14)/hfg)-x(2);

F3=@(x) (qlin*L3*x(14)/hfg)-x(3);

F4=@(x) (qlin*L3*x(14)/hfg)-x(4);

F5=@(x) (qlin*L3*x(14)/hfg)-x(5);

F6=@(x) (qlin*L3*x(14)/hfg)-x(6);

F7=@(x) (qlin*L3*x(14)/hfg)-x(7);

F8=@(x) (qlin*L3*x(14)/hfg)-x(8);

F9=@(x) (qlin*L3*x(14)/hfg)-x(9);

F10=@(x) (qlin*L3*x(14)/hfg)-x(10);

F11=@(x) (qlin*L3*x(14)/hfg)-x(11);

F12=@(x) (qlin*L3*x(15)/hfg)-x(12);

F13=@(x) x(15)+(x(14)/2)-1;

F14=@(x) (x(13)*L1/Ltotal)+(x(14)*L2/Ltotal)+(x(15)*L3/Ltotal)-alpha;

F15=@(x) (x(14)/2)-x(13)+0.5;

% derivatives

dF1dx1=-1;

dF1dx13=qlin*L1/hfg;

dF2dx2=-1;

dF2dx14=qlin*L3/hfg;

dF3dx3=-1;

dF3dx14=qlin*L3/hfg;

dF4dx4=-1;

dF4dx14=qlin*L3/hfg;

dF5dx5=-1;

dF5dx14=qlin*L3/hfg;

dF6dx6=-1;

dF6dx14=qlin*L3/hfg;

dF7dx7=-1;

dF7dx14=qlin*L3/hfg;

dF8dx8=-1;

dF8dx14=qlin*L3/hfg;

dF9dx9=-1;

dF9dx14=qlin*L3/hfg;

dF10dx10=-1;

dF10dx14=qlin*L3/hfg;

dF11dx11=-1;

dF11dx14=qlin*L3/hfg;

dF12dx12=-1;

dF12dx15=qlin*L3/hfg;

dF13dx14=1/2;

dF13dx15=1;

dF14dx13=L1/Ltotal;

dF14dx14=L2/Ltotal;

dF14dx15=L3/Ltotal;

dF15dx13=-1;

dF15dx14=1/2;

%% matrix creation

x=zeros(15,1);

for i=1:12
    
    x(i,1)=Mflow/12;
    
end

x(13,1)=alpha1;

x(14,1)=alpha;

x(15,1)=alpha3;

F=zeros(15,1);

Jac=zeros(15,15);

for j=1:100000
    
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
    
    Jac(1,1)=dF1dx1;
    
    Jac(1,13)=dF1dx13;
    
    Jac(2,2)=dF2dx2;
    
    Jac(2,14)=dF2dx14;
    
    Jac(3,3)=dF3dx3;
    
    Jac(3,14)=dF3dx14;
    
    Jac(4,4)=dF4dx4;
    
    Jac(4,14)=dF4dx14;
    
    Jac(5,5)=dF5dx5;
    
    Jac(5,14)=dF5dx14;
    
    Jac(6,6)=dF6dx6;
    
    Jac(6,14)=dF6dx14;
    
    Jac(7,7)=dF7dx7;
    
    Jac(7,14)=dF7dx14;
    
    Jac(8,8)=dF8dx8;
    
    Jac(8,14)=dF8dx14;
    
    Jac(9,9)=dF9dx9;
    
    Jac(9,14)=dF9dx14;
    
    Jac(10,10)=dF10dx10;
    
    Jac(10,14)=dF10dx14;
    
    Jac(11,11)=dF11dx11;
    
    Jac(11,14)=dF11dx14;
    
    Jac(12,12)=dF12dx12;
    
    Jac(12,15)=dF12dx15;
    
    Jac(13,14)=dF13dx14;
    
    Jac(13,15)=dF13dx15;
       
    Jac(14,13)=dF14dx13;
    
    Jac(14,14)=dF14dx14;
    
    Jac(14,15)=dF14dx15;
    
    Jac(15,13)=dF15dx13;
    
    Jac(15,14)=dF15dx14;
    
    deltas=Jac\-F;
    
    x=x+deltas;
    
    err=max(abs(deltas));
    
    if err<1e-20
        
        display('Converged!', num2str(j));
                 
        send_mail_message('izaak.lea','CONVERGED!','CONVERGED','R1.m') 
          
        break
    end
    
    
    display('Count:', num2str(j));
end
    
    
send_mail_message('izaak.lea','NOT CONVERGED!','NOT CONVERGED','R1.m') 


