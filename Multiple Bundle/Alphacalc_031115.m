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

x=zeros(length(L)+1,1);

F1=@(x) 1-(1/2/DPT*(x(1)+x(2)));

F2=@(x) 1-(1/2/DPT*(x(2)+x(3)));

F3=@(x) 1-(1/2/DPT*(x(3)+x(4)));

F4=@(x) 1-(1/2/DPT*(x(4)+x(5)));

F5=@(x) 1-(1/2/DPT*(x(5)+x(6)));

F6=@(x) 1-(1/2/DPT*(x(6)+x(7)));

F7=@(x) 1-(1/2/DPT*(x(7)+x(8)));

F8=@(x) 1-(1/2/DPT*(x(8)+x(9)));

F9=@(x) 1-(1/2/DPT*(x(9)+x(10)));

F10=@(x) 1-(1/2/DPT*(x(10)+x(11)));

F11=@(x) 1-(1/2/DPT*(x(11)+x(12)));

F12=@(x) 1-(1/2/DPT*(x(12)+x(13)));

F13=@(x) Qchannel*1000/hfg*((L2p/Lbund)-(L1/L2p*x(1)/DPT)-((L1/L2p+1)*x(2)/2/DPT)-x(3)/DPT-x(4)/DPT-x(5)/DPT-x(6)/DPT-x(7)/DPT-x(8)/DPT-x(9)/DPT-x(10)/DPT-x(11)/DPT-x(12)/DPT-x(13)/DPT)-Mflow;

dF1dx1=-1/2/DPT;

dF1dx2=-1/2/DPT;

dF2dx2=-1/2/DPT;

dF2dx3=-1/2/DPT;

dF3dx3=-1/2/DPT;

dF3dx4=-1/2/DPT;

dF4dx4=-1/2/DPT;

dF4dx5=-1/2/DPT;

dF5dx5=-1/2/DPT;

dF5dx6=-1/2/DPT;

dF6dx6=-1/2/DPT;

dF6dx7=-1/2/DPT;

dF7dx7=-1/2/DPT;

dF7dx8=-1/2/DPT;

dF8dx8=-1/2/DPT;

dF8dx9=-1/2/DPT;

dF9dx9=-1/2/DPT;

dF9dx10=-1/2/DPT;

dF10dx10=-1/2/DPT;

dF10dx11=-1/2/DPT;

dF11dx11=-1/2/DPT;

dF11dx12=-1/2/DPT;

dF12dx12=-1/2/DPT;

dF12dx13=-1/2/DPT;

dF13dx1=-Qchannel*1000/hfg*L1/L2p/DPT;

dF13dx2=-Qchannel*1000/hfg*(L1/L2p+1)/2/DPT;

dF13dx3=-Qchannel*1000/hfg/DPT;

dF13dx4=-Qchannel*1000/hfg/DPT;

dF13dx5=-Qchannel*1000/hfg/DPT;

dF13dx6=-Qchannel*1000/hfg/DPT;

dF13dx7=-Qchannel*1000/hfg/DPT;

dF13dx8=-Qchannel*1000/hfg/DPT;

dF13dx9=-Qchannel*1000/hfg/DPT;

dF13dx10=-Qchannel*1000/hfg/DPT;

dF13dx11=-Qchannel*1000/hfg/DPT;

dF13dx12=-Qchannel*1000/hfg/DPT;

dF13dx13=-Qchannel*1000/hfg/DPT/2;


%% Creation of Jacobian matrix and Function Matrix

Jac=zeros(length(L)+1,length(L)+1);

F=ones(length(L)+1,1);

iters=100000;

err=zeros(length(x),iters);

abserr=zeros(length(x),iters);

for j=1:iters
    
    Jac(1,1)=dF1dx1;
    
    Jac(1,2)=dF1dx2;
    
    Jac(2,2)=dF2dx2;
    
    Jac(2,3)=dF2dx3;
    
    Jac(3,3)=dF3dx3;
    
    Jac(3,4)=dF3dx4;
    
    Jac(4,4)=dF4dx4;
    
    Jac(4,5)=dF4dx5;
    
    Jac(5,5)=dF5dx5;
    
    Jac(5,6)=dF5dx6;
    
    Jac(6,6)=dF6dx6;
    
    Jac(6,7)=dF6dx7;
    
    Jac(7,7)=dF7dx7;
    
    Jac(7,8)=dF7dx8;
    
    Jac(8,8)=dF8dx8;
    
    Jac(8,9)=dF8dx9;
    
    Jac(9,9)=dF9dx9;
    
    Jac(9,10)=dF9dx10;
    
    Jac(10,10)=dF10dx10;
    
    Jac(10,11)=dF10dx11;
    
    Jac(11,11)=dF11dx11;
    
    Jac(11,12)=dF11dx12;
    
    Jac(12,12)=dF12dx12;
    
    Jac(12,13)=dF12dx13;
    
    Jac(13,1)=dF13dx1;
    
    Jac(13,2)=dF13dx2;
    
    Jac(13,3)=dF13dx3;
    
    Jac(13,4)=dF13dx4;
    
    Jac(13,5)=dF13dx5;
    
    Jac(13,6)=dF13dx6;
    
    Jac(13,7)=dF13dx7;
    
    Jac(13,8)=dF13dx8;
    
    Jac(13,9)=dF13dx9;
    
    Jac(13,10)=dF13dx10;
    
    Jac(13,11)=dF13dx11;
    
    Jac(13,12)=dF13dx12;
    
    Jac(13,13)=dF13dx13;
    
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
    
    deltas=-Jac\F;
    
    xnew=x+0.01*(deltas);
    
    err(1:length(x),j)=(xnew-x)./2;
    
    x=xnew;
    
   
    
    abserr(1:length(x),j)=abs(err(1:length(x),j));
    
    if abserr(1,j)<0.0001
        
        display('Converged');
        
        break
    end
end

if abserr(1,j)>0.0001
    
    display('Not Converged');
end
    





    