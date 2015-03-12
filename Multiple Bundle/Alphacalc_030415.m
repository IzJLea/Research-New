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
hin=XSteam('hL_T',Tenter);
hlsat=XSteam('hL_p',Peval);
hvsat=XSteam('hV_p',Peval);
hfg=hvsat-hlsat;
Lsat=Mflow*(hlsat-hin)/qlin;
K=(Lchannel-Lsat)/Lbund;
Qboil=qlin*(Lchannel-Lsat);
Ao=Qboil/Mflow/hfg;
B=Qboil/Mflow/hfg;



%% Functions and derivatives
% functions

F1=@(x) (Ao*x(1))+(B*x(1)^2)-x(2);

F2=@(x) (x(1)*Qboil/hfg/x(2))-1;

F3=@(x) (Ao*x(3))+(B*x(3)^2)-x(4);

F4=@(x) (x(3)*Qboil/hfg/x(4))-1;

F5=@(x) (Ao*x(5))+(B*x(5)^2)-x(6);

F6=@(x) (x(5)*Qboil/hfg/x(6))-1;

F7=@(x) (Ao*x(7))+(B*x(7)^2)-x(8);

F8=@(x) (x(7)*Qboil/hfg/x(8))-1;

F9=@(x) (Ao*x(9))+(B*x(9)^2)-x(10);

F10=@(x) (x(9)*Qboil/hfg/x(10))-1;

F11=@(x) (Ao*x(11))+(B*x(11)^2)-x(12);

F12=@(x) (x(11)*Qboil/hfg/x(12))-1;

F13=@(x) (Ao*x(13))+(B*x(13)^2)-x(14);

F14=@(x) (x(13)*Qboil/hfg/x(14))-1;

F15=@(x) (Ao*x(15))+(B*x(15)^2)-x(16);

F16=@(x) (x(15)*Qboil/hfg/x(16))-1;

F17=@(x) (Ao*x(17))+(B*x(17)^2)-x(18);

F18=@(x) (x(17)*Qboil/hfg/x(18))-1;

F19=@(x) (Ao*x(19))+(B*x(19)^2)-x(20);

F20=@(x) (x(19)*Qboil/hfg/x(20))-1;

F21=@(x) (Ao*x(21))+(B*x(21)^2)-x(22);

F22=@(x) (x(21)*Qboil/hfg/x(22))-1;

F23=@(x) (Ao*x(23))+(B*x(23)^2)-x(24);

F24=@(x) (x(23)*Qboil/hfg/x(24))-1;

dF1dx1=@(x) Ao+(B*x(1));

dF1dx2=@(x)-1;

dF2dx1=@(x) Qboil/x(2)/hfg;

dF2dx2=@(x) -Qboil*x(1)/x(2)^2/hfg;

dF3dx3=@(x) Ao+(B*x(3));

dF3dx4=@(x)-1;

dF4dx3=@(x) Qboil/x(4)/hfg;

dF4dx4=@(x) -Qboil*x(3)/x(4)^2/hfg;

dF5dx5=@(x) Ao+(B*x(5));

dF5dx6=@(x)-1;

dF6dx5=@(x) Qboil/x(6)/hfg;

dF6dx6=@(x) -Qboil*x(5)/x(6)^2/hfg;

dF7dx7=@(x) Ao+(B*x(7));

dF7dx8=@(x)-1;

dF8dx7=@(x) Qboil/x(8)/hfg;

dF8dx8=@(x) -Qboil*x(7)/x(8)^2/hfg;

dF9dx9=@(x) Ao+(B*x(9));

dF9dx10=@(x)-1;

dF10dx9=@(x) Qboil/x(10)/hfg;

dF10dx10=@(x) -Qboil*x(9)/x(10)^2/hfg;

dF11dx11=@(x) Ao+(B*x(10));

dF11dx12=@(x)-1;

dF12dx11=@(x) Qboil/x(12)/hfg;

dF12dx12=@(x) -Qboil*x(11)/x(12)^2/hfg;

dF13dx13=@(x) Ao+(B*x(13));

dF13dx14=@(x)-1;

dF14dx13=@(x) Qboil/x(14)/hfg;

dF14dx14=@(x) -Qboil*x(13)/x(14)^2/hfg;

dF15dx15=@(x) Ao+(B*x(15));

dF15dx16=@(x)-1;

dF16dx15=@(x) Qboil/x(16)/hfg;

dF16dx16=@(x) -Qboil*x(15)/x(16)^2/hfg;

dF17dx17=@(x) Ao+(B*x(17));

dF17dx18=@(x)-1;

dF18dx17=@(x) Qboil/x(18)/hfg;

dF18dx18=@(x) -Qboil*x(17)/x(18)^2/hfg;

dF19dx19=@(x) Ao+(B*x(19));

dF19dx20=@(x)-1;

dF20dx19=@(x) Qboil/x(20)/hfg;

dF20dx20=@(x) -Qboil*x(19)/x(20)^2/hfg;

dF21dx21=@(x) Ao+(B*x(21));

dF21dx22=@(x)-1;

dF22dx21=@(x) Qboil/x(22)/hfg;

dF22dx22=@(x) -Qboil*x(21)/x(22)^2/hfg;

dF23dx23=@(x) Ao+(B*x(23));

dF23dx24=@(x)-1;

dF24dx23=@(x) Qboil/x(24)/hfg;

dF24dx24=@(x) -Qboil*x(23)/x(24)^2/hfg;

%% Creation of Jacobian and solving matrices

% x matrix and initial guess

x=ones(24,1)*alpha;

Jac=zeros(24,24);

F=zeros(24,1);

for j=1:24
    
    if mod(j,2)
        
        x(j)=alpha;
    else
        x(j)=Mflow/K;
    end
end



for i=1:1
    
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
    
    Jac(1,1)=dF1dx1(x);
    
    Jac(1,2)=dF1dx2(x);
    
    Jac(2,1)=dF2dx1(x);
    
    Jac(2,2)=dF2dx2(x);
    
    Jac(3,3)=dF3dx3(x);
    
    Jac(3,4)=dF3dx4(x);
    
    Jac(4,3)=dF4dx3(x);
    
    Jac(4,4)=dF4dx4(x);
    
    Jac(5,5)=dF5dx5(x);
    
    Jac(5,6)=dF5dx6(x);
    
    Jac(6,5)=dF6dx5(x);
    
    Jac(6,6)=dF6dx6(x);
    
    Jac(7,7)=dF7dx7(x);
    
    Jac(7,8)=dF7dx8(x);

    Jac(8,7)=dF8dx8(x);
    
    Jac(8,8)=dF8dx8(x);
    
    Jac(9,9)=dF9dx9(x);
    
    Jac(9,10)=dF9dx10(x);
    
    Jac(10,9)=dF10dx9(x);
    
    Jac(10,10)=dF10dx10(x);
    
    Jac(11,11)=dF11dx11(x);
    
    Jac(11,12)=dF11dx12(x);
    
    Jac(12,11)=dF12dx11(x);
    
    Jac(12,12)=dF12dx12(x);
    
    Jac(13,13)=dF13dx13(x);
    
    Jac(13,14)=dF13dx14(x);
    
    Jac(14,14)=dF14dx14(x);
    
    Jac(15,15)=dF15dx15(x);
    
    Jac(15,16)=dF15dx16(x);
    
    Jac(16,15)=dF16dx15(x);

    Jac(16,16)=dF16dx16(x);
    
    Jac(17,17)=dF17dx17(x);
    
    Jac(17,18)=dF17dx18(x);
    
    Jac(18,17)=dF18dx17(x);
    
    Jac(18,18)=dF18dx18(x);
    
    Jac(19,19)=dF19dx19(x);
    
    Jac(19,20)=dF19dx20(x);
    
    Jac(20,19)=dF20dx19(x);
    
    Jac(20,20)=dF20dx20(x);
    
    Jac(21,21)=dF21dx21(x);
    
    Jac(21,22)=dF21dx22(x);
    
    Jac(22,21)=dF22dx21(x);
    
    Jac(22,22)=dF22dx22(x);
    
    Jac(23,23)=dF23dx23(x);
    
    Jac(23,24)=dF23dx24(x);
    
    Jac(24,23)=dF24dx23(x);
    
    Jac(24,24)=dF24dx24(x);
    
    deltas=F\Jac;
    
    deltas1=rot90(deltas,-1);
    
    
end
