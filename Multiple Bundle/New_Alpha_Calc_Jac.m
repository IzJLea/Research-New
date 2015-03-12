
Tenter=300; % Entering coolant temperature

Peval=110; % Evaluation pressure

mflow=0.097243680263615; % kg/s mass flow rate

Boil=3; % initial boiling position (12-Boil+1)*2=# of equations (20)

Functions=12-Boil+1;

Qchannel=0.150; % MW 

Bpower=cosPower(12,Qchannel*1000000); %Per bundle power values

mmax=Bpower/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval))/1000;

hin=1.445e+06; % Entering enthalpy

%%
% Main Function creation;

F1=@(x,y) -2*x(1)*y(1)+x(1)+y(1);

F2=@(x,y) -x(1)*y(1)*mflow+((1-x(1))*Bpower(1,Boil))/((XSteam('hV_p',Peval)*1000)-hin);

F3=@(x,y) -2*x(2)*y(2)+x(2)+y(2);

F4=@(x,y) x(1)*y(1)*mflow-x(2)*y(2)*mflow+(1-x(2))*mmax(1,Boil+1);

F5=@(x,y) -2*x(3)*y(3)+x(3)+y(3);

F6=@(x,y) x(2)*y(2)*mflow-x(3)*y(3)*mflow+(1-x(3))*mmax(1,Boil+2);

F7=@(x,y) -2*x(4)*y(4)+x(4)+y(4);

F8=@(x,y) x(3)*y(3)*mflow-x(4)*y(4)*mflow+(1-x(4))*mmax(1,Boil+3);

F9=@(x,y) -2*x(5)*y(5)+x(5)+y(5);

F10=@(x,y) x(4)*y(4)*mflow-x(5)*y(5)*mflow+(1-x(5))*mmax(1,Boil+4);

F11=@(x,y) -2*x(6)*y(6)+x(6)+y(6);

F12=@(x,y) x(5)*y(5)*mflow-x(6)*y(6)*mflow+(1-x(6))*mmax(1,Boil+5);

F13=@(x,y) -2*x(7)*y(7)+x(7)+y(7);

F14=@(x,y) x(6)*y(6)*mflow-x(7)*y(7)*mflow+(1-x(7))*mmax(1,Boil+6);

F15=@(x,y) -2*x(8)*y(8)+x(8)+y(8);

F16=@(x,y) x(7)*y(7)*mflow-x(8)*y(8)*mflow+(1-x(8))*mmax(1,Boil+7);

F17=@(x,y) -2*x(9)*y(9)+x(9)+y(9);

F18=@(x,y) x(8)*y(8)*mflow-x(9)*y(9)*mflow+(1-x(9))*mmax(1,Boil+8);

F19=@(x,y) 1-x(10)*y(10);

F20=@(x,y) x(9)*y(9)*mflow-x(10)*y(10)*mflow+(1-x(10))*mmax(1,Boil+9);

%%
%Derivatives for Jacobian

dF1dx1=@(y) -2*y(1)+1;
 
dF1dy1=@(x) -2*x(1)+1;

dF2dx1=@(y) -y(1)*mflow-mmax(1,Boil);

dF2dy1=@(x) -x(1)*mflow;

dF3dx1=@(y) 0;

dF3dy1=@(x) 0;

dF3dx2=@(y) -2*y(2)+1;

dF3dy2=@(x) -2*x(2)+1;

dF4dx1=@(y) y(1)*mflow;

dF4dy1=@(x) x(1)*mflow;

dF4dx2=@(y) -y(2)*mflow-mmax(1,Boil+1);

dF4dy2=@(x) -x(2)*mflow;

dF5dx2=@(y) 0;

dF5dy2=@(x) 0;

dF5dx3=@(y) -2*y(3)+1;

dF5dy3=@(x) -2*x(3)+1;

dF6dx2=@(y) y(2)*mflow;

dF6dy2=@(x) x(2)*mflow;

dF6dx3=@(y) -y(3)*mflow-mmax(1,Boil+2);

dF6dy3=@(x) -x(3)*mflow;

dF7dx3=@(y) 0;

dF7dy3=@(x) 0;

dF7dx4=@(y) -2*y(4)+1;

dF7dy4=@(x) -2*x(4)+1;

dF8dx3=@(y) y(3)*mflow;

dF8dy3=@(x) x(3)*mflow;

dF8dx4=@(y) -y(4)*mflow-mmax(1,Boil+3);

dF8dy4=@(x) -x(4)*mflow;

dF9dx4=@(y) 0;

dF9dy4=@(x) 0;

dF9dx5=@(y) -2*y(5)+1;

dF9dy5=@(x) -2*x(5)+1;

dF10dx4=@(y) y(4)*mflow;

dF10dy4=@(x) x(4)*mflow;

dF10dx5=@(y) -y(5)*mflow-mmax(1,Boil+4);

dF10dy5=@(x) -x(5)*mflow;

dF11dx5=@(y) 0;

dF11dy5=@(x) 0;

dF11dx6=@(y) -2*y(6)+1;

dF11dy6=@(x) -2*x(6)+1;

dF12dx5=@(y) y(5)*mflow;

dF12dy5=@(x) x(5)*mflow;

dF12dx6=@(y) -y(6)*mflow-mmax(1,Boil+5);

dF12dy6=@(x) -x(6)*mflow;

dF13dx6=@(y) 0;

dF13dy6=@(x) 0;

dF13dx7=@(y) -2*y(7)+1;

dF13dy7=@(x) -2*x(7)+1;

dF14dx6=@(y) y(6)*mflow;

dF14dy6=@(x) x(6)*mflow;

dF14dx7=@(y) -y(7)*mflow-mmax(1,Boil+6);

dF14dy7=@(x) -x(7)*mflow;

dF15dx7=@(y) 0;

dF15dy7=@(x) 0;

dF15dx8=@(y) -2*y(8)+1;

dF15dy8=@(x) -2*x(8)+1;

dF16dx7=@(y) y(7)*mflow;

dF16dy7=@(x) x(7)*mflow;

dF16dx8=@(y) -y(8)*mflow-mmax(1,Boil+7);

dF16dy8=@(x) -x(8)*mflow;

dF17dx8=@(y) 0;

dF17dy8=@(x) 0;

dF17dx9=@(y) -2*y(9)+1;

dF17dy9=@(x) -2*x(9)+1;

dF18dx8=@(y) y(8)*mflow;

dF18dy8=@(x) x(8)*mflow;

dF18dx9=@(y) -y(9)*mflow-mmax(1,Boil+8);

dF18dy9=@(x) -x(9)*mflow;

dF19dx10=@(y) -y(10);

dF19dy10=@(x) -x(10);

dF20dx9=@(y) y(9)*mflow;

dF20dy9=@(x) x(9)*mflow;

dF20dx10=@(y) -y(10)*mflow-mmax(1,Boil+9);

dF20dy10=@(x) -x(10)*mflow;

%%
%Solving system

x=ones(12-Boil+1,1)*0.1;

y=ones(12-Boil+1,1)*0.9;

Iterations=200000000;

for iter=1:Iterations
   
    
    Jac=zeros(Functions,Functions);
    
    F=zeros(Functions,1);
    
    dx=zeros(Functions,1);
    
    dy=zeros(Functions,1);
    
    Jac(1,1)=dF1dx1(y);
    
    Jac(1,2)=dF1dy1(x);
    
    Jac(2,1)=dF2dx1(y);
    
    Jac(2,2)=dF2dy1(x);
    
    Jac(3,1)=dF3dx1(y);
    
    Jac(3,2)=dF3dy1(x);
    
    Jac(3,3)=dF3dx2(y);
    
    Jac(3,4)=dF3dy2(x);
    
    Jac(4,1)=dF4dx1(y);
    
    Jac(4,2)=dF4dy1(x);
    
    Jac(4,3)=dF4dx2(y);
        
    Jac(4,4)=dF4dy2(y);
    
    Jac(5,3)=dF5dx2(y);
    
    Jac(5,4)=dF5dy2(x);
    
    Jac(5,5)=dF5dx3(y);
    
    Jac(5,6)=dF5dy3(x);
    
    Jac(6,3)=dF6dx2(y);
    
    Jac(6,4)=dF6dy2(x);
    
    Jac(6,5)=dF6dx3(y);
        
    Jac(6,6)=dF6dy3(y);
    
    Jac(7,5)=dF7dx3(y);
    
    Jac(7,6)=dF7dy3(x);
    
    Jac(7,7)=dF7dx4(y);
    
    Jac(7,8)=dF7dy4(x);
    
    Jac(8,5)=dF8dx3(y);
    
    Jac(8,6)=dF8dy3(x);
    
    Jac(8,7)=dF8dx4(y);
        
    Jac(8,8)=dF8dy4(y);
    
    Jac(9,7)=dF9dx4(y);
    
    Jac(9,8)=dF9dy4(x);
    
    Jac(9,9)=dF9dx5(y);
    
    Jac(9,10)=dF9dy5(x);
    
    Jac(10,7)=dF10dx4(y);
    
    Jac(10,8)=dF10dy4(x);
    
    Jac(10,9)=dF10dx5(y);
        
    Jac(10,10)=dF10dy5(y);
    
    Jac(11,9)=dF11dx5(y);
    
    Jac(11,10)=dF11dy5(x);
    
    Jac(11,11)=dF11dx6(y);
    
    Jac(11,12)=dF11dy6(x);
    
    Jac(12,9)=dF12dx5(y);
    
    Jac(12,10)=dF12dy5(x);
    
    Jac(12,11)=dF12dx6(y);
        
    Jac(12,12)=dF12dy6(y);
    
    Jac(13,11)=dF13dx6(y);
    
    Jac(13,12)=dF13dy6(x);
    
    Jac(13,13)=dF13dx7(y);
    
    Jac(13,14)=dF13dy7(x);
    
    Jac(14,11)=dF14dx6(y);
    
    Jac(14,12)=dF14dy6(x);
    
    Jac(14,13)=dF14dx7(y);
        
    Jac(14,14)=dF14dy7(y);
    
    Jac(15,13)=dF15dx7(y);
    
    Jac(15,14)=dF15dy7(x);
    
    Jac(15,15)=dF15dx8(y);
    
    Jac(15,16)=dF15dy8(x);
    
    Jac(16,13)=dF16dx7(y);
    
    Jac(16,14)=dF16dy7(x);
    
    Jac(16,15)=dF16dx8(y);
        
    Jac(16,16)=dF16dy8(y);
    
    Jac(17,15)=dF17dx8(y);
    
    Jac(17,16)=dF17dy8(x);
    
    Jac(17,17)=dF17dx9(y);
    
    Jac(17,18)=dF17dy9(x);
    
    Jac(18,15)=dF18dx8(y);
    
    Jac(18,16)=dF18dy8(x);
    
    Jac(18,17)=dF18dx9(y);
        
    Jac(18,18)=dF18dy9(y);
    
    Jac(19,19)=dF19dx10(y);
    
    Jac(19,20)=dF19dy10(x);
    
    Jac(20,17)=dF20dx9(y);
    
    Jac(20,18)=dF20dy9(x);
    
    Jac(20,19)=dF20dx10(y);
        
    Jac(20,20)=dF20dy10(y);
    
    F(1,1)=F1(x,y);
    
    F(2,1)=F2(x,y);
    
    F(3,1)=F3(x,y);
    
    F(4,1)=F4(x,y);
    
    F(5,1)=F5(x,y);
    
    F(6,1)=F6(x,y);
    
    F(7,1)=F7(x,y);
    
    F(8,1)=F8(x,y);
    
    F(9,1)=F9(x,y);
    
    F(10,1)=F10(x,y);
    
    F(11,1)=F11(x,y);
    
    F(12,1)=F12(x,y);
    
    F(13,1)=F13(x,y);
    
    F(14,1)=F14(x,y);
    
    F(15,1)=F15(x,y);
    
    F(16,1)=F16(x,y);
    
    F(17,1)=F17(x,y);
    
    F(18,1)=F18(x,y);
    
    F(19,1)=F19(x,y);
    
    F(20,1)=F20(x,y);
    
    Deltas=-F\Jac;
    
    for sep=1:2*Functions
        
        if mod(sep,2)
            
            dx((sep+1)/2)=Deltas(1,sep);
            
        else
            
            dy(sep/2)=Deltas(1,sep);
            
        end
    end
    
    xnew=x+dx;
    
    ynew=y+dy;
    
    errx=abs((xnew-x)./x);
    
    erry=abs((ynew-y)./y);
    
    x=xnew;
    
    y=ynew;
    
    maxerrx=abs(max(errx));
    
    maxerry=abs(max(erry));
    
    if maxerrx<0.0001 
        if maxerry<0.0001
            display('Converged!');
                 
            send_mail_message('izaak.lea','CONVERGED!','CONVERGED','R1.m') 
            
            break
        end
    end
    
    display('Count:', num2str(iter));
   
    
end


send_mail_message('izaak.lea','NOT CONVERGED!','NOT CONVERGED','R1.m')     
     
