function res=Alphacalc(Qchannel,mflow,Tenter,Peval)

Bundpower=cosPower(12,Qchannel); %cosine distrubution of power in core

mmax=Bundpower/1000/(XSteam('hV_p',Peval)-XSteam('hL_p',Peval)); %maximum generated mass flow per bundle

%enthalpies 

hout=zeros(1,12); 

hin=zeros(1,12);

hin(1,1)=XSteam('h_pT',Peval,Tenter)*1000;

%determining the initial point of vaporization

for p=1:12
    hout(1,p)=hin(1,p)+Bundpower(1,p);
    
    if hout(1,p)>XSteam('hL_p',Peval)
        
        start=p;
        
        break
    else
        
        hin(1,p+1)=hout(1,p);
    end
end

% number of equations
r=(12-(start-1))*2;

%Initiation of functions

for p=2:2:r
    if p~=r
        eval(['F' num2str(p-1) '= @(x1,y1) 2*x1*y1-x1-y1']);
    else
        
        eval(['F' num2str(p-1) '= @(x' num2str(p/2) ',y' num2str(p/2) ') 1-x' num2str(p/2) '*y' num2str(p/2)]);
        
    end
    
    if r==2
        eval(['F' num2str(p-1) '= @(x1,y1) (1-x1)*mmax(1,p/2)-x1*y1*mflow']);
        
       
    else
        
        eval(['F' num2str(p) '= @(x' num2str(p/2-1) ',y' num2str(p/2-1) ',x' num2str(p/2) ',y' num2str(p/2) ') x' num2str(p/2-1) '*y' num2str(p/2-1) '*mflow-x' num2str(p/2) '*y' num2str(p/2) '*mflow+(1-x' num2str(p/2) ')*mmax(1,p/2)']);
        
    end
end


%derivative equation generation


df1dx1=@(y1) -2*y1+1;

df1dy1=@(x1) -2*x1+1;

df2dx1=@(y1) -y1*mflow-mmax(1,1);

df2dy1=@(x1) -x1*mflow;

for p=4:2:r
    
    
    if p==r
        
        eval(['df' num2str(p) 'dx' num2str(p/2-1) '=@(y' num2str(p/2-1) ') y' num2str(p/2-1) '*mflow']);
        eval(['df' num2str(p) 'dy' num2str(p/2-1) '=@(x' num2str(p/2-1) ') x' num2str(p/2-1) '*mflow']);
        eval(['df' num2str(p-1) 'dx' num2str(p/2) '=@(y' num2str(p/2) ') -y' num2str(p/2)]);
        eval(['df' num2str(p-1) 'dy' num2str(p/2) '=@(x' num2str(p/2) ') -x' num2str(p/2)]);
        eval(['df' num2str(p) 'dx' num2str(p/2) '=@(y' num2str(p/2) ') -y' num2str(p/2) '*mflow-mmax(1,' num2str(p/2) ')']);
        eval(['df' num2str(p) 'dy' num2str(p/2) '=@(x' num2str(p/2) ') -x' num2str(p/2) '*mflow']);
    else
            
            eval(['df' num2str(p-1) 'dx' num2str(p/2-1) '=@(y' num2str(p/2-1) ')= 2*y' num2str(p/2-1) '-1']);
            eval(['df' num2str(p) 'dx' num2str(p/2-1) '=@(y' num2str(p/2-1) ')= y' num2str(p/2-1) '*mflow']);
   
            eval(['df' num2str(p-1) 'dy' num2str(p/2-1) '=@(x' num2str(p/2-1) ')= 2*x' num2str(p/2-1) '-1']);
            eval(['df' num2str(p) 'dy' num2str(p/2-1) '=@(x' num2str(p/2-1) ')= x' num2str(p/2-1) '*mflow']);

            eval(['df' num2str(p-1) 'dx' num2str(p/2) '=@(y' num2str(p/2) ')= -2*y' num2str(p/2) '+1']);
            eval(['df' num2str(p) 'dx' num2str(p/2) '=@(y' num2str(p/2) ')= -y' num2str(p/2) '*mflow-mmax(1,' num2str(p/2-start+1) ')']);

            eval(['df' num2str(p-1) 'dy' num2str(p/2) '=@(x' num2str(p/2) ')= -2*x' num2str(p/2) '+1']);
            eval(['df' num2str(p) 'dy' num2str(p/2) '=@(x' num2str(p/2) ')=-x' num2str(p/2) '*mflow']);
        
    end


end

% % Rewrite to make single matrices for x and y
% for l=1:r/2
%     
%     eval(['x' num2str(1) '=0.5']);
%     
%     eval(['y' num2str(1) '=0.5']);
%     
% end

Jac=zeros(r);
F=zeros(r,1);
while errx>0.00001 && erry>0.00001
    
    Jac(1,1)=df1dx1(y1);
    Jac(1,2)=df1dy1(x1);
    Jac(2,1)=df2dx1(y1);
    Jac(2,2)=df2dy1(x1);
    
    for t=3:r
        if mod(t,2)
            
            eval(['Jac(t,t-2)=df' num2str(t) 'dx' num2str((t+1)/2-1) '(y' num2str((t+1)/2-1) ')' ]);
            
            eval(['Jac(t,t-1)=df' num2str(t) 'dy' num2str((t+1)/2-1) '(x' num2str((t+1)/2-1) ')' ]);
            
            eval(['Jac(t,t)=df' num2str(t) 'dx' num2str((t+1)/2) '(y' num2str((t+1)/2) ')' ]);
            
            eval(['Jac(t,t+1)=df' num2str(t) 'dy' num2str((t+1)/2) '(x' num2str((t+1)/2) ')' ]);
        else
            
            eval(['Jac(t,t-3)=df' num2str(t) 'dx' num2str(t/2-1) '(y' num2str(t/2-1) ')' ]);
            
            eval(['Jac(t,t-2)=df' num2str(t) 'dy' num2str(t/2-1) '(x' num2str(t/2-1) ')' ]);
            
            eval(['Jac(t,t-1)=df' num2str(t) 'dx' num2str(t/2) '(y' num2str(t/2) ')' ]);
            
            eval(['Jac(t,t)=df' num2str(t) 'dx' num2str(t/2) '(y' num2str(t/2) ')' ]);
            
        end
    end
    
    
    F(1,1)=F1(x1,y1);
    
    F(2,1)=F2(x1,y1);
    
    for m=3:r
        
        if mod(m,2)
        
            eval(['F(m,1)=F' num2str(m) '(x' num2str((t+1)/2-1) ', y' num2str((t+1)/2-1) ',x' num2str((t+1)/2) ',y' num2str((t+1)/2) ')']);
        else
            eval(['F(m,1)=F' num2str(m) '(x' num2str(t/2-1) ', y' num2str(t/2-1) '(x' num2str(t/2) ', y' num2str(t/2) ')']);
        end
    end
    
    Delta=-F\Jac;
    
    for k=1:r
        
        if mod(k,2) 
            
            eval(['dx' num2str((k+1)/2) '=Delta(' num2str(k) ',1)']);
        else
            eval(['dy' num2str(k/2) '=Delta(' num2str(k) ',1)']);
        end
    end
    
    eval([X

