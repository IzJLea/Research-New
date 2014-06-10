%%  Basic Lagrange Interpolation
%for use with steam tables.
% general function for any data set

%% Input

%This will interpolate for the saturation pressure 

function Lag= Lagint(x,y,xx)

xl=length(x);
%check that length of both are the same

if xl~=length(y)
    error('length of x and y must be equal');
end

%set up summation

s=0;

for i=1:xl
    
    funx=y(i);
    
    for p=1:xl
        
        if p~=i
            
            funx=(xx-x(p))/(x(i)-x(p))*funx;
        end
        
    end
    s=s+funx;
end

Lag = s;



