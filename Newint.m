function Newt = Newint(x,y,xx)

n=length(x);
%% set up b matrix
b(:,1)=y(:);

if n~=length(y)
    error('x and y must be equal')
end
for j=2:n
    for i= 1:n-j+1
        b(i,j)=((b(i+1,j-1)-(b(i,j-1))))/(x(i+j-1)-x(i));
    end
end

Newt=b(1,1);

xt=1;

for k=1:n-1
    
    xt=xt*(xx-x(k));
    
    Newt=Newt+b(1,k+1)*xt;
    
end
