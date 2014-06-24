function Look= RevLinLook(x,y,xx)

n=length(x);

i=1;
while (1)
    if xx>=x(i+1)
        break
    end
    i=i+1;
end

Look=y(i)+(y(i)-y(i+1))/(x(i)-x(i+1))*(x(i)-xx);