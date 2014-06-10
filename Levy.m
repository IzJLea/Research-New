function F = Levy(a)

a=a(1);

x=a(2);

Dens=a(3);

F = (((1-x)^2)/(1-a))+(x^2*Dens/a)-(0.5*((1-x)^2)*((1-a)^2))-0.5;

display(a, 'void fraction');

