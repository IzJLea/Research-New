

xi=0.1;

fun=((1-0.1)^2/(1-x))+(0.1^2/x*55/709)-(0.5*(1-0.1)^2/(1-x)^2);

al=fzero(fun,0.1);


Display(al);

