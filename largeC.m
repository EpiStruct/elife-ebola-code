function C = largeC(M1,M2,p,q)
% See numerics.pdf for an explanation of this file

latilde = ((1-p)/(1-q));
Y = zeros(1,M2);
for x=M1:M2
    yy = (((4*q*(1-q))^(x-1))/(sqrt(pi*(x-1))));
    Y(x) = yy;
    for a=2:x
        yy = (latilde*(x-a)/((2*x)-a-1))*yy;
        Y(x) = Y(x)+(a*yy);
    end
    Y(x) = Y(x)*(p*q/x);
end
C = sum(Y);

end
