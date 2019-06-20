function c=fai(L,N)
% this function is to
f=10^7;
C=3*10^8;
lambda=C/f;
K=2*pi/lambda;
l=L/(N+1);
a=.001*lambda;
Reff=0;
for m=1:N+1
        d=a;
        delta=(m-2)*l;
        alpha=2*l+delta;
        beta=l+delta;
        gama=l+delta;
        Reff=l^2/(alpha*asinh(alpha/d)-beta*asinh(beta/d)-gama*asinh(gama/d)....
        +delta*asinh(delta/d)-sqrt(alpha^2+d^2)+sqrt(beta^2+d^2)+sqrt(gama^2+d^2)-sqrt(delta^2+d^2));
        c(m)=exp(complex(0,-K*Reff))/4/pi/Reff;
    end
end
    