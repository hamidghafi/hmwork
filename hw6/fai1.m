function c=fai1(L1,L2,d)
f=10^8;
C=3*10^8;
lambda=C/f;
K=2*pi/lambda;
N1=99;
N2=N1;
l1=L1/N1;
l2=L2/N2;
for m=1:N1+1
    for n=1:N2+1
        delta=abs(((n-1)*l2-m*l1));
        alpha=l1+l2+delta;
        beta=l1+delta;
        gama=l2+delta;
        Reff=l1*l2/(alpha*asinh(alpha/d)-beta*asinh(beta/d)-gama*asinh(gama/d)....
        +delta*asinh(delta/d)-sqrt(alpha^2+d^2)+sqrt(beta^2+d^2)+sqrt(gama^2+d^2)-sqrt(delta^2+d^2));
        c(m,n)=exp(complex(0,-K*Reff))/4/pi/Reff;
    end
end
