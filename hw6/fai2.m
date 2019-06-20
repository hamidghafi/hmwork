function c=fai2(L1,L2);
f=10^8;
C=3*10^8;
lambda=C/f;
K=2*pi/lambda;
N1=99;
N2=N1;
d=.01*lambda;
l1=L1/(N1+1);
l2=L2/(N2+1);
Reff=0;
for m=1:N1+1
    for n=1:N2+1
        if (abs(m-3)<((N1+1)/2))
            g=@(x,z)(1./(sqrt(x.^2+z.^2+d^2)))
            Reff=l1*l2/(dblquad(g,(m-N1/2-1/2)*l1,(m+1-N1/2-1/2)*l1,(n-1)*l2,n*l2))
            c(m,n)=exp(complex(0,-K*Reff))/4/pi/Reff;
        else
            Reff=sqrt(((m-N1/2)*l1)^2+((n-1/2)*l2)^2+d^2);
            c(m,n)=exp(complex(0,-K*Reff))/4/pi/Reff;
        end
    end
end