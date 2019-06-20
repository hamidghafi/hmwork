clc
clear
f=10^10;
C=3*10^8;
omega=2*pi*f;
mue0=4*pi*10^-7;
eps0=8.85*10^-12;
lambda=C/f;
K=2*pi/lambda;
N1=135;
N2=N1;
N3=N1;
L1=.475*lambda;
L2=.46*lambda;
L3=.44*lambda;
l1=L1/(N1+1);
l2=L2/(N2+1);
l3=L3/(N3+1);
Z11=[];
Z12=[];
Z13=[];
Z14=[];
Z22=[];
Z23=[];
Z24=[];
Z33=[];
Z34=[];
Z44=[];
a=.002*lambda;
p=0;
sai=[];
sai1=[];
sai2=[];
sai3=[];
I=[];

 V=zeros(N1+N2+N3,1);
 V(N1+(N2+1)/2,1)=1;
 sai=fai3(L1,N1);
for m=1:N1
        for n=1:N1
           
            if (n==m)
            Z11(m,n)=omega*mue0*complex(0,1)*l1^2*sai(1)+1/(omega*eps0*complex(0,1))*(2*sai(1)-2*sai(2));
            else
            Z11(m,n)=omega*mue0*complex(0,1)*l1^2*sai(abs(m-n)+1)+1/(omega*eps0*complex(0,1))*(2*sai(abs(m-n)+1)-sai(abs(m-n))-sai(abs(m-n)+2));
                
            end
        end
end
%
sai=fai3(L2,N2);
for m=1:N2
        for n=1:N2
            
            if (n==m)
            Z22(m,n)=omega*mue0*complex(0,1)*l2^2*sai(1)+1/(omega*eps0*complex(0,1))*(2*sai(1)-2*sai(2));
            else
                Z22(m,n)=omega*mue0*complex(0,1)*l2^2*sai(abs(m-n)+1)+1/(omega*eps0*complex(0,1))*(2*sai(abs(m-n)+1)-sai(abs(m-n))-sai(abs(m-n)+2));
                
            end
        end
end
%
sai=fai3(L3,N3);
for m=1:N3
        for n=1:N3

            if (n==m)
            Z33(m,n)=omega*mue0*complex(0,1)*l3^2*sai(1)+1/(omega*eps0*complex(0,1))*(2*sai(1)-2*sai(2));
            else
            Z33(m,n)=omega*mue0*complex(0,1)*l3^2*sai(abs(m-n)+1)+1/(omega*eps0*complex(0,1))*(2*sai(abs(m-n)+1)-sai(abs(m-n))-sai(abs(m-n)+2));
                
            end
        end
end
%



for m=1:N1
        for n=1:N2
            dr=.0075*lambda;
            d=.25*lambda;
            delta1=(n-1/2)*l2+dr-(m+1/2)*l1;
            R1=Reff(l1,l2,delta1,l1+l2+delta1,l1+delta1,l2+delta1,d);
            phai1=exp(-j*K*R1)/4/pi/R1;
            %m+,n+
            delta2=(n)*l2+dr-(m+1)*l1;
            R2=Reff(l1,l2,delta2,l1+l2+delta2,l1+delta2,l2+delta2,d);
            phai2=exp(-j*K*R2)/4/pi/R2;
            %m+,n-
            delta3=(n-1)*l2+dr-(m+1)*l1;
            R3=Reff(l1,l2,delta3,l1+l2+delta3,l1+delta3,l2+delta3,d);
            phai3=exp(-j*K*R3)/4/pi/R3;
            %m-,n+
            
            delta4=(n)*l2+dr-(m)*l1;
             R4=Reff(l1,l2,delta4,l1+l2+delta4,l1+delta4,l2+delta4,d);
            phai4=exp(-j*K*R4)/4/pi/R4;
            %m-,n-
            delta5=(n-1)*l2+dr-(m)*l1;
            R5=Reff(l1,l2,delta5,l1+l2+delta5,l1+delta5,l2+delta5,d);
            phai5=exp(-j*K*R5)/4/pi/R5;
            %
             Z12(m,n)=omega*mue0*l1*l2*complex(0,1)*phai1+1/(omega*eps0*complex(0,1))*(phai2-phai3-phai4+phai5);
            
            
            
           
           
        end
end

for m=1:N1
        for n=1:N3
              dr=.0175*lambda;
             d=.56*lambda;
            delta1=(n-1/2)*l2+dr-(m+1/2)*l1;
            R1=Reff(l1,l2,delta1,l1+l2+delta1,l1+delta1,l2+delta1,d);
            phai1=exp(-j*K*R1)/4/pi/R1;
            %m+,n+
            delta2=(n)*l2+dr-(m+1)*l1;
            R2=Reff(l1,l2,delta2,l1+l2+delta2,l1+delta2,l2+delta2,d);
            phai2=exp(-j*K*R2)/4/pi/R2;
            %m+,n-
            delta3=(n-1)*l2+dr-(m+1)*l1;
            R3=Reff(l1,l2,delta3,l1+l2+delta3,l1+delta3,l2+delta3,d);
            phai3=exp(-j*K*R3)/4/pi/R3;
            %m-,n+
            
            delta4=(n)*l2+dr-(m)*l1;
             R4=Reff(l1,l2,delta4,l1+l2+delta4,l1+delta4,l2+delta4,d);
            phai4=exp(-j*K*R4)/4/pi/R4;
            %m-,n-
            delta5=(n-1)*l2+dr-(m)*l1;
            R5=Reff(l1,l2,delta5,l1+l2+delta5,l1+delta5,l2+delta5,d);
            phai5=exp(-j*K*R5)/4/pi/R5;
            %
             Z13(m,n)=omega*mue0*l1*l2*complex(0,1)*phai1+1/(omega*eps0*complex(0,1))*(phai2-phai3-phai4+phai5);
            
        end
end
for m=1:N2
        for n=1:N3
           
             dr=.01*lambda;
            d=.31*lambda;
            delta1=(n-1/2)*l2+dr-(m+1/2)*l1;
            R1=Reff(l1,l2,delta1,l1+l2+delta1,l1+delta1,l2+delta1,d);
            phai1=exp(-j*K*R1)/4/pi/R1;
            %m+,n+
            delta2=(n)*l2+dr-(m+1)*l1;
            R2=Reff(l1,l2,delta2,l1+l2+delta2,l1+delta2,l2+delta2,d);
            phai2=exp(-j*K*R2)/4/pi/R2;
            %m+,n-
            delta3=(n-1)*l2+dr-(m+1)*l1;
            R3=Reff(l1,l2,delta3,l1+l2+delta3,l1+delta3,l2+delta3,d);
            phai3=exp(-j*K*R3)/4/pi/R3;
            %m-,n+
            
            delta4=(n)*l2+dr-(m)*l1;
             R4=Reff(l1,l2,delta4,l1+l2+delta4,l1+delta4,l2+delta4,d);
            phai4=exp(-j*K*R4)/4/pi/R4;
            %m-,n-
            delta5=(n-1)*l2+dr-(m)*l1;
            R5=Reff(l1,l2,delta5,l1+l2+delta5,l1+delta5,l2+delta5,d);
            phai5=exp(-j*K*R5)/4/pi/R5;
            Z23(m,n)=omega*mue0*l1*l2*complex(0,1)*phai1+1/(omega*eps0*complex(0,1))*(phai2-phai3-phai4+phai5);
            
        end
end
%

Z=[Z11,Z12,Z13;Z12,Z22,Z23;Z13,Z23,Z33];
I=inv(Z)*V;
I1=I(1:N1);
I2=I(N1+1:N1+N2);
I3=I(N1+N2+1:N1+N2+N3);     
Zin=1/I(N1+(N2+1)/2);
X1=[-L1/2:l1:L1/2];
X2=[-L2/2:l2:L2/2];
X3=[-L3/2:l3:L3/2];
figure(1)
subplot(3,1,1)
plot(X1,[0;abs(I1)*1000;0])
subplot(3,1,2)
plot(X2,[0;abs(I2)*1000;0])
subplot(3,1,3)
plot(X3,[0;abs(I3)*1000;0])
teta=[0:.1:pi];
g=zeros(1,length(teta));
for m=1:N1+N2+N3
        for q=1:length(g)
        g(q)=g(q)+(I(m)*exp(j*K*(m*l1-L1/2)*cos(teta(q))))*sin(abs(teta(q)));
        end    
end
g0=abs(g);
g1=(max(g0));
g2=g./g1;
figure(2)
polar(teta,g2)
title(['pattern for phai=0 and L=',num2str(L2)])

