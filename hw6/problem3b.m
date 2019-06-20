clc
clear
f=10^10;
C=3*10^8;
omega=2*pi*f;
mue0=4*pi*10^-7;
eps0=8.85*10^-12;
lambda=C/f;
K=2*pi/lambda;
N1=123;
N2=N1;
N3=N1;
N4=N1;
L1=.475*lambda;
L2=.46*lambda;
L3=.44*lambda;
L4=.56*lambda;
l1=L1/(N1+1);
l2=L2/(N2+1);
l3=L3/(N3+1);
l4=L4/(N4+1);
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
I=[];
 
 V=zeros(N1+N2+N3+N4,1);
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
sai=fai3(L4,N4);
for m=1:N4
        for n=1:N4
            
            if (n==m)
            Z44(m,n)=omega*mue0*complex(0,1)*l4^2*sai(1)+1/(omega*eps0*complex(0,1))*(2*sai(1)-2*sai(2));
            else
                Z44(m,n)=omega*mue0*complex(0,1)*l4^2*sai(abs(m-n)+1)+1/(omega*eps0*complex(0,1))*(2*sai(abs(m-n)+1)-sai(abs(m-n))-sai(abs(m-n)+2));
                
            end
        end
end


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
             Z13(m,n)=omega*mue0*l1*l3*complex(0,1)*phai1+1/(omega*eps0*complex(0,1))*(phai2-phai3-phai4+phai5);
            
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
            Z23(m,n)=omega*mue0*l2*l3*complex(0,1)*phai1+1/(omega*eps0*complex(0,1))*(phai2-phai3-phai4+phai5);
            
        end
end
p=zeros(N1,N1);
for m=1:N1
        for n=1:N4
            if (abs(m-3)<(N1+1)/2& abs(n)<3)
                g=@(p,t)(1./sqrt(p.^2+t.^2));
                M1=1/l1/l4*dblquad(g,L1/2-(m-1/2)*l1,L1/2-(m+1/2)*l1,(n-1/2)*l4,(n+1/2)*l4);
                M2=1/l1/l4*dblquad(g,L1/2-(m)*l1,L1/2-(m+1)*l1,(n)*l4,(n+1)*l4);
                M3=1/l1/l4*dblquad(g,L1/2-(m)*l2,L1/2-(m+1)*l2,(n-1)*l4,(n)*l4);
                M4=1/l1/l4*dblquad(g,L1/2-(m-1)*l1,L1/2-(m)*l2,(n)*l4,(n+1)*l4);
                M5=1/l1/l4*dblquad(g,L1/2-(m-1)*l1,L1/2-(m)*l2,(n-1)*l4,(n)*l4);
                phai1=exp(-j*K/M1)/4/pi*M1;
                phai2=exp(-j*K/M2)/4/pi*M2;
                phai3=exp(-j*K/M3)/4/pi*M3;
                phai4=exp(-j*K/M4)/4/pi*M4;
                phai5=exp(-j*K/M5)/4/pi*M5;
               
            else
               
            z1=(L1/2-m*l1);
            x1=n*l4;
            R1=sqrt(x1^2+z1^2);
            phai1=exp(-j*K*R1)/4/pi/R1;
            %m+,n+
            z2=(L1/2-(m+1/2)*l1);
            x2=(n+1/2)*l4;
            R2=sqrt(x2^2+z2^2);
            phai2=exp(-j*K*R2)/4/pi/R2;
            %m+,n-
            z3=(L1/2-(m+1/2)*l1);
            x3=(n-1/2)*l4;
            R3=sqrt(x3^2+z3^2);
            phai3=exp(-j*K*R3)/4/pi/R3;
            %m-,n+
            
            z4=(L1/2-(m-1/2)*l1);
            x4=(n+1/2)*l4;
            R4=sqrt(x4^2+z4^2);
            phai4=exp(-j*K*R4)/4/pi/R4;
            %m-,n-
            z5=(L1/2-(m-1/2)*l1);
            x5=(n-1/2)*l4;
            R5=sqrt(x5^2+z5^2);
            phai5=exp(-j*K*R5)/4/pi/R5;
            end
            Z14(m,n)=omega*mue0*l1*l4*complex(0,1)*phai1+1/(omega*eps0*complex(0,1))*(phai2-phai3-phai4+phai5);
            
        end
end

%

for m=1:N2
        for n=1:N4
           if (abs(m-3)<(N2+1)/2& abs(n-3)<(N4+1)/2)
                g=@(p,t)(1./sqrt(p.^2+t.^2));
                M1=1/l2/l4*dblquad(g,L2/2-(m-1/2)*l2,L2/2-(m+1/2)*l2,L4/2-(n-1/2)*l4,L4/2-(n+1/2)*l4);
                M2=1/l2/l4*dblquad(g,L2/2-(m)*l2,L2/2-(m+1)*l2,L4/2-(n)*l4,L4/2-(n+1)*l4);
                M3=1/l2/l4*dblquad(g,L2/2-(m)*l2,L2/2-(m+1)*l2,L4/2-(n-1)*l4,L4/2-(n)*l4);
                M4=1/l2/l4*dblquad(g,L2/2-(m-1)*l2,L2/2-(m)*l2,L4/2-(n)*l4,L4/2-(n+1)*l4);
                M5=1/l2/l4*dblquad(g,L2/2-(m-1)*l2,L2/2-(m)*l2,L4/2-(n-1)*l4,L4/2-(n)*l4);
                phai1=exp(-j*K/M1)/4/pi*M1;
                phai2=exp(-j*K/M2)/4/pi*M2;
                phai3=exp(-j*K/M3)/4/pi*M3;
                phai4=exp(-j*K/M4)/4/pi*M4;
                phai5=exp(-j*K/M5)/4/pi*M5;
               
            else
            
            
            z1=(L2/2-m*l2);
            x1=L4/2-n*l4;
            R1=sqrt(x1^2+z1^2);
            phai1=exp(-j*K*R1)/4/pi/R1;
            save(m,n)=phai1;
           
           
            %m+,n+
            z2=(L2/2-(m+1/2)*l2);
            x2=L4/2-(n+1/2)*l4;
            R2=sqrt(x2^2+z2^2);
            phai2=exp(-j*K*R2)/4/pi/R2;
            %m+,n-
            z3=(L2/2-(m+1/2)*l2);
            x3=L4/2-(n-1/2)*l4;
            R3=sqrt(x3^2+z3^2);
            phai3=exp(-j*K*R3)/4/pi/R3;
            %m-,n+
            
            z4=(L2/2-(m-1/2)*l2);
            x4=L4/2-(n+1/2)*l4;
            R4=sqrt(x4^2+z4^2);
            phai4=exp(-j*K*R4)/4/pi/R4;
            %m-,n-
            z5=(L2/2-(m-1/2)*l2);
            x5=L4/2-(n-1/2)*l4;
            R5=sqrt(x5^2+z5^2);
            phai5=exp(-j*K*R5)/4/pi/R5;
           end
            Z24(m,n)=omega*mue0*l2*l4*complex(0,1)*phai1+1/(omega*eps0*complex(0,1))*(phai2-phai3-phai4+phai5);
            
        end
end
%

for m=1:N3
        for n=1:N4
            if (abs(m-3)<(N3+1)/2& abs(n-N4)<(3))
                g=@(p,t)(1./sqrt(p.^2+t.^2));
                M1=1/l3/l4*dblquad(g,L3/2-(m-1/2)*l3,L3/2-(m+1/2)*l3,L4-(n-1/2)*l4,L4-(n+1/2)*l4);
                M2=1/l3/l4*dblquad(g,L3/2-(m)*l3,L3/2-(m+1)*l3,L4-(n)*l4,L4-(n+1)*l4);
                M3=1/l3/l4*dblquad(g,L3/2-(m)*l3,L3/2-(m+1)*l3,L4-(n-1)*l4,L4-(n)*l4);
                M4=1/l3/l4*dblquad(g,L3/2-(m-1)*l3,L3/2-(m)*l3,L4-(n)*l4,L4-(n+1)*l4);
                M5=1/l3/l4*dblquad(g,L3/2-(m-1)*l3,L3/2-(m)*l3,L4-(n-1)*l4,L4-(n)*l4);
                phai1=exp(-j*K/M1)/4/pi*M1;
                phai2=exp(-j*K/M2)/4/pi*M2;
                phai3=exp(-j*K/M3)/4/pi*M3;
                phai4=exp(-j*K/M4)/4/pi*M4;
                phai5=exp(-j*K/M5)/4/pi*M5;
               
            else
            z1=(L3/2-m*l3);
            x1=L4-n*l4;
            R1=sqrt(x1^2+z1^2);
            phai1=exp(-j*K*R1)/4/pi/R1;
            %m+,n+
            z2=(L3/2-(m+1/2)*l3);
            x2=L4-(n+1/2)*l4;
            R2=sqrt(x2^2+z2^2);
            phai2=exp(-j*K*R2)/4/pi/R2;
            %m+,n-
            z3=(L3/2-(m+1/2)*l3);
            x3=L4-(n-1/2)*l4;
            R3=sqrt(x3^2+z3^2);
            phai3=exp(-j*K*R3)/4/pi/R3;
            %m-,n+
            
            z4=(L3/2-(m-1/2)*l3);
            x4=L4-(n+1/2)*l4;
            R4=sqrt(x4^2+z4^2);
            phai4=exp(-j*K*R4)/4/pi/R4;
            %m-,n-
            z5=(L3/2-(m-1/2)*l3);
            x5=L4-(n-1/2)*l4;
            R5=sqrt(x5^2+z5^2);
            phai5=exp(-j*K*R5)/4/pi/R5;
            end
            Z34(m,n)=omega*mue0*l3*l4*complex(0,1)*phai1+1/(omega*eps0*complex(0,1))*(phai2-phai3-phai4+phai5);
            
        end
end

Z=[Z11,Z12,Z13,Z14;Z12,Z22,Z23,Z24;Z13,Z23,Z33,Z34;Z14,Z24,Z34,Z44];
I=inv(Z)*V;
I1=I(1:N1);
I2=I(N1+1:N1+N2);
I3=I(N1+N2+1:N1+N2+N3);
I4=I(N1+N2+N3+1:N1+N2+N3+N4);     
Zin=1/I(N1+(N2+1)/2);
X1=[-L1/2:l1:L1/2];
X2=[-L2/2:l2:L2/2];
X3=[-L3/2:l3:L3/2];
X4=[-L4/2:l4:L4/2];
figure(1)
subplot(4,1,1)
plot(X1,[0;abs(I1)*1000;0])
subplot(4,1,2)
plot(X2,[0;abs(I2)*1000;0])
subplot(4,1,3)
plot(X3,[0;abs(I3)*1000;0])
subplot(4,1,4)
plot(X4,[0;abs(I4)*1000;0])
teta=[0:.1:2*pi];
g=zeros(1,length(teta));
for m=1:N2+N1+N3
        for q=1:length(g)
        g(q)=g(q)+(I(m)*exp(j*K*(m*l1-L1/2)*cos(teta(q))))*sin(abs(teta(q)));
        end    
end
g0=(abs(g));
g1=max(g0);
g2=g./g1;
figure(2)
polar(teta,g2)
title(['pattern for phai=0 and L=',num2str(L2)])
    
