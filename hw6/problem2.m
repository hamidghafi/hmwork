clc
clear
f=10^7;
C=3*10^8;
omega=2*pi*f;
mue0=4*pi*10^-7;
eps0=8.85*10^-12;
lambda=C/f;
K=2*pi/lambda;
N=123;
L=lambda*[1/4 1/2 1 3/2];
Z=[];
a=.001*lambda;
Reff=0;
p=0;
sai=[];
I=[];
V=zeros(N,1);
V(N/3)=1;
Yin=[];
    Gin=[];
    Bin=[];
 teta=[-pi:.01:pi];
     g=zeros(1,length(teta));
     h=zeros(1,length(teta));
for i=1:4
    p=p+1;
    l=L(i)/(N+1);
    X=[-L(i)/2:l:L(i)/2];
    
    sai=fai(L(i),N);
    for m=1:N
        for n=1:N
            if (n==m)
            Z(m,n)=omega*mue0*complex(0,1)*l^2*sai(1)+1/(omega*eps0*complex(0,1))*(2*sai(1)-2*sai(2));
            else
            Z(m,n)=omega*mue0*complex(0,1)*l^2*sai(abs(m-n)+1)+1/(omega*eps0*complex(0,1))*(2*sai(abs(m-n)+1)-sai(abs(m-n))-sai(abs(m-n)+2));
            end
        end
    end
    I(:,i)=inv(Z)*V;
    I1=[0;abs(I(:,i)*1000);0];
    figure (i)
    plot (X,I1')
    xlabel ('z axis')
    ylabel ('current distebution')
    title(['current disterbution for L=',num2str(L(i))])
     for m=1:N
        for q=1:length(g)
        g(q)=g(q)+(I(m,i)*exp(j*K*(m*l-L(i)/2)*cos(teta(q))))*sin(abs(teta(q)));
        end
        for o=1:length(g)
        h(o)=h(o)+(I(m,i)*exp(j*K*(m*l-L(i)/2)*cos(pi/2)))*sin(abs(pi/2));
        end
        h1=max(h);
        h2=h./h1;
        g1=(max(g));
        g2=g./g1;
    end
    figure(i+40)
    subplot(1,2,1)
    polar(teta,g2)
    title(['pattern for phai=0 and L=',num2str(L(i))])
    subplot(1,2,2)
    polar(teta,h2)
    title(['pattern for teta=pi/2 and L=',num2str(L(i))])
end
Yin=I((N+1)/2,:);
Zin1=1./Yin;
Gin=real(Zin1);
Bin=imag(Zin1);
figure (40)
X1=L./lambda;
plot(X1,Gin,'--r')
title ('          real part of impedenace by -- graph and imaginary by continuse graph')
xlabel ('L/lambda')
ylabel ('real and imaginary part of  impedance')
hold on
plot(X1,Bin,'b')
hold off
