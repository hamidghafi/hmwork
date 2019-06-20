%% now this file is changed
%% this coment exist in second branch.
clc
f=1e10;u0=4*pi*1e-7;e0=1/36/pi*1e-9;
w=2*pi*f;
Imped=[];
lambda=3e8/f;
k0=w/3e8;
%first the antenna should be discretize in N1 N2 N3 part for every
%parasitic element
%For Reflector
for N1=25:8:65

L1=0.475*lambda;
dl1=L1/(N1+2);
X=zeros(N1,1);
Y=zeros(N1,1);
Z=(linspace(dl1,L1-dl1,N1))';
DL=ones(N1,1)*dl1;
%For Dipole
N2=N1;
L2=0.46*lambda;
dl2=L2/(N2+2);
X=[X;ones(N2,1)*lambda/4];
Y=[Y;zeros(N2,1)];
Z=[Z;(linspace(dl2+.015*lambda/2,L2-dl2+.015*lambda/2,N2))'];
DL=[DL;ones(N2,1)*dl2];
%For Director
N3=N1;
L3=0.44*lambda;
dl3=L3/(N3+2);
X=[X;ones(N3,1)*lambda*0.56];
Y=[Y;zeros(N3,1)];
Z=[Z;(linspace(dl3+.035*lambda/2,L3-dl3+.035*lambda/2,N3))'];
DL=[DL;ones(N3,1)*dl3];
% %------------------------------------------------------------------
N=N1+N2+N3;
d=lambda/1000;
for m=1:N
    for n=1:N
        Reff=Psi_1(Z(n),Z(m),DL(m),DL(n),d+abs(X(n)-X(m)));
        psi(m,n)=DL(m)*DL(n)*1/4/pi/Reff*exp(-i*k0*Reff);

        Reffpp=Psi_1(Z(n)+DL(n),Z(m)+DL(m),DL(m),DL(n),d+abs(X(n)-X(m)));
        Psipp(m,n)=1/4/pi/Reffpp*exp(-i*k0*Reffpp);

        Reffpn=Psi_1(Z(n)+DL(n),Z(m),DL(m),DL(n),d+abs(X(n)-X(m)));
        Psipn(m,n)=1/4/pi/Reffpn*exp(-i*k0*Reffpn);

        Reffnp=Psi_1(Z(n),Z(m)+DL(m),DL(m),DL(n),d+abs(X(n)-X(m)));
        Psinp(m,n)=1/4/pi/Reffnp*exp(-i*k0*Reffnp);

        Reffnn=Psi_1(Z(n),Z(m),DL(m),DL(n),d+abs(X(n)-X(m)));
        Psinn(m,n)=1/4/pi/Reffnn*exp(-i*k0*Reffnn);
        Z_an(m,n)=i*w*u0*psi(m,n)+(Psipp(m,n)+Psinn(m,n)-Psinp(m,n)-Psipn(m,n))/i/w/e0;
    end
end
V=[zeros(1,N1),zeros(1,(N2-1)/2),1,zeros(1,(N2-1)/2),zeros(1,N3)];
% V=[zeros(1,N1),zeros(1,(N2-1)/2),1,zeros(1,(N2-1)/2)];
% V=[zeros(1,((N1-1)/2)),1,zeros(1,((N1-1)/2))];
Y=inv(Z_an);
J=Y*V';
%****************************************************************
J1=[0;J(1:N1);0];
x=linspace(0,L1,length(J1))/lambda;
figure(1)
subplot(3,1,1)
hold on
plot(x,abs(J1)*1000,'b')
J2=[0;J(N1+1:N2+N1);0];
x=linspace(0,L2,length(J2))/lambda;
title(['Current distribution for Yagi antenna d=\lambda/1000 L = .475\lambda Reflector'],'FontSize',13)
xlabel('l/\lambda','FontSize',12)
ylabel('Milliamperes','FontSize',12)
subplot(3,1,2)
hold on
plot(x,abs(J2)*1000,'b-')
title(['Current distribution for Yagi antenna d=\lambda/1000 L = .46\lambda Dipole'],'FontSize',13)
xlabel('l/\lambda','FontSize',12)
ylabel('Milliamperes','FontSize',12)
subplot(3,1,3)
hold on
J3=[0;J(N1+N2+1:N);0];
x=linspace(0,L3,length(J3))/lambda;
plot(x,abs(J3)*1000,'c-.')
title(['Current distribution for Yagi antenna d=\lambda/1000 L = .44\lambda Director'],'FontSize',13)
xlabel('l/\lambda','FontSize',12)
ylabel('Milliamperes','FontSize',12)
%****************************************************************
g=@(Theta,Phi)sin(Theta)*sum(J.*exp(i*k0*(cos(Theta)*Z+sin(Theta)*cos(Phi)*X)));
Patternt1=[];
theta=-pi:.005:pi;
phi=linspace(0,2*pi,100);
theta1=linspace(0,pi,100);
P=0;
for th=theta
    Patternt1=[Patternt1,20*log10(abs(g(th,0)))];
end
for th=theta1
    for ph=phi
        P=abs(g(th,ph))^2*(2*pi/100)*(pi/100)*sin(th)+P;
    end
end
figure(2)
        G=10*log10(4*pi/P)+max(Patternt1);
        disp(['Max Gain for Yagi antenna : ',num2str(G)])
        Temp=find(Patternt1<(max(Patternt1)-15));
        Patternt1(Temp)=[];
        theta(Temp)=[];
        Polar(theta,Patternt1-min(Patternt1),[-10,-5,0],[5,10,15])
        title('Raidiation Pattern for 3 Element Yagi Antenna','FontSize',13)
        hold on
        disp(['Input Impedace of Yagi antenna L = ',num2str(1/J(N1+(N2-1)/2))])        

end