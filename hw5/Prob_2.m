function Prob_2()
clc
close all
f=1e7;u0=4*pi*1e-7;e0=1/36/pi*1e-9;
w=2*pi*f;
Imped=[];
lambda=3e8/f;
k0=w/3e8;
F=0;
c=[.25,.5,1,1.5];
Fp=1/3;
for C=c
    F=F+1;
    for N=45:6:57
        L=lambda*C;
        d=lambda/1000;
        dl=L/(N+2);
        Z=[];
        psi=[];
        Psipp=[];Psinn=[];Psinp=[];Psipn=[];
        Z=[];
        for m=1:N
            n=1;
            Reff=Psi(m*dl,dl,dl,d);
            psi(m,n)=dl^2*1/4/pi/Reff*exp(-i*k0*Reff);
            Reffpp=Psi(m*dl+dl,2*dl,dl,d);
            Psipp(m,n)=1/4/pi/Reffpp*exp(-i*k0*Reffpp);

            Reffpn=Psi(m*dl+dl,dl,dl,d);
            Psipn(m,n)=1/4/pi/Reffpn*exp(-i*k0*Reffpn);

            Reffnp=Psi(m*dl,2*dl,dl,d);
            Psinp(m,n)=1/4/pi/Reffnp*exp(-i*k0*Reffnp);

            Reffnn=Psi(m*dl,dl,dl,d);
            Psinn(m,n)=1/4/pi/Reffnn*exp(-i*k0*Reffnn);
            Z(m,n)=i*w*u0*psi(m,n)+(Psipp(m,n)+Psinn(m,n)-Psinp(m,n)-Psipn(m,n))/i/w/e0;

        end
        V=[];
        J=[];
        for m=1:N
            for n=1:N
                Z(m,n)=Z(abs(m-n)+1,1);
            end
        end
        V=[zeros(1,floor((N-1)*Fp)),1,zeros(1,floor(N*(1-Fp)))];
%         V=[zeros(1,floor((N-1)/2)),1,zeros(1,floor((N-1)/2))];
        Y=inv(Z);
        J=Y*V';
        J=[0;J;0];
        figure(F)
        x=linspace(0,L,length(J));
        plot(x,abs(J)*1000)
        title(['Current distribution for Linear antenna d=lambda/1000 L = ',num2str(C),' \lambda'])
        xlabel('l')
        ylabel('Milliamperes')
        g=@(Theta)sin(Theta)*sum(J'.*exp(i*k0*cos(Theta)*x));
        Patternt=[];
        Patternp=[];
        theta=-pi:.005:pi;        
        for th=theta
            Patternt=[Patternt,20*log10(abs(g(th)))];
            Patternp=[Patternp,20*log10(abs(g(pi/2)))];            
        end
        figure(10+F)
        subplot(1,2,1)
        Temp=find(Patternt<(max(Patternt)-15));
        Patternt(Temp)=[];
        theta(Temp)=[];
        Polar(theta,Patternt-min(Patternt),[-10,-5,0],[5,10,15])
        title(['Raidiation Pattenr for Linear antenna \Phi=0 L = ',num2str(C),' \lambda'])                
        subplot(1,2,2)
        theta=-pi:.005:pi;        
        Polar(theta,Patternp,[0],[Patternp(1)])
        title(['Raidiation Pattenr for Linear antenna \theta=\pi/2 L = ',num2str(C),' \lambda'])                
        disp(['Input Impedance of Linear antenna L = ',num2str(C),' lambda N = ',num2str(N),' : ',num2str(1/Y((N+1)/2,(N+1)/2))]);

    end
        Imped=[Imped,J((N+1)/2)];
end
% figure(3)
% plot(c,abs(real(Imped))*1000)
% title('Input Impedance of Current fed linear antenna d=lambda/1000')
% xlabel('L/lamda')
% ylabel('G(millimhos)')
% figure(4)
% plot(c,imag(Imped)*1000)
% title('Input Impedance of Current fed linear antenna d=lambda/1000')
% xlabel('L/lamda')
% ylabel('B(millimhos)')
end