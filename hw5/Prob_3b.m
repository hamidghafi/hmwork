function Prob_3b()
clc
close all
f=1e9;u0=4*pi*1e-7;e0=1/36/pi*1e-9;
w=2*pi*f;
Imped=[];
lambda=3e8/f;
k0=w/3e8;
%first the antenna should be discretize in N1 N2 N3 part for every
%parasitic element
%For Reflector
for N1=35:8:65

    L1=0.475*lambda;
    dl1=L1/(N1+2);
    X=zeros(N1,1);
    Y=zeros(N1,1);
    Z=(linspace(dl1,L1-dl1,N1))';
    DL=ones(N1,1)*dl1;
    Num=ones(N1,1);
    %For Dipole
%     N21=floor(N1/2);
%     L2=0.23*lambda;
%     dl21=L2/(N21+1);
%     X=[X;ones(N21,1)*lambda/4];
%     Y=[Y;zeros(N21,1)];
%     Z=[Z;(linspace(.24*lambda+.015/2*lambda,L2-dl21+.23*lambda+.015/2*lambda,N21))'];
%     DL=[DL;ones(N21,1)*dl21];
% 
%     N22=floor(N1/2);
%     dl22=L2/(N22+1);
%     X=[X;ones(N22,1)*lambda/4];
%     Y=[Y;zeros(N22,1)];
%     Z=[Z;(linspace(dl22+.015/2*lambda,.22*lambda+.015/2*lambda,N22))'];
%     DL=[DL;ones(N22,1)*dl22];
%     N2=N21+N22;

    N2=N1;
    L2=0.46*lambda;
    dl2=L2/(N2+2);
    X=[X;ones(N2,1)*lambda/4];
    Y=[Y;zeros(N2,1)];
    Z=[Z;(linspace(dl2+.015*lambda/2,L2-dl2+.015*lambda/2,N2))'];
    DL=[DL;ones(N2,1)*dl2];
    Num=[Num;ones(N2,1)*2];
    %For Dirctor
    N3=N1;
    L3=0.44*lambda;
    dl3=L3/(N3+2);
    X=[X;ones(N3,1)*lambda*0.56];
    Y=[Y;zeros(N3,1)];
    Z=[Z;(linspace(dl3+.035/2*lambda,L3-dl3+.035/2*lambda,N3))'];
    DL=[DL;ones(N3,1)*dl3];
    Num=[Num;ones(N3,1)*3];
    %For Boom
    N4=N1;
    L4=.56*lambda;
    dl4=L4/(N4+2);
    X=[X;(linspace(0,L4,N4))'];
    Y=[Y;ones(N4,1)*0];
    Z=[Z;ones(N4,1)*(.23*lambda+.015/2*lambda)];
    DL=[DL;ones(N4,1)*dl4];
    Num=[Num;ones(N4,1)*4];
    figure(4)
    plot(X,Z,'.')
    % %------------------------------------------------------------------
    N=N1+N2+N3+N4;
    d=lambda/1000;
    d1=d;
    for m=1:N
        for n=1:N
            if (Num(n)<4 & Num(m)<4)
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
                
                
            elseif(Num(m)==4 & Num(n)==4)
                Reff=Psi_1(X(n),X(m),DL(m),DL(n),d);
                psi(m,n)=DL(m)*DL(n)*1/4/pi/Reff*exp(-i*k0*Reff);
                
                Reffpp=Psi_1(X(n)+DL(n)/2,X(m)+DL(m)/2,DL(m),DL(n),d);
                Psipp(m,n)=1/4/pi/Reffpp*exp(-i*k0*Reffpp);

                Reffpn=Psi_1(X(n)+DL(n)/2,X(m)-DL(n)/2,DL(m),DL(n),d);
                Psipn(m,n)=1/4/pi/Reffpn*exp(-i*k0*Reffpn);

                Reffnp=Psi_1(X(n)-DL(m)/2,X(m)+DL(m)/2,DL(m),DL(n),d);
                Psinp(m,n)=1/4/pi/Reffnp*exp(-i*k0*Reffnp);
                
                Reffnn=Psi_1(X(n),X(m),DL(m),DL(n),d);
                Psinn(m,n)=1/4/pi/Reffnn*exp(-i*k0*Reffnn);

            else
%                  Reff=Psi_2(X(n),Z(n),X(m),Z(m),DL(n),DL(m),d1);
%                  psi(m,n)=DL(m)*DL(n)*1/4/pi/Reff*exp(-i*k0*Reff);          
                psi(m,n)=0;
                Reffpp=Psi_2(X(n)+DL(n)/2,Z(n),X(m),Z(m)+DL(m)/2,DL(n),DL(m),d1);
                Psipp(m,n)=1/4/pi/Reffpp*exp(-i*k0*Reffpp);
                Reffpn=Psi_2(X(n)+DL(n)/2,Z(n),X(m),Z(m)-DL(m)/2,DL(n),DL(m),d1);
                Psipn(m,n)=1/4/pi/Reffpn*exp(-i*k0*Reffpn);

                Reffnp=Psi_2(X(n)-DL(n)/2,Z(n),X(m),Z(m)+DL(m)/2,DL(n),DL(m),d1);
                Psinp(m,n)=1/4/pi/Reffnp*exp(-i*k0*Reffnp);
                
                Reffnn=Psi_2(X(n)-DL(n)/2,Z(n),X(m),Z(m)-DL(m)/2,DL(n),DL(m),d1);
                Psinn(m,n)=1/4/pi/Reffnn*exp(-i*k0*Reffnn);
            end

            Z_an(m,n)=i*w*u0*psi(m,n)+(Psipp(m,n)+Psinn(m,n)-Psinp(m,n)-Psipn(m,n))/i/w/e0;
        end
    end
%     V=[zeros(1,N1),zeros(1,(N1-1)/2),1,zeros(1,(N1-1)/2),zeros(1,N3),zeros(1,N4)];
    V=[zeros(1,N1),zeros(1,(N1-1)/2),1,zeros(1,(N1-1)/2),zeros(1,N3),zeros(1,N4)];
    % V=[zeros(1,N1),zeros(1,(N2-1)/2),1,zeros(1,(N2-1)/2)];
    % V=[zeros(1,((N1-1)/2)),1,zeros(1,((N1-1)/2))];
    Y=inv(Z_an);
    J=Y*V';
    %****************************************************************
    J1=[0;J(1:N1);0];
    x=linspace(0,L1,length(J1))/lambda;
    figure(1)
    subplot(4,1,1)
    hold on
    plot(x,abs(J1)*1000,'b')
    J2=[0;J(N1+1:N2+N1);0];
    x=linspace(0,L2,length(J2))/lambda;
    title(['Current distribution for Yagi antenna d=\lambda/1000 L = .475\lambda Reflector'],'FontSize',13)
    xlabel('l/\labmda','FontSize',12)
    ylabel('Milliamperes','FontSize',12)
    subplot(4,1,2)
    hold on
    plot(x,abs(J2)*1000,'b-')
    title(['Current distribution for Yagi antenna d=\lambda/1000 L = .46\lambda Dipole'],'FontSize',13)
    xlabel('l/\labmda','FontSize',12)
    ylabel('Milliamperes','FontSize',12)
    subplot(4,1,3)
    hold on
    J3=[0;J(N1+N2+1:N1+N2+N3);0];
    x=linspace(0,L3,length(J3))/lambda;
    plot(x,abs(J3)*1000,'c-.')
    title(['Current distribution for Yagi antenna d=\lambda/1000 L = .44\lambda Director'],'FontSize',13)
    xlabel('l/\labmda','FontSize',12)
    ylabel('Milliamperes','FontSize',12)
    subplot(4,1,4)
    hold on
    J4=[0;J(N1+N2+N3+1:N);0];
    x=linspace(0,L4,length(J4))/lambda;
    plot(x,abs(J4)*1000,'c-.')
    title(['Current distribution for Yagi antenna d=\lambda/1000 L = .56\lambda Director'],'FontSize',13)
    xlabel('l/\labmda','FontSize',12)
    ylabel('Milliamperes','FontSize',12)
    %****************************************************************
    X1=X(1:N-N4);
    X2=X(N-N4+1:N);
    Z1=Z(1:N-N4);
    Z2=Z(N-N4+1:N);
    J1=J(1:N-N4);
    J2=J(N-N4+1:N);
    g=@(Theta,Phi)sin(Theta)*sum(J1.*exp(i*k0*(cos(Theta)*Z1+sin(Theta)*cos(Phi)*X1)))+cos(Theta)*cos(Phi)*sum(J2.*exp(i*k0*(cos(Theta)*Z2+sin(Theta)*cos(Phi)*X2)));
    %for X direction of Element
    %g1=@(Theta,Phi)cos(Theta)*cos(Phi)*sum(J2.*exp(i*k0*(cos(Theta)*Z2+sin(Theta)*cos(Phi)*X2)));
    g1=@(Theta,Phi)-sin(Phi)*sum(J2.*exp(i*k0*(cos(Theta)*Z2+sin(Theta)*cos(Phi)*X2)));
    Patternt1=[];
    Patternt2=[];
    theta=-pi:.005:pi;
    phi=linspace(0,2*pi,100);
    theta1=linspace(0,pi,100);
    P=0;
    for th=theta
        Patternt1=[Patternt1,10*log10((abs(g(th,0)))^2)];
        Patternt2=[Patternt2,10*log10((abs(g(pi/2,th)))^2)+10*log10((abs(g1(pi/2,th)))^2)];        
    end
    for th=theta1
        for ph=phi
            P=((abs(g(th,ph)))^2+(abs(g1(th,ph)))^2)*(2*pi/100)*(pi/100)*sin(th)+P;
        end
    end
    figure(2)
    G=10*log10(4*pi/P)+max(Patternt1);
    disp(['Max Gain for Yagi antenna N = ',num2str(N),': ',num2str(G)])
    Temp=find(Patternt1<(max(Patternt1)-45));
    Patternt1(Temp)=[];
    theta(Temp)=[];
    Polar(theta,Patternt1-min(Patternt1),[-35,-25,-15,-5,0],[5,15,25,40,45])
    title('Raidiation Pattern for 3 Element Yagi Antenna','FontSize',13)
    hold on
    disp(['Input Impedace of Yagi antenna L = ',num2str(1/J(N1+(N1+1)/2))])    
    figure(3)
    theta=[-pi:.005:pi];
    Temp=find(Patternt2<(max(Patternt2)-15));
    Patternt2(Temp)=[];
    theta(Temp)=[];
    Polar(theta,Patternt2-min(Patternt2),[-10,-5,0],[5,10,15])
    hold on
end
end