function Reff=Psi_2(X1,Z1,X2,Z2,dl1,dl2,d)
% R3=norm([Z1+dl1/2-Z2,X1-X2+dl2/2])
% R2=norm([Z1+dl1/2-Z2,X1-X2-dl2/2])
% R4=norm([Z1-dl1/2-Z2,X1-X2+dl2/2])
% R1=norm([Z1-dl1/2-Z2,X1-X2-dl2/2])
% m=dl1
% l=dl2
% A=R4^2-R3^2+R2^2-R1^2
% e=acos(A/2/l/m)
% v=m*2*l^2*((R4^2-R3^2-m^2)+A*(R2^2-R3^2-l^2)/(4*l^2*m^2-A^2))
% u=l*2*m^2*((R2^2-R3^2-l^2)+A*(R4^2-R3^2-m^2)/(4*l^2*m^2-A^2))
% D=R3^2-u^2-v^2+2*u*v*cos(e)
% d=sqrt(D);
% O=atan((D*cos(e)+(u+l)*(v+m)*(sin(e))^2)/d/R1/sin(e))+...
%                 -atan((D*cos(e)+u*v*(sin(e))^2)/d/R2/sin(e))+...
%                             atan((D*cos(e)+u*v*(sin(e))^2)/d/R3/sin(e))+...
%                                         atan((D*cos(e)+u*(v+m)*(sin(e))^2)/d/R4/sin(e))
% M=2*cos(e)*((u+l)*atanh(m/(R2+R1))+(v+m)*atanh(1/(R1+R4))-u*atanh(m/(R3+R4))-v*atanh(1/(R2+R3)))^(-O*d*cot(e));
% Reff=dl1*dl2/M
if (sqrt((X1-X2)^2+(Z1-Z2)^2)<3*sqrt(dl1^2+dl2^2))
    f=@(x,z)1./sqrt((X1-X2-x).^2+(Z1+z-Z2).^2+d^2);
    Reff=dl1*dl2/dblquad(f,-dl2/2,dl2/2,-dl1/2,dl1/2);
else
    Reff=sqrt((X1-X2)^2+(Z1-Z2)^2+d^2);
end
end