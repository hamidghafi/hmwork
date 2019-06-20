% X1 , X2 the location of two charge node around a current node for
% calculating PsiF
%Or the location of two current node around a charge point for PsiC+or-
%Z are the location of current or charge node
function Reff=Psi(Z1,Z2,dl,d)
    delta=abs(Z2-Z1)-dl;
    alpha=delta+2*dl;
    beta=dl+delta;
    gama=dl+delta;
    if Z2==Z1
        Reff=dl^2/(alpha*asinh(alpha/d)-beta*asinh(beta/d)-gama*asinh(gama/d)+delta*asinh(delta/d)-sqrt(alpha^2+d^2)+sqrt(beta^2+d^2)+sqrt(gama^2+d^2)-sqrt(delta^2+d^2));
    else
        Reff=dl^2/(alpha*asinh(alpha/d)-beta*asinh(beta/d)-gama*asinh(gama/d)+delta*asinh(delta/d)-sqrt(alpha^2+d^2)+sqrt(beta^2+d^2)+sqrt(gama^2+d^2)-sqrt(delta^2+d^2));
%         Reff=abs(Z2-Z1);
    end
%     f=@(x,y)1./sqrt(d^2+(x-y).^2);
%     Reff=1/dblquad(f,Z1-dl/2,Z1+dl/2,Z2-dl/2,Z2+dl/2);
end