% X1 , X2 the location of two charge node around a current node for
% calculating PsiF
%Or the location of two current node around a charge point for PsiC+or-
%Z are the location of current or charge node
function Reff=Psi_1(Z1,Z2,dlm,dln,d)
    delta=abs(Z2-Z1)-(dlm+dln)/2;
    alpha=delta+dlm+dln;
    beta=dlm+delta;
    gama=dln+delta;
 Reff=dlm*dln/(alpha*asinh(alpha/d)-beta*asinh(beta/d)-gama*asinh(gama/d)+delta*asinh(delta/d)-sqrt(alpha^2+d^2)+sqrt(beta^2+d^2)+sqrt(gama^2+d^2)-sqrt(delta^2+d^2));    
%     if ((d<.003)|(abs(Z2-Z1)))
%         Reff=dlm*dln/(alpha*asinh(alpha/d)-beta*asinh(beta/d)-gama*asinh(gama/d)+delta*asinh(delta/d)-sqrt(alpha^2+d^2)+sqrt(beta^2+d^2)+sqrt(gama^2+d^2)-sqrt(delta^2+d^2));
% %         f=@(x,y)1./sqrt(d^2+(x-y).^2);
% %         Reff=dlm*dln/dblquad(f,Z1-dlm/2,Z1+dlm/2,Z2-dln/2,Z2+dln/2);
%     else
% %         Reff=dlm*dln/(alpha*asinh(alpha/d)-beta*asinh(beta/d)-gama*asinh(gama/d)+delta*asinh(delta/d)-sqrt(alpha^2+d^2)+sqrt(beta^2+d^2)+sqrt(gama^2+d^2)-sqrt(delta^2+d^2));
%          Reff=sqrt((Z2-Z1)^2+d^2);
%     end
%     f=@(x,y)1./sqrt(d^2+(x-y).^2);
%     Reff=dlm*dln/dblquad(f,Z1-dlm/2,Z1+dlm/2,Z2-dln/2,Z2+dln/2);
end