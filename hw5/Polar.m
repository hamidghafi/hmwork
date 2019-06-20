function Polar(T,R,P,S)
plot(R.*cos(T),R.*sin(T),'LineWidth',2)
hold on
t=[-pi:.01:pi];
for j=1:length(S)
    plot(S(j)*cos(t),S(j)*sin(t),'b')
    axis off
    text(S(j),S(j)/20,num2str(P(j)),'HorizontalAlignment','right')
end
r=max(S);
T=[0,pi/8,pi/4,3*pi/8,pi/2,5*pi/8,6*pi/8,7*pi/8,pi];
for t=T
    plot([-r*cos(t),r*cos(t)],[-r*sin(t),r*sin(t)],'b')
    if (t==pi/2 | t==0 | t==pi)
        text(r*cos(t),r*sin(t),num2str(t/pi*180),'HorizontalAlignment','left')    
    end
end
hold off
end