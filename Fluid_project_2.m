clear all
clc
%학번 뒤 3자리
P=0;
Q=4;
R=9;
N=10; %시간 측정(N초)

density_material=(P+1)*100; %고체의 밀도
Diameter=(Q+1)*0.01; %지름
density_air=1.25; %공기의 밀도
viscosity=1.5*10^-5; %공기의 점성계수
mass=(pi/6)*Diameter^3*density_material; %고체의 질량
g=9.81;%중력가속도

%Re수 (0.1에서 10^6까지)의 변화에 따른 CD값의 변화를 그려본다.
Re_1=0.1:10^6;
C_D_1=(24./Re_1)+...
    (2.6*(Re_1/5)./(1+(Re_1/5).^1.52))+...
    0.411*((Re_1/(2.63*10^5)).^(-7.94))./(1+(Re_1/(2.63*10^5)).^(-8.00))+...
    0.25*(Re_1/10^6)./(1+(Re_1/10^6));
loglog(Re_1,C_D_1),hold on,title('항력계수'),xlabel('Re수'),ylabel('항력계수')
axis([0.1 10^6 0.06 400])

U=10^-20;%초기속도(0으로 입력하면 오류발생)
dt=0.0001%델타 t
N=N*10^4;% 범위설정
velocity=zeros(1,(N/10^4)/dt);
n=1;

while n<=(N/10^4)/dt %N초
Re=(density_air*U*Diameter)/viscosity;
C_D=(24./Re)+...
    (2.6*(Re/5)/(1+(Re/5)^1.52))+...
    0.411*((Re/(2.63*10^5))^(-7.94))/(1+(Re/(2.63*10^5))^(-8.00))+...
    0.25*(Re/10^6)/(1+(Re/10^6));
F_D=C_D*(1/2)*(density_air)*U^2*(pi*Diameter^2/4);
U=U+dt*(g-F_D/mass);
velocity(1,n)=U;
n=n+1
end

time=0:dt:(N-1)/10^4;
figure(2)
terminal_velocity=velocity(1,(N/10^4)/dt) %최종 낙하속도
plot(time,velocity),title('속도'),xlabel('시간'),ylabel('속력')

[A,B]=min(abs(velocity-0.99*terminal_velocity));
Time=time(1,B)%속도가 최종 낙하속도의 0.99%일 때의 낙하시간
