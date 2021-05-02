clear all
clc
[x,y]=meshgrid(-15:0.02:15,0:0.02:15);

%학번 B617052
U=1;% 0+1
a=8;% 5+1
m=4;% 2+1

theta_1=x+y*i-a;
theta_2=x+y*i+a;
Theta_1=angle(theta_1);%세타1
Theta_2=angle(theta_2);%세타2
Theta_0=Theta_1-Theta_2;%세타1-세타

psi=U*y-(m/(2*pi))*Theta_0;%유동함수  psi=y-(3/(2*pi))*(세타1-세타2)
Psi=(psi+abs(psi))/2;
%유동함수가 0인 표면 아래의 지점은 음수값을 가지므로 그 값들을 0으로 만들기 위해

contour(x,y,Psi,70,'b'),hold on
contour(x,y,psi,[0 0],'r','ShowText','on'),hold on %유동함수가 0인 지점
axis([-10 10 0 10]),title('유선')

figure(2);
contour(x,y,psi,[0.0000000000000000000001 0.0000000000000000000001],'r')
axis([-10 10 0 10]),title('유동함수가 0인 유선 _(_둔_덕_)')

syms q w
n=1;
N=ceil(sqrt((m*a)/(pi*U)+a^2));%계산범위(-N,N)
psi_0_x=ones(1,N*2000-1);%유동함수가 0인 지점들의 x좌표
psi_0_y=zeros(1,N*2000-1);%유동함수가 0인 지점들의 y좌표
while -N+0.001*n<N
    S=vpasolve(U*w-(m/(2*pi))*(angle(q+w*i-a)-angle(q+w*i+a))==0,[q w],[-N+0.001*n,1]);   
    psi_0_x(1,n)=S.q;
    psi_0_y(1,n)=S.w;
    n=n+1
end
psi_0_y=round(psi_0_y*10^30)/10^30;%10^-30 이하의 수는 0 취급
u=U+(m/(4*pi))*(2*(psi_0_x+a)./((psi_0_x+a).^2+psi_0_y.^2)-2*(psi_0_x-a)./((psi_0_x-a).^2+psi_0_y.^2));% x방향 속도
u=u.*logical(psi_0_y)%x방향 속도중에 표면에 접한부분의 x방향속도
v=(m/(4*pi))*(2*psi_0_y./((psi_0_x+a).^2+psi_0_y.^2)-2*psi_0_y./((psi_0_x-a).^2+psi_0_y.^2));% y방향 속도
velocity=sqrt(u.^2+v.^2);%속력
pressure=101.3*10^3+(1/2)*(U^2-velocity.^2);%절대압력
Pressure=pressure.*logical(pressure-0.5);%표면에 접한부분의 압력
st=-N+0.001:0.001:N-0.001;

figure(3)
plot(st,Pressure);%압력분포 그림
title('둔덕에서의 압력분포')

u_1=u(:,((2000*N-1)-length(find(u)))/2+1:(2000*N-1)-((2000*N-1)-length(find(u)))/2);%106~13894행 
v_1=v(:,((2000*N-1)-length(find(v))-1)/2+1:(2000*N-1)-((2000*N-1)-length(find(v))-1)/2);%106~13894행
velocity_1=velocity(:,((2000*N-1)-length(find(velocity)))/2+1:(2000*N-1)-((2000*N-1)-length(find(velocity)))/2);%106~13894행
Pressure_1=Pressure(:,((2000*N-1)-length(find(u)))/2+1:(2000*N-1)-((2000*N-1)-length(find(u)))/2);%106~13894행
psi_0_x_1=psi_0_x(:,((2000*N-1)-length(find(psi_0_y)))/2+1:(2000*N-1)-((2000*N-1)-length(find(psi_0_y)))/2);%106~13894행
psi_0_y_1=psi_0_y(:,((2000*N-1)-length(find(psi_0_y)))/2+1:(2000*N-1)-((2000*N-1)-length(find(psi_0_y)))/2);%106~13894행

unit_vector=zeros(2,length(find(u)));%속도벡터와 직교하는 단위벡터(둔덕과 수직)
unit_vector(1,:)=v_1./velocity_1;% X성분
unit_vector(2,:)=-u_1./velocity_1;% Y성분
P_x=unit_vector(1,:).*Pressure_1;%압력의 X방향성분
P_y=unit_vector(2,:).*Pressure_1;%압력의 Y방향성분
dx=abs(gradient(psi_0_x_1));
dy=abs(gradient(psi_0_y_1));

F_x=round(100*sum(P_x.*dy))/100 %표면을 따라 작용하는 X방향 압력의 합 = X방향 힘, 소수점 3자리에서 반올림
F_y=round(sum(P_y.*dx)/100)*100 %표면을 따라 작용하는 Y방향 압력의 합 = Y방향 힘, 소수점 3자리에서 반올림

figure(4)
plot(st,velocity)
title('둔덕에서의 속도분포')