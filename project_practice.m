clear all
clc
[x,y]=meshgrid(-15:0.02:15,0:0.02:15);

%�й� B617052
U=1;% 0+1
a=8;% 5+1
m=4;% 2+1

theta_1=x+y*i-a;
theta_2=x+y*i+a;
Theta_1=angle(theta_1);%��Ÿ1
Theta_2=angle(theta_2);%��Ÿ2
Theta_0=Theta_1-Theta_2;%��Ÿ1-��Ÿ

psi=U*y-(m/(2*pi))*Theta_0;%�����Լ�  psi=y-(3/(2*pi))*(��Ÿ1-��Ÿ2)
Psi=(psi+abs(psi))/2;
%�����Լ��� 0�� ǥ�� �Ʒ��� ������ �������� �����Ƿ� �� ������ 0���� ����� ����

contour(x,y,Psi,70,'b'),hold on
contour(x,y,psi,[0 0],'r','ShowText','on'),hold on %�����Լ��� 0�� ����
axis([-10 10 0 10]),title('����')

figure(2);
contour(x,y,psi,[0.0000000000000000000001 0.0000000000000000000001],'r')
axis([-10 10 0 10]),title('�����Լ��� 0�� ���� _(_��_��_)')

syms q w
n=1;
N=ceil(sqrt((m*a)/(pi*U)+a^2));%������(-N,N)
psi_0_x=ones(1,N*2000-1);%�����Լ��� 0�� �������� x��ǥ
psi_0_y=zeros(1,N*2000-1);%�����Լ��� 0�� �������� y��ǥ
while -N+0.001*n<N
    S=vpasolve(U*w-(m/(2*pi))*(angle(q+w*i-a)-angle(q+w*i+a))==0,[q w],[-N+0.001*n,1]);   
    psi_0_x(1,n)=S.q;
    psi_0_y(1,n)=S.w;
    n=n+1
end
psi_0_y=round(psi_0_y*10^30)/10^30;%10^-30 ������ ���� 0 ���
u=U+(m/(4*pi))*(2*(psi_0_x+a)./((psi_0_x+a).^2+psi_0_y.^2)-2*(psi_0_x-a)./((psi_0_x-a).^2+psi_0_y.^2));% x���� �ӵ�
u=u.*logical(psi_0_y)%x���� �ӵ��߿� ǥ�鿡 ���Ѻκ��� x����ӵ�
v=(m/(4*pi))*(2*psi_0_y./((psi_0_x+a).^2+psi_0_y.^2)-2*psi_0_y./((psi_0_x-a).^2+psi_0_y.^2));% y���� �ӵ�
velocity=sqrt(u.^2+v.^2);%�ӷ�
pressure=101.3*10^3+(1/2)*(U^2-velocity.^2);%����з�
Pressure=pressure.*logical(pressure-0.5);%ǥ�鿡 ���Ѻκ��� �з�
st=-N+0.001:0.001:N-0.001;

figure(3)
plot(st,Pressure);%�зº��� �׸�
title('�д������� �зº���')

u_1=u(:,((2000*N-1)-length(find(u)))/2+1:(2000*N-1)-((2000*N-1)-length(find(u)))/2);%106~13894�� 
v_1=v(:,((2000*N-1)-length(find(v))-1)/2+1:(2000*N-1)-((2000*N-1)-length(find(v))-1)/2);%106~13894��
velocity_1=velocity(:,((2000*N-1)-length(find(velocity)))/2+1:(2000*N-1)-((2000*N-1)-length(find(velocity)))/2);%106~13894��
Pressure_1=Pressure(:,((2000*N-1)-length(find(u)))/2+1:(2000*N-1)-((2000*N-1)-length(find(u)))/2);%106~13894��
psi_0_x_1=psi_0_x(:,((2000*N-1)-length(find(psi_0_y)))/2+1:(2000*N-1)-((2000*N-1)-length(find(psi_0_y)))/2);%106~13894��
psi_0_y_1=psi_0_y(:,((2000*N-1)-length(find(psi_0_y)))/2+1:(2000*N-1)-((2000*N-1)-length(find(psi_0_y)))/2);%106~13894��

unit_vector=zeros(2,length(find(u)));%�ӵ����Ϳ� �����ϴ� ��������(�д��� ����)
unit_vector(1,:)=v_1./velocity_1;% X����
unit_vector(2,:)=-u_1./velocity_1;% Y����
P_x=unit_vector(1,:).*Pressure_1;%�з��� X���⼺��
P_y=unit_vector(2,:).*Pressure_1;%�з��� Y���⼺��
dx=abs(gradient(psi_0_x_1));
dy=abs(gradient(psi_0_y_1));

F_x=round(100*sum(P_x.*dy))/100 %ǥ���� ���� �ۿ��ϴ� X���� �з��� �� = X���� ��, �Ҽ��� 3�ڸ����� �ݿø�
F_y=round(sum(P_y.*dx)/100)*100 %ǥ���� ���� �ۿ��ϴ� Y���� �з��� �� = Y���� ��, �Ҽ��� 3�ڸ����� �ݿø�

figure(4)
plot(st,velocity)
title('�д������� �ӵ�����')