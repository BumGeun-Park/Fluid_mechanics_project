clear all
clc
%�й� �� 3�ڸ�
P=0;
Q=4;
R=9;
N=10; %�ð� ����(N��)

density_material=(P+1)*100; %��ü�� �е�
Diameter=(Q+1)*0.01; %����
density_air=1.25; %������ �е�
viscosity=1.5*10^-5; %������ �������
mass=(pi/6)*Diameter^3*density_material; %��ü�� ����
g=9.81;%�߷°��ӵ�

%Re�� (0.1���� 10^6����)�� ��ȭ�� ���� CD���� ��ȭ�� �׷�����.
Re_1=0.1:10^6;
C_D_1=(24./Re_1)+...
    (2.6*(Re_1/5)./(1+(Re_1/5).^1.52))+...
    0.411*((Re_1/(2.63*10^5)).^(-7.94))./(1+(Re_1/(2.63*10^5)).^(-8.00))+...
    0.25*(Re_1/10^6)./(1+(Re_1/10^6));
loglog(Re_1,C_D_1),hold on,title('�׷°��'),xlabel('Re��'),ylabel('�׷°��')
axis([0.1 10^6 0.06 400])

U=10^-20;%�ʱ�ӵ�(0���� �Է��ϸ� �����߻�)
dt=0.0001%��Ÿ t
N=N*10^4;% ��������
velocity=zeros(1,(N/10^4)/dt);
n=1;

while n<=(N/10^4)/dt %N��
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
terminal_velocity=velocity(1,(N/10^4)/dt) %���� ���ϼӵ�
plot(time,velocity),title('�ӵ�'),xlabel('�ð�'),ylabel('�ӷ�')

[A,B]=min(abs(velocity-0.99*terminal_velocity));
Time=time(1,B)%�ӵ��� ���� ���ϼӵ��� 0.99%�� ���� ���Ͻð�
