%% 2D-MUSIC 针对Q1，L=1
clear;
tic
% 基本参数
global Na
global L
global T
global Ts
global gamma
global  f0
global Nt
global fs
global c

c=3e8;
Na=86;
L=0.0815;
T=3.2e-5;
Ts=1.25e-7;
gamma=78.986e12;
f0=78.8e9;
Nt=256;
fs=1/Ts;
d=L/85;
%% 平滑操作 增加秩
path="./data_q1.mat";
data=load(path);
% X0= data.Z_noisy;  %86,256
X0= data.Z;  %86,256
X=X0.';
% 窗口大小，自定义参数  64 32

wa=32;
wt=32;
R=linspace(0,10,400);
theta0=linspace(-50,50,400);

% 寻找第一个物体的参数
% R=linspace(6.95,7.05,400);
% theta0=linspace(-1,0.5,400);
% R=linspace(6.98,7.025,400);
% theta0=linspace(-0.7,0,400);
% R=linspace(6.99,7.015,400);
% theta0=linspace(-0.55,-0.15,400);
% R=linspace(6.996,7.012,400);
% theta0=linspace(-0.45,-0.2,400);


% R=linspace(8.14,8.24,400);
% theta0=linspace(-0.5,1,400);
% R=linspace(8.165,8.21,400);
% theta0=linspace(0,0.65,400);
% R=linspace(8.175,8.2,400);
% theta0=linspace(0.15,0.5,400);



na=Na-wa+1;
nt=Nt-wt+1;
Y=zeros(wt*wa,nt*na);
t=1;
for i=1:nt
    for j=1:na
        temp1=X(i:i+wt-1,j:j+wa-1);
        temp2=reshape(temp1,wa*wt,1);
        Y(:,t)=temp2;
        t=t+1;
    end
end

J=flip(eye(wa*wt));
Rxx=( Y*Y' +  J*(Y*Y')*J )   ./(2*nt*na);
[S,V,~]=svd(Rxx);

%%
% 判断几个障碍物   部分组  和完备集  的结论应当一致
%%

U_noise=S(:,3:wa*wt);

% 
% aa=vrtheta(1,pi/2,wt,wa)


idx1=(0:wt-1)';
idx2=(0:wa-1)';

% v_R=exp(1j*4*pi*gamma*T*R.*idx1./(c*T*fs) )  ;
% v_theta=exp(1j*2*pi*(L/85)*sin(theta)/(c/fs).*idx2);
% y=kron(v_R, v_theta);

V_R=@(R) exp(1j*4*pi*gamma*T*R.*idx1./(c*T*fs) )  ;
V_theta=@(theta0) exp(1j*2*pi/(c/f0)*d*sind(theta0).*idx2);
% tao_mk=@(R,theta) 2*(R + )  /c
% 
% f1=


% theta0=linspace(-1.5,1,20);
% R=linspace(6.96,7.04,20);



q1=length(theta0);
q2=length(R);

result=zeros(q1,q2);
parfor i =1:q1
    fprintf('当前进度：%d / %d \n', i,q1)
    for j=1:q2
        tt=1./( kron( V_theta( theta0(i)), V_R( R(j)) )' *( U_noise * U_noise') *  kron( V_theta( theta0(i)), V_R( R(j)) ) );
        result(i,j)=abs(tt);
    end
end
figure
imagesc(R,theta0,result)
colorbar
xlabel('距离/ m')
ylabel('角度/ °')
toc

% save q1_4_200 R theta0 result

%%
tmax=max(max(result));
[R_idx,theta_idx] = find(result==tmax);
Rf=zeros(length(R_idx),1);
thetaf=zeros(length(R_idx),1);
for i=1:length(R_idx)
    Rf(i)=R(R_idx(i));
    thetaf(i)=theta0(theta_idx(i));
end


% 转化为二维坐标，求几何平均
x=zeros(length(R_idx),1);
y=zeros(length(R_idx),1);
for i=1:length(R_idx)
    x(i)=Rf(i)*sind(thetaf(i));
    y(i)=Rf(i)*cosd(thetaf(i));
end
xx=mean(x);
yy=mean(y);

% [t,r]=cart2pol(xx,yy);

[Rf,thetaf]

save q1_32_32 R theta0 result Rf thetaf
