% x exist vs x do not exist
% y = S*x + w
clear;clc;clf;close all;
n = 800;
m = 100;
index = load('n_permindex.mat');
idx0 = index.idx0;
idx = sort(idx0(1:m));
MC = 1000;
% P_fa = 0.05:0.05:1;
P_fa = linspace(5e-4,0.001,10);
P_fa = [5e-4,6e-4,7e-4,8e-4,1e-3,2e-3,3e-3,1e-2,2e-2,3e-2,3e-2,1e-1,2e-1,3e-1];
P_fa = [1e-3;1e-2;5e-2;8e-2;1e-1;5e-1;8e-1;1];
P_fa = [6e-4;8e-4;1e-3;1e-2;5e-2;8e-2;1e-1;5e-1;8e-1;1];
P_fa = [5e-4;6e-4;8e-4;1e-3;6e-3;1e-2;5e-2;8e-2;1e-1;5e-1;8e-1;1];
P_fa = [5e-4;1e-2;8e-2;1e-1;2e-1;3e-1;4e-1;5e-1;8e-1;1];

% d_c = 4;
rho = m/n;
H = -rho*log(rho) - (1-rho)*log(1-rho);

%% Standard plot the figures 
alw = 0.75;    % AxesLineWidth
fsz = 10;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
figure(1)
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties


%% ramp siganl
x = linspace(0,1,n);
x = x'-0.5;
z0 = x(idx);
t = -19;
SNR_eq = exp(t);
dmin2 = calcu_dmin2(z0,x,idx);
% sigma_w = z0'*z0/d_c;
% kappa = dmin2/(n*H*sigma_w)
sigma_w = dmin2/(SNR_eq*n*H);
fprintf('ramp signal\n');    
P_d = calcu_Pd(x,idx,P_fa,sigma_w,MC);
subplot(2,2,1)
plot(P_fa, P_d(:,1),'--r*','LineWidth',lw,...
                       'MarkerSize',msz)
                    hold on;
plot(P_fa, P_d(:,2),'-bo','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on;
plot(P_fa, P_d(:,3),'-r+','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on;                   
plot(P_fa, P_d(:,4),'-kd','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on; 
grid
xlabel('P_{\rm FA}')
ylabel('P_{\rm D}')
text(0.5,0.3,'ramp signal')
text(0.4,0.7,['\eta = ',num2str(t)])


%% sinusoidal signal
f = 1;
xt = linspace(0,1,n);
x = sin(2*pi*f*xt');
z0 = x(idx);
t = -30;
SNR_eq = exp(t); 
dmin2 = calcu_dmin2(z0,x,idx);
sigma_w = dmin2/(SNR_eq*n*H);
% sigma_w = z0'*z0/d_c;
% kappa = dmin2/(n*H*sigma_w)
fprintf('sinusoidal signal\n');   
P_d = calcu_Pd(x,idx,P_fa,sigma_w,MC);
subplot(2,2,2)
plot(P_fa, P_d(:,1),'--r*','LineWidth',lw,...
                       'MarkerSize',msz)
                    hold on;
plot(P_fa, P_d(:,2),'-bo','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on;
plot(P_fa, P_d(:,3),'-r+','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on;                   
plot(P_fa, P_d(:,4),'-kd','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on; 
grid
xlabel('P_{\rm FA}')
ylabel('P_{\rm D}')
text(0.4,0.2,'sinusoidal signal')
text(0.4,0.7,['\eta = ',num2str(t)])


%% duobinary signal pulse
W = 1;
xt = linspace(0,1,n);
x = sinc(2*pi*W*xt') + sinc(2*pi*(W*xt'-0.5));
z0 = x(idx);
t = -28;
SNR_eq = exp(t);
dmin2 = calcu_dmin2(z0,x,idx);
sigma_w = dmin2/(SNR_eq*n*H);
% sigma_w = z0'*z0/d_c;
% kappa = dmin2/(n*H*sigma_w)
fprintf('duobinary signal pulse\n');   
P_d = calcu_Pd(x,idx,P_fa,sigma_w,MC);
subplot(2,2,3)
plot(P_fa, P_d(:,1),'--r*','LineWidth',lw,...
                       'MarkerSize',msz)
                    hold on;
plot(P_fa, P_d(:,2),'-bo','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on;
plot(P_fa, P_d(:,3),'-r+','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on;                   
plot(P_fa, P_d(:,4),'-kd','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on; 
grid
axis([0 1 0 1]);
xlabel('P_{\rm FA}')
ylabel('P_{\rm D}')
text(0.2,0.4,'duobinary signal pulse')
text(0.4,0.7,['\eta = ',num2str(t)])


%% chirp signal
f = 1;
xt = linspace(0,1,n);
x = sin(2*pi*f*(xt').^2);
z0 = x(idx);
t = -26;
SNR_eq = exp(t); 
dmin2 = calcu_dmin2(z0,x,idx);
sigma_w = dmin2/(SNR_eq*n*H);
% sigma_w = z0'*z0/d_c;
% kappa = dmin2/(n*H*sigma_w)
fprintf('chirp signal\n');        
P_d = calcu_Pd(x,idx,P_fa,sigma_w,MC);
subplot(2,2,4)
p(1)=plot(P_fa, P_d(:,1),'--r*','LineWidth',lw,...
                       'MarkerSize',msz)
                    hold on;
p(2)=plot(P_fa, P_d(:,2),'-bo','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on;
p(3)=plot(P_fa, P_d(:,3),'-r+','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on;                   
p(4)=plot(P_fa, P_d(:,4),'-kd','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on; 
grid
axis([0,1,0,1])
xlabel('P_{\rm FA}')
ylabel('P_{\rm D}')
text(0.4,0.2,'chirp signal')
text(0.4,0.7,['\eta = ',num2str(t)])


%%
legend(p(1:2),'GLRT','clairvoyant')
ah = axes('position',get(gca,'position'),...
            'visible','off');
legend(ah,p(3:4),'energy','mean')        
        
