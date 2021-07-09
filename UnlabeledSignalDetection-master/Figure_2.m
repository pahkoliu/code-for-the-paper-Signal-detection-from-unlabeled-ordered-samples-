% clear;clc;clf;close all;
n = 800;
m1 = 100;
m2 = 200;
m3 = 300;
m4 = 400;
mc = 1000;
t = [-27;-24;-21;-18;-15;-12;-9;-6;-3;0;1;2;5;10];
SNR_eq = exp(t);
P_fa = [1e-3;5e-3;8e-3;1e-2;5e-2;8e-2;1e-1;5e-1;8e-1;1];

% index = load('n_permindex.mat');
% idx0 = index.idx0;

idx0 = randperm(n);
save('n_permindex.mat','idx0');



%% ramp siganl
x = linspace(0,1,n);
x = x'-0.5;
P_success_ramp = zeros(length(SNR_eq),4);

fprintf('ramp 1\n');
m = m1;
P_success_ramp(:,1) =  calcu_P_success(idx0,mc,x,SNR_eq,m);

fprintf('ramp 2\n');           
m = m2;
P_success_ramp(:,2) =  calcu_P_success(idx0,mc,x,SNR_eq,m);

fprintf('ramp 3\n');
m = m3;
P_success_ramp(:,3) =  calcu_P_success(idx0,mc,x,SNR_eq,m);

fprintf('ramp 4\n');
m = m4;
P_success_ramp(:,4) =  calcu_P_success(idx0,mc,x,SNR_eq,m);
save('P_success_ramp.mat','P_success_ramp');

%% sinusoidal signal
f = 1;
xt = linspace(0,1,n);
x = sin(2*pi*f*xt');
P_success_sinusoidal = zeros(length(SNR_eq),4);

fprintf('sinusoidal 1\n');
m = m1;
P_success_sinusoidal(:,1) =  calcu_P_success(idx0,mc,x,SNR_eq,m);

fprintf('sinusoidal 2\n');           
m = m2;
P_success_sinusoidal(:,2) =  calcu_P_success(idx0,mc,x,SNR_eq,m);

fprintf('sinusoidal 3\n');
m = m3;
P_success_sinusoidal(:,3) =  calcu_P_success(idx0,mc,x,SNR_eq,m);

fprintf('sinusoidal 4\n');
m = m4;
P_success_sinusoidal(:,4) =  calcu_P_success(idx0,mc,x,SNR_eq,m);
save('P_success_sinusoidal.mat','P_success_sinusoidal');

%% duobinary signal pulse
W = 1;
xt = linspace(0,1,n);
x = sinc(2*pi*W*xt') + sinc(2*pi*(W*xt'-0.5));

P_success_PRS = zeros(length(SNR_eq),4);

fprintf('duobinary signal pulse 1\n');
m = m1;
P_success_PRS(:,1) =  calcu_P_success(idx0,mc,x,SNR_eq,m);

fprintf('duobinary signal pulse 2\n');           
m = m2;
P_success_PRS(:,2) =  calcu_P_success(idx0,mc,x,SNR_eq,m);

fprintf('duobinary signal pulse 3\n');
m = m3;
P_success_PRS(:,3) =  calcu_P_success(idx0,mc,x,SNR_eq,m);

fprintf('duobinary signal pulse 4\n');
m = m4;
P_success_PRS(:,4) =  calcu_P_success(idx0,mc,x,SNR_eq,m);
save('P_success_PRS.mat','P_success_PRS');

%% chirp signal
f = 1;
xt = linspace(0,1,n);
x = sin(2*pi*f*(xt').^2);

P_success_chirp = zeros(length(SNR_eq),4);

fprintf('chirp 1\n');
m = m1;
P_success_chirp(:,1) =  calcu_P_success(idx0,mc,x,SNR_eq,m);

fprintf('chirp 2\n');           
m = m2;
P_success_chirp(:,2) =  calcu_P_success(idx0,mc,x,SNR_eq,m);

fprintf('chirp 3\n');
m = m3;
P_success_chirp(:,3) =  calcu_P_success(idx0,mc,x,SNR_eq,m);

fprintf('chirp 4\n');
m = m4;
P_success_chirp(:,4) =  calcu_P_success(idx0,mc,x,SNR_eq,m);
save('P_success_chirp.mat','P_success_chirp');


%% Standard plot the figures
alw = 0.75;    % AxesLineWidth
fsz = 10;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
figure(1)
set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties

subplot(2,2,1)
p(1)=plot(log(SNR_eq),P_success_ramp(:,1),'--r*','LineWidth',lw,...
                       'MarkerSize',msz)
                    hold on;
p(2)=plot(log(SNR_eq),P_success_ramp(:,2),'-ro','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on;
p(3)=plot(log(SNR_eq),P_success_ramp(:,3),'-b+','LineWidth',lw,...
                       'MarkerSize',msz)
                    hold on;
p(4)=plot(log(SNR_eq),P_success_ramp(:,4),'-kd','LineWidth',lw,...
                       'MarkerSize',msz)              
xlabel('${\rm log}(\frac{d^2_{\rm min}}{n\sigma^2 H(\rho)})$'...
    ,'interpreter','latex','Fontsize',fsz+2)

ylabel('${\rm Pr}(\bf {\hat S}={\bf S})$'...
    ,'interpreter','latex','Fontsize',fsz+2)
text(-18,0.9,'ramp signal')

subplot(2,2,2)
plot(log(SNR_eq),P_success_sinusoidal(:,1),'--r*','LineWidth',lw,...
                       'MarkerSize',msz)
                    hold on;
plot(log(SNR_eq),P_success_sinusoidal(:,2),'-ro','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on;
plot(log(SNR_eq),P_success_sinusoidal(:,3),'-b+','LineWidth',lw,...
                       'MarkerSize',msz)
                    hold on;
plot(log(SNR_eq),P_success_sinusoidal(:,4),'-kd','LineWidth',lw,...
                       'MarkerSize',msz)              
xlabel('${\rm log}(\frac{d^2_{\rm min}}{n\sigma^2 H(\rho)})$'...
    ,'interpreter','latex','Fontsize',fsz+2)
ylabel('${\rm Pr}(\bf {\hat S}={\bf S})$'...
    ,'interpreter','latex','Fontsize',fsz+2)
text(-18,0.9,'sinusoidal signal')


subplot(2,2,3)
plot(log(SNR_eq),P_success_PRS(:,1),'--r*','LineWidth',lw,...
                       'MarkerSize',msz)
                    hold on;
plot(log(SNR_eq),P_success_PRS(:,2),'-ro','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on;
plot(log(SNR_eq),P_success_PRS(:,3),'-b+','LineWidth',lw,...
                       'MarkerSize',msz)
                    hold on;
plot(log(SNR_eq),P_success_PRS(:,4),'-kd','LineWidth',lw,...
                       'MarkerSize',msz)              
xlabel('${\rm log}(\frac{d^2_{\rm min}}{n\sigma^2 H(\rho)})$'...
    ,'interpreter','latex','Fontsize',fsz+2)
ylabel('${\rm Pr}(\bf {\hat S}={\bf S})$'...
    ,'interpreter','latex','Fontsize',fsz+2)
text(-18,0.9,'duobinary signal pulse')


subplot(2,2,4)
plot(log(SNR_eq),P_success_chirp(:,1),'--r*','LineWidth',lw,...
                       'MarkerSize',msz)
                    hold on;
plot(log(SNR_eq),P_success_chirp(:,2),'-ro','LineWidth',lw,...
                       'MarkerSize',msz)
                   hold on;
plot(log(SNR_eq),P_success_chirp(:,3),'-b+','LineWidth',lw,...
                       'MarkerSize',msz)
                    hold on;
plot(log(SNR_eq),P_success_chirp(:,4),'-kd','LineWidth',lw,...
                       'MarkerSize',msz)              
xlabel('${\rm log}(\frac{d^2_{\rm min}}{n\sigma^2 H(\rho)})$'...
    ,'interpreter','latex','Fontsize',fsz+2)
ylabel('${\rm Pr}(\bf {\hat S}={\bf S})$'...
    ,'interpreter','latex','Fontsize',fsz+2)
text(-18,0.9,'chirp signal')


% 32 is  the ASCII of blank
str1 = strcat('m = ',32,num2str(m1));
str2 = strcat('m = ',32,num2str(m2));
str3 = strcat('m = ',32,num2str(m3));
str4 = strcat('m = ',32,num2str(m4));

legend(p(1:2),str1,str2)
ah = axes('position',get(gca,'position'),...
            'visible','off');
legend(ah,p(3:4),str3,str4)


