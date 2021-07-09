function P_d = calcu_Pd(x,idx,P_fa,sigma_w,MC)
z0 = x(idx);
m = length(z0);
n = length(x);
rho = m/n;
H = -rho*log(rho) - (1-rho)*log(1-rho);
P_d = zeros(length(P_fa),4);
px = (norm(z0).^2)/m;
miu = sum(z0)/m;

% GLRT_montecarlo
for i = 1:MC
    waitbar(i/MC)
    y0 = sqrt(sigma_w).*randn(m,1);
    y1 = sqrt(sigma_w).*randn(m,1) + z0;
    S0 = proS(y0,x);
    S1 = proS(y1,x);
    x_altmin0 = S0*x;             % MLE of S*x
    x_altmin1 = S1*x;             % MLE of S*x
    T0(i) = 2*(y0'*x_altmin0) - norm(x_altmin0).^2;
    T1(i) = 2*(y1'*x_altmin1) - norm(x_altmin1).^2;
end
T0_sort = sort(T0,'descend');
for i = 1:length(P_fa)
    threshold = T0_sort(round(MC*P_fa(i)));
    P_d(i,1) = sum(T1>threshold)/MC;
end

% NP_theoretical
P_d(:,2) = Q(Qinv(P_fa) - sqrt(m*px./sigma_w));

% energy detector
for i = 1:MC
    waitbar(i/MC)
    y0 = sqrt(sigma_w).*randn(m,1);
    y1 = sqrt(sigma_w).*randn(m,1) + z0;
    T0(i) = (norm(y0).^2)/m;
    T1(i) = (norm(y1).^2)/m;
end
T0_sort = sort(T0,'descend');
for i = 1:length(P_fa)
    threshold = T0_sort(round(MC*P_fa(i)));
    P_d(i,3) = sum(T1>threshold)/MC;
end


% mean detector
P_d(:,4) = Q(Qinv(P_fa) - sqrt(m*(miu.^2)./sigma_w));

end
