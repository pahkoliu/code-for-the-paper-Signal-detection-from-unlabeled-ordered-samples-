
function P_success = calcu_P_success(idx0,mc,x,SNR_eq,m)
n = length(x);
test = zeros(mc,1);
ID_m = eye(n);
idx = sort(idx0(1:m));
z0 = x(idx);
S0 = ID_m(idx,:);
rho = m/n;
H = -rho*log(rho) - (1-rho)*log(1-rho);
dmin2 = calcu_dmin2(z0,x,idx);
P_success = zeros(length(SNR_eq),1);

for j = 1:length(SNR_eq)
    waitbar(j/length(SNR_eq))
    sigma = sqrt(dmin2/(SNR_eq(j)*n*H));
    for i = 1:mc
        y = sigma.*randn(m,1) + z0;
        S = proS(y,x);
        test(i) = norm(S-S0);
    end
    P_success(j) = 1 - sum(test>0)/mc;
end

end