% calculate dmin^2
% Input: the original n*1 signal x, the specific m*1 signal z0, the index
% of z0 in x
% Output: the minimum squared distance of dmin^2 between z0 and other
% possible m*1 vector z that is selected from the signal x
function dmin2 = calcu_dmin2(z0,x,idx)
m = length(z0);
n = length(x);
for i=1:m
    if i==1
        k1=0;
    else
        k1=idx(i-1);
    end
    
    if i==m
        k2=n+1;
    else
        k2 = idx(i+1);
    end
    k3=0;
    for j=(k1+1):(k2-1)
        k3=k3+1;
        d1(k3)=abs(z0(i)-x(j)).^2;
    end
    indexgq0=find(d1>0);
    if(isempty(indexgq0))
        dmin2_possible(i) = 0;
    else
        dmin2_possible(i) = min(d1(indexgq0));
    end
end
indexgq0 = find(dmin2_possible>0);
dmin2 = min(dmin2_possible(indexgq0));
end

