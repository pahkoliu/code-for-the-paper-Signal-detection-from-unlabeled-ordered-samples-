%projection on S
function S=proS(y,z)
m = size(y,1);
n = size(z,1);
%fill the dynamic programming table
T = 1./zeros(m,n);
    for i=1:m
        T(i,i) = norm(y(1:i)-z(1:i)).^2;
    end
    
    for i=2:n
%         yz_norm = zeros(i,1);
         yz_norm = (y(1).*ones(i,1)-z(1:i)).^2;
%         for j=1:i
%             yz_norm(j) = norm(y(1)-z(j)).^2;
%         end
        T(1,i) = min(yz_norm);
    end
    
    for i=2:m
        for j=i+1:n
            T(i,j) = min( ( T(i-1,j-1) + norm(y(i)-z(j)).^2 ) , T(i,j-1) );
        end
    end
    
    %find the indices of those m elements of z that are matched to y
    S = zeros(m,n);
    idx = n+1;
    for i=m:-1:1
        idx0 = idx-1;
        for j=i:idx0
            if T(i,j)==T(i,idx0)
                idx = j;
                S(i,idx) = 1;
                break;
            end
        end
    end
    