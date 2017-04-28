function N = removenan( L,filter,Num)
% Apply filers on a sample L.
if nargin < 3 % if less that 3 parameters use all sample
    Num = length(L);
end 
N=[];
for i=1:Num
    x=L(i,:);
    N_col=length(x);
    if ~isnan(x(2)) && filter(x(3),x(4),x(6),x(8))
        if N_col>7 && isnan(x(8))
            continue;
        elseif N_col>8 && isnan(x(9))
            continue;
        end
        N=[N;x];
    end
end

end

