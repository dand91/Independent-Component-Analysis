function [E1] = EOne(P)

[s1,s2] = size(P);

Etemp1 = zeros(s1,s2);
Etemp2 = Etemp1;

n = 1;
for i = 1:s1;
    for j = 1:s2;
        e1(n) = abs(P(i,j))/max(abs( P(:,j)) );
        e2(n) = abs(P(i,j))/max(abs( P(i,:)) );
        n = n + 1;
    end
end

E1 = ((sum(e1)) - 2) + ((sum(e2)) - 2) ;

end