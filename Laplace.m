
function x = Laplace(mean, variance, Nbr)

u = rand(1,Nbr);
u = u - 0.5;
b = sqrt(variance/2);

x = mean - b * sign(u).*log(1-2*abs(u));

end

