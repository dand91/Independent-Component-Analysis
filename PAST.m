function [W] = PAST(x,nbrIterations,beta)

disp('NbrIterations:');
disp(nbrIterations);
[~,s2] = size(x);

W = eye(s2) + randn(s2);
P = eye(s2);

x = x.';

for k = 1:nbrIterations;    
        
    % PAST algorithm:
    
   y = W'*x;
   h = P*y;
   g = h / (beta + h' * y);
   P = (1/beta)*(P - g*h');
   e = x - W*y;
   W = W + e*g';
   W = W./norm(W); 
   
end 
end