function W = GHA(x,nbrIterations,alpha)

disp('NbrIterations');
disp('alpha');
disp(nbrIterations);
disp(alpha);


% Initiatee W.

[~,x2] = size(x);
W = eye(x2);

% GHA algorithm.

k = 0;
x = x.';
min = 10;

while k < nbrIterations;
    
    y = W * x;
    W = W + alpha * ( y * x.' - tril( y * y.', 0 ) * W );
    
    % Normalize vectors
    
    for i = 1:x2
        W(i,:) = W(i,:)/norm(W(i,:));
    end
    k = k + 1;
    
end
disp('min');
disp(min);
W = W;

end