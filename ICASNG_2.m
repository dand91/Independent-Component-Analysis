function [A,W,xout] = ICASNG()

warning('off','all')

nbr = 2;

play = 1;

[e(1,:), Fs1] = wavread('3.wav'); %#ok<DWVRD>
[e(2,:), Fs2] = wavread('5.wav'); %#ok<DWVRD>

%[e, Fs1] = wavread('2x2_Conv_Mixture.wav'); %#ok<DWVRD>
%e = e.';

%SampleSize = 1000;
%e(1,:) = Laplace(0,2,SampleSize);
%e(2,:) = Laplace(0,2,SampleSize);


% Mix data

A = randn(nbr);
Aref = inv(A);
x = A * e;

% Center data.

[Es1,Es2] = size(e);
MeanS = zeros(1,Es1);

for n = 1:Es1;
    
    MeanS(n) = mean(x(n,:));
    x(n,:) = (x(n,:) - MeanS(n));
    
end

%Calculate covariance matrix.

C = (1 / (Es2-1)) * (x * x.');

%Find eigenvectors and eigenvalues.

[U,D] = eig(C);

% Whiten data

z = U.'*x;

% Initiate variables

MmaxIterations = 10^5;
NmaxIterations = 10^4;

m = 1000;
  

while(m < MmaxIterations)

    % Change step size

    mu = 1/m;
        
        n = 1;
        W = eye(Es1);

        while(n < NmaxIterations);
            
           %1. Algorithm

            Wold = W;
            
            y = W * z;

            g = -2 * tanh(y);
             
            W = W + mu * (eye(Es1) + g * y.') * W;
                 
            % Dot product 
            
            d1 =  dot(Wold(:,1),W(:,1));
            d2 =  dot(Wold(:,2),W(:,2));

            d1V(n) = d1;
            d2V(n) = d2;
            
            % Error function
            
            sumE(1) =  (Aref(1,1) - W(1,1));
            sumE(2) =  (Aref(2,1) - W(2,1));
            sumE(3) =  (Aref(1,2) - W(1,2));
            sumE(4) =  (Aref(2,2) - W(2,2));

            se = sum(sumE);
            seV(n) = se;
    
            % Matrix condition
            try
                P = W*A;
                p(n) = cond(P); %#ok<AGROW>
                
                % If diverged, restart with new step size
            catch
                
                n = NmaxIterations;
            end
            
            % yi and gj(yj) Uncorralated 
            
            try
                
            I = mean(g*y.',4);
            i = cond(abs(I));
            pi(n) = i;
            catch
                
                n = NmaxIterations;
            end
              % Convergence criteria 
              
              %if n ~= NmaxIterations && abs(i-1) < 0.01;
              if n ~= NmaxIterations && n > 100 && abs(p(n)/p(n-100) - 1) < 0.0001 
              %if n > 50 && abs(seV(n)/seV(n-2) - 1 ) < 0.00001 || abs(seV(n)) < 0.0001 
              %tol = 0.001;
              %if n > 50 && abs(abs(d1) - 1) < tol &&  abs(abs(d2) - 1) < tol 
                    disp('Found');    
                     disp('n');
                     disp(n);
                     disp('m');
                     disp(m);
                     n = NmaxIterations;
                     m = MmaxIterations;
                
              end     
            
            n = n + 1;
            
        end
    
    m = m + 1;
end

% Information printing


figure(2);
plot(pi, 'g-')
figure(3);
plot(seV, 'b-');
figure(4)
subplot(2,1,1);
plot(d2V, 'g-');
subplot(2,1,2);
plot(d1V, 'r-');
figure(5);
plot(p, 'y-');

disp('P:')
disp(W*A);

disp('W:');
disp(W);


disp('A:');
disp(A);
disp('rank:');
disp(rank(A));

xout = W*x*10;

for n = 1:Es1;
    
    xout(n,:) = (xout(n,:) + MeanS(n));
    
end

figure(6);

subplot(2,2,1);
plot([0:size(e,2)-1],e(1,:),'-g');
hold on
subplot(2,2,2);
plot([0:size(e,2)-1],e(2,:),'-g');

subplot(2,2,3);
plot([0:size(xout,2)-1],xout(1,:),'-r');
subplot(2,2,4);
plot([0:size(xout,2)-1],xout(2,:),'-r');


figure(7);
subplot(2,2,1);
scatter(xout(1,:),x(1,:));
xlabel('1,1');
subplot(2,2,2);
scatter(xout(2,:),x(2,:));
xlabel('2,2');
subplot(2,2,3);
scatter(xout(1,:),x(2,:));
xlabel('1,2');
subplot(2,2,4);
scatter(xout(2,:),x(1,:));
xlabel('2,1');


if(abs(det(W)) > 0.01)
    
    disp('The components are independent.');
    
else
    
    disp('The components are not independent, sorry');
end
disp('Determinant:')
disp(det(W));

disp('mu')
disp(mu);
    
if play == 1
    sound(xout(1,:),Fs1);
    pause(5);
    sound(xout(2,:),Fs1);
end

end
