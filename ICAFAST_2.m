function [A,W,xout] = ICAFAST_wav()

warning('off','all')
warning

nMax = 1000;

play = 1;
SampleSize = 10000;

[e(1,:), Fs1] = wavread('1.wav'); %#ok<DWVRD>
[e(2,:), Fs2] = wavread('2.wav'); %#ok<DWVRD>

%play = 1;
%[e, Fs1] = wavread('2x2_Conv_Mixture.wav'); %#ok<DWVRD>
%e = e.';

%e(1,:) = Laplace(0,1,SampleSize);
%e(2,:) = Laplace(0,1,SampleSize);

%f1 = [1 0.8 0.3];
%f2 = [1 -0.3 -0.8];

%e(1,:) = filter(f1,1,Laplace(0,2,SampleSize));
%e(2,:) = filter(f2,1,Laplace(0,2,SampleSize));

%e(1,:) = filter(f1,1,normrnd(0,2,1,SampleSize));
%e(2,:) = filter(f2,1,normrnd(0,2,1,SampleSize));


%Mix signal
A = randn(2,2); 
Aref = inv(A);
x = A*e;

%Subtract mran
meanX(1) = mean(x(1,:));
meanX(2) = mean(x(2,:));

x(1,:) = x(1,:) - meanX(1);
x(2,:) = x(2,:) - meanX(2);

%Calculate covariance matrix
c = (1/size(x,2))*(x*x');

%Calculate eigenvectors 
[U,D] = eig(c);

%Whiten data
z = U.'*x;
    
%Initiate variables

n = 1;
m = 1;
alpha = 1 ;
p(1) = 0;

%Initiate W 

W = eye(2);

while n < nMax

    % Algorithm 
    
    Wold = W;
    
    w1 = W(:,1);
    w2 = W(:,2);

    w1 = mean(z*G(w1'*z,alpha)',2) - mean(Gprim(w1'*z,alpha), 2)*w1;
    w2 = mean(z*G(w2'*z,alpha)',2) - mean(Gprim(w2'*z,alpha), 2)*w2;
    
    W = [w1 w2];
    
    W = inv(sqrtm(W*W.'))*W;   

    %Error function
    
    sumE(1) =  (W(1,1) - Aref(1,1))^2;
    sumE(2) =  (W(2,1) - Aref(2,1))^2;
    sumE(3) =  (W(1,2) - Aref(1,2))^2;
    sumE(4) =  (W(2,2) - Aref(2,2))^2;
    
    se = sum(sumE);
    seV(n) = se;

    % Dot product
    
    d1 =  dot(Wold(:,1),W(:,1));
    d2 =  dot(Wold(:,2),W(:,2));
    
    d1V(n) = d1;
    d2V(n) = d2;
    
    %Matrix condition
    
    p(n) = cond(W*inv(A));

    %Convergence criteria 
    
    %if n ~= nMax && n > 100 && abs(p(n)/p(n-100) - 1) < 0.000001 
    %if n > 50 && abs(seV(n)/seV(n-2) - 1 ) < 0.00001 || abs(seV(n)) < 0.0001 
    tol = 0.00001;
    if n > 50 && abs(abs(d1) - 1) < tol &&  abs(abs(d2) - 1) < tol 
    
        
    disp('Found');    
    disp('n:');
    disp(n);
    n = nMax;

    end      
        n = n + 1;
end

%Information printing 

figure(1);
plot(seV);
figure(2)
subplot(2,1,1);
plot(d2V, 'g-');
subplot(2,1,2);
plot(d1V, 'r-');
figure(5);
plot(p);
figure(3);
disp('W:');
disp(W);

if(abs(det(W)) > 0.01)
    
    disp('The components are independent.');
    
else
    
    disp('The components are not independent, sorry');
end

disp('Determinant:')
disp(det(W));

xout = W*x;

xout(1,:) = xout(1,:) + meanX(1);
xout(2,:) = xout(2,:) + meanX(2);

subplot(2,2,1);
plot(e(1,:),'b');
subplot(2,2,2);
plot(e(2,:),'b');

subplot(2,2,3);
plot(xout(1,:),'r');
subplot(2,2,4);
plot(xout(2,:),'r');


figure(10);
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

if play == 1
    
    sound(xout(1,:),Fs1);
    pause(8);
    sound(xout(2,:),Fs1);
    
end

end

function [out1] = G(arg1,a1)

    out1 = tanh(a1 * arg1);

end

function [out2] = Gprim(arg2,a2)

    out2 = a2 - a2 * tanh(a2 * arg2).^2;

end
