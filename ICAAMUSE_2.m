function [A,W] = ICAAMUSE_wav()

play = 1;
SampleSize = 1000;
f1 = [1 0.8 0.3];
f2 = [1 -0.3 -0.8];

nbr = 2;

[e(1,:), Fs1] = wavread('1.wav'); %#ok<DWVRD>
[e(2,:), Fs2] = wavread('2.wav'); %#ok<DWVRD>

%[e, Fs1] = wavread('2x2_Conv_Mixture.wav'); %#ok<DWVRD>
%e = e.';

%e(1,:) = filter(f1,1,Laplace(0,2,SampleSize));
%e(2,:) = filter(f2,1,Laplace(0,2,SampleSize));

%e(1,:) = filter(f1,1,normrnd(0,2,1,SampleSize));
%e(2,:) = filter(f2,1,normrnd(0,2,1,SampleSize));


% Mix signal

A = rand(nbr);
x = A * e;

%1. Center data.

[Es1,Es2] = size(e);
MeanS = zeros(1,Es1);

for n = 1:Es1;
    
    MeanS(n) = mean(e(n,:));
    e(n,:) = (e(n,:) - MeanS(n));
    
end


% Calculate covariance matrix.

C = (1 / (Es2-1)) * (x * x.');

% Find eigenvectors and eigenvalues.

[D,U] = eig(C);

% Whiten data

zW = D.'*x;

% Initiate variables
tau = 2;

z(1,:) = zW(1,1:end-tau);
z(2,:) = zW(2,1:end-tau);

ztau(1,:) = zW(1,tau + 1:end);
ztau(2,:) = zW(2,tau + 1:end);

%5. Run algorithm

    Ctau = mean(z*ztau.',4);
    Chat = (1/2) * (Ctau + Ctau.');
    [W,z] = eig(Chat);

% Information printing 

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

for n = 1:2;
    
    xout(n,:) = (xout(n,:) + MeanS(n));
    
end

disp('A^-1:');
disp(inv(A));
disp('rank:');
disp(rank(A));

subplot(2,3,1);
plot(0:size(e,2)-1,e(1,:),'-g');
hold on
subplot(2,3,2);
plot(0:size(e,2)-1,e(2,:),'-g');

subplot(2,3,4);
plot(0:size(x,2)-1,xout(1,:),'-r');
subplot(2,3,5);
plot(0:size(x,2)-1,xout(2,:),'-r');

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
