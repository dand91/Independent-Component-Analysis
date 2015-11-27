function [PCM,U,Mout] =  PCAPICeig(NbrEig)

load picture.mat;

clock1 = clock;

%1. Break picture.

Ms1 = 588;
Ms2 = 200;
X = 280;
Y = 420;

Min = zeros(Ms1,Ms2);

n = 1;

    for i = 1:10:X
        for j = 1:20:Y
        
            T = G(i:(i+10 - 1) , j:(j+20 - 1));
            Min(n,:) = reshape(T,1,Ms2);
            n = n + 1;
        end   
    end

%2. Subtract mean.

MeanS = zeros(1,Ms1);

    for n = 1:Ms1;

        MeanS(n) = mean(Min(n,:));
        Min(n,:) = (Min(n,:) - MeanS(n));

    end    
    
%3. Calculate covariance matrix.

%C = cov(Min);

C = (1 / (Ms1-1)) * (Min.' * Min);

%4. Find eigenvectors and eigenvalues.

[U,D] = eig(C);

 
%5. Sort eigenvectors.

[Us1,Us2] = size(U);
Diag = diag(D);
min = -sum(abs(Diag(:)));
Temp = zeros(Us1,Us2);

    for n = 1:Us2
    
        [~, index] = max(Diag);
        Diag(index) = min;
        Temp(:,n) = U(:,index);
    
    end
    
Temp = Temp(:,1:NbrEig);
U = Temp;
    
%6. Project data onto vectors.

PCM = Min * U;

%7. Inverse projection.

 PCMOut = PCM * pinv(U);

%8. Add mean.

    for n = 1:Ms1;

        PCMOut(n,:) = (PCMOut(n,:) + MeanS(n));

    end

%9. Undo break.

Mout = zeros(X,Y);
n = 1;

    for i = 1:10:X
        for j = 1:20:Y

            T = PCMOut(n,:);
            Mout(i:(i+10-1),j:(j+20-1)) = reshape(T,10,20);
            n = n + 1;
            
        end   
    end
figure(300);
imshow(Mout);

%10. Information printing. 

clock2 = clock;

[c1,c2] = size(Min);
[c3,c4] = size(U);
[c5,c6] = size(PCM);

disp('');
disp('Procent of original size:');
disp( ((c3 * c4)+(c5 * c6)) / (c1 * c2) );
disp('');
disp('Time to process:');
disp( etime(clock2, clock1) );
%imwrite(G,'pic.bmp','bmp');

end
