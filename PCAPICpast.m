function [PCM,U,Mout] =  PCAPICpast(NbrEig)

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

%2. Subtracting mean.

[Ms1,~] = size(Min);
MeanS = zeros(1,Ms1);

    for n = 1:Ms1;

        MeanS(n) = mean(Min(n,:));
        Min(n,:) = (Min(n,:) - MeanS(n));

    end

%4. Calculate eigenvectors.

U = PAST(Min,100,1);

U = U(:,1:NbrEig);

%5. Project data onto vectors.

PCM = Min * U;

%6. Inverse projection.


 PCMOut = PCM * pinv(U);

%7. Add mean.

    for n = 1:Ms1;

        PCMOut(n,:) = (PCMOut(n,:) + MeanS(n));

    end

%8. Undo break.

Mout = zeros(X,Y);
n = 1;

    for i = 1:10:X
        for j = 1:20:Y

            T = PCMOut(n,:);
            Mout(i:(i+10-1),j:(j+20-1)) = reshape(T,10,20);
            n = n + 1;
        end   
    end
    figure(200);
    imshow(Mout);
%9. Information printing. 

clock2 = clock;

[c1,c2] = size(Min);
[c3,c4] = size(U);
[c5,c6] = size(PCM);

disp('');
disp('Procent of original size:');
disp( ((c3*c4)+(c5 * c6)) / (c1 * c2) );
disp('');
disp('Time to process:');
disp( etime(clock2, clock1) );

%imwrite(G,'pic.bmp','bmp');

end
