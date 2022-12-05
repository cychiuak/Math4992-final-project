clc;
clear;
close all;
T = double(imread('miss.PNG'));
normT = norm(T(:));
Groundtruth = Unfold(T, size(T), 1);
normg = norm(Groundtruth,"fro");
dim = ndims(T)
alpha = [1/dim,1/dim,1/dim];
mu = [0.01,0.01,0.01];
[m,n] = size( T( :, :, 1) );
k = randi([1,5]);
for i = 1:m
    for j = 1:n
        if(mod(k,5) == 2 || mod(k,5) == 1 )
        T(i,j,:) = 255;
        end



        k = randi([1,5]);
    end
end
epsilon = 1e-6;
 beta = 0.005*ones(1, ndims(T));
 Omega = (T < 254);
 maxIteration = 500;
 row = 1e-4;
 L = 1e-5;
 C =  0.6;
  tic;
 [Si_results, difference_S] = SiLRTC(T,Omega,alpha,beta,maxIteration,epsilon);
 toc;
 tic;
 [Ha_results, difference_H] = HaLRTC(T,Omega,alpha,row,maxIteration,epsilon);
 toc;
 tic;
 [Fa_results,difference_F] = FaLRTC(T,Omega,alpha,mu,L,C,maxIteration,epsilon);
 toc;
 subplot(2,4,1);
imshow(uint8(T));
title('Original');
subplot(2,4,2);
imshow(double(1-Omega));
title('missing values');
subplot(2,4,3);
imshow(uint8(Si_results));
title('siLRTC');
subplot(2,4,4);
imshow(uint8(Fa_results));
title('FaLRTC');
subplot(2,4,5);
imshow(uint8(Ha_results));
title('HaLRTC');
h = subplot(2,4,8);
plot(1:length(difference_S), -log(difference_S), '-.y', 'linewidth', 1.5); hold on;
plot(1:length(difference_H), -log(difference_H), '--g', 'linewidth', 1.5); hold on;
plot(1:length(difference_F), -log(difference_F), '--b', 'linewidth', 1.5); hold on;
title('Convergence Curves');
ylabel('-log(criteria)');
xlabel('iterations');