% This demo code reproduces the overcomplete transform learning result presented in Figure 1 of the following paper:
% 1) S. Ravishankar and Y. Bresler, "Learning overcomplete sparsifying transforms for signal processing," in Proc. IEEE ICASSP, 2013, pp. 3088-3092.

clear,clc;
close all;

image   = imread('../data/lena.jpg');

[h,w]   = size(image);
image   = double(image);

n       = 10*10;
M       = 256;
thr     = 0.6;
l2      = 4e5;  % weight of -log(det(W'*W))
l4      = l2;   % weight of incoherence
l3      = l2;   % weight of weight decay
p       = 20;   % power of incoherence
step    = 1e-7; % CG step
iter    = 130;  % iter times
cg_iter = 100;  % CG iter times
debug   = 0;
visual  = 1;
saved_path  = '../result/transform_exp/';

blocks  = my_im2col(image,[sqrt(n),sqrt(n)],1);
Y       = blocks - (ones(n,1)*mean(blocks));

Y       = whitening(Y);

W0      = 0.1*rand(M,n);

tic
[W,X]   = TransformLearning(W0,Y,iter,l2,l3,l4,p,step,cg_iter,thr,debug,visual,saved_path);
toc
Dic     = displayDictionaryElementsAsImage(W', sqrt(M), sqrt(M),sqrt(n),sqrt(n));