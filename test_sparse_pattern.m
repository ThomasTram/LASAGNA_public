clear;clc;
load('west0479.mat')
A = west0479;

A(find(A)) = 1;
ML_C = A + A';
ML_C(find(ML_C)) = 1;

I = load('I.dat');
J = load('J.dat');

C = sparse(I,J,ones(size(I)));
nnz(C-C')
