function test_linear_system
clear;
clc;
load('west0479.mat')
A = west0479;
x = [(1:479)', (-479:1:-1)'];
z = (1-1i)*x;
rng(11);
Aimag = 2*(rand(nnz(A),1)-0.5);
Az_mat = A;
Az_mat(find(Az_mat)) = Az_mat(find(Az_mat))+Aimag*1i;

b = A*x;
bz = Az_mat*z;

xsol = A\b;
zsol = Az_mat\bz;

relreal = abs(1-xsol./x);
relcmp = abs(1-zsol./z);
semilogy(1:479,relreal,1:479,relcmp)
[i,j,s] = find(Az_mat);

[~,m,~] = unique(j,'rows','last'); %Only works because of none empty columns

Ap = [0; m];
Ai = i-1;
Az = s;

savevector('testlin_Ap.dat',Ap,'int32');
savevector('testlin_Ai.dat',Ai,'int32');
savevector('testlin_Ax.dat',real(Az),'double');
savevector('testlin_Azim.dat',imag(Az),'double');
savevector('testlin_b.dat',b,'double');
savevector('testlin_bzr.dat',real(bz),'double');
savevector('testlin_bzi.dat',imag(bz),'double');

%For dense method
Adense_row_major = transpose(full(Az_mat));
Adense_row_major = [42; Adense_row_major(:)];
bdense_row_major = [42; b(:)];
bzrdense_row_major = [42; real(bz(:))];
bzidense_row_major = [42; imag(bz(:))];
savevector('testlin_Adense.dat',real(Adense_row_major),'double');
savevector('testlin_Adense_imag.dat',imag(Adense_row_major),'double');
savevector('testlin_bdense.dat',bdense_row_major,'double');
savevector('testlin_bzrdense.dat',bzrdense_row_major,'double');
savevector('testlin_bzidense.dat',bzidense_row_major,'double');


function savevector(filename,variable,precision)
fid = fopen(filename,'w');
fwrite(fid,variable,precision);
fclose(fid);
