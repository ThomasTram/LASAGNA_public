function test_linear_system
clear;
clc;
load('west0479.mat')
A = west0479;
x = [(1:479)', (-479:1:-1)'];
z = (1-1i)*x;

b = A*x;
bz = A*z;

xsol = A\b;
zsol = A\bz;

zsolr = A\real(bz);
zsoli = A\imag(bz);

relreal = abs(1-xsol./x);
relcmp = abs(1-zsol./z);
relcmp2 =abs(1-(zsolr+1i*zsoli)./z);
semilogy(1:479,relreal,1:479,relcmp,1:479,relcmp2)
[i,j,s] = find(A);

[~,m,~] = unique(j,'rows','last'); %Only works because of none empty columns

Ap = [0; m];
Ai = i-1;
Ax = s;

savevector('testlin_Ap.dat',Ap,'int32');
savevector('testlin_Ai.dat',Ai,'int32');
savevector('testlin_Ax.dat',Ax,'double');
savevector('testlin_b.dat',b,'double');
savevector('testlin_bzr.dat',real(bz),'double');
savevector('testlin_bzi.dat',imag(bz),'double');

%For dense method
Adense_row_major = full(A)';
Adense_row_major = [42; Adense_row_major(:)];
bdense_row_major = [42; b(:)];
savevector('testlin_Adense.dat',Adense_row_major,'double');
savevector('testlin_bdense.dat',bdense_row_major,'double');

function savevector(filename,variable,precision)
fid = fopen(filename,'w');
fwrite(fid,variable,precision);
fclose(fid);
