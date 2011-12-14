clear; clc;
close all
%dirname = 'd:\Shared\lasagna_radau\output\run1\';
%dirname = 'd:\Shared\lasagna_radau\output\';
%dirname = 'd:\Shared\lasagna_new\output\';
dirname = 'd:\Shared\lasagna_svn\';

S = dir([dirname,'*.mat']);
for i=1:length(S)
    disp([num2str(i),': ',S(i).name])
end
%mask = [1,5,6,7,9,10];
%S = S(mask);
for i=1:length(S)
    disp([num2str(i),': ',S(i).name])
end

for i=1:length(S)
    dat = load([dirname,S(i).name],'T','L','Tres_vres','alpha_rs');
    mask = dat.T~=0;
    T{i} = 1e3*dat.T(mask);
    L{i} = dat.L(mask);
    L_mod{i} = sign(L{i}).*abs(L{i}).^(1/10);
    mylegend{i} = ['vres: ',num2str(dat.Tres_vres(2)),', \alpha = ',...
        num2str(dat.alpha_rs(1)),'.'];
end
plotstring = 'T{1},L{1}';
for i=2:length(S)
    plotstring = [plotstring,',T{',num2str(i),'},L{',num2str(i),'}'];
end
eval(['plot(',plotstring,')'])
set(gca,'xdir','reverse')
legend(mylegend,'Location','southwest')
xlabel('T (MeV)')
ylabel('L')

figure
plotstring = 'T{1},abs(L{1})';
for i=2:length(S)
    plotstring = [plotstring,',T{',num2str(i),'},abs(L{',num2str(i),'})'];
end
eval(['semilogy(',plotstring,')'])
set(gca,'xdir','reverse')
legend(mylegend,'Location','southwest')
xlabel('T (MeV)')
ylabel('L')

figure

plotstring = 'T{1},L_mod{1}';
for i=2:length(S)
    plotstring = [plotstring,',T{',num2str(i),'},L_mod{',num2str(i),'}'];
end
eval(['plot(',plotstring,')'])
set(gca,'xdir','reverse')
legend(mylegend,'Location','southwest')
xlabel('T (MeV)')
ylabel('L^{1/10}')
