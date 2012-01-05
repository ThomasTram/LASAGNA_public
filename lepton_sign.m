clear; clc;
close all
dirname = 'd:\Shared\chaos_dat\loop_line\';
S = dir([dirname,'*.mat']);

for i=1:length(S)
    dat=load([dirname,S(i).name],'T','L','Tres_vres','delta_m2_theta_zero',...
        'alpha_rs');
    mask = dat.T~=0;
    T{i} = 1e3*dat.T(mask);
    L{i} = dat.L(mask);
    final_L(i) = L{i}(end);
    delta_m2(i) = dat.delta_m2_theta_zero(1)*1e18;    
    sinsq2theta(i) = sin(2*dat.delta_m2_theta_zero(2))^2;
    mylegend{i} = ['vres: ',num2str(dat.Tres_vres(2)),', \alpha = ',...
        num2str(dat.alpha_rs(1)),'.'];
end
[sinsq2theta, idx] = sort(sinsq2theta);
final_L = final_L(idx);
plot(log10(sinsq2theta),final_L)
%pcolor(log10(-delta_m2),log10(sinsq2theta),final_L);