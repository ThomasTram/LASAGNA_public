clear; clc;
close all
%dirname = 'd:\Shared\lasagna_radau\output\run1\';
%dirname = 'd:\Shared\lasagna_radau\output\';
dirname = '/users/ire/Desktop/RADAUS5/version_0.4/lasagnaLcheck/output/zero_nh/';


S = dir([dirname,'*.mat']);
for i=1:length(S)
    disp([num2str(i),': ',S(i).name])
end
mask_S = length(S)-3:length(S);
%S = S(mask_S);
plotstring = 'plot(';

plotstring2 = 'plot(';
count = 1;
deltaNeff_cell =cell(length(S),1);
for i=1:length(S)
    dat = load([dirname,S(i).name],'T','L','x_grid','Ps_plus',...
        'Tres_vres','alpha_rs','delta_m2_theta_zero','xi');
    mask = dat.T~=0;
    T{i} = 1e3*dat.T(mask);
    L{i} = dat.L(mask);
    xi_cell{i}=dat.xi(:,mask);
    idx_final = max(find(mask));
    if isempty(idx_final)
        T{i} = 40;
        idx_final = 1;
    end
    delta_m2(i,1)=dat.delta_m2_theta_zero(1);
    sinsq_theta(i,1)=dat.delta_m2_theta_zero(2);
    
    dvec = zeros(1,idx_final);
    for j=1:idx_final
        xvec = dat.x_grid(:,j);
        Ps_plus_vec = dat.Ps_plus(:,j);
        I = trapz(xvec,0.5./(1+exp(xvec)).*Ps_plus_vec.*xvec.^3);
        I = 3.046*I/(7/20*pi^4);
        dvec(j) = I;
    end
    deltaNeff_cell{i} = dvec;
    deltaNeff_final(i)=dvec(end);
    loglog(xvec,Ps_plus_vec)
    hold on
    if any(i==mask_S)
    plotstring = [plotstring,'T{',num2str(count),'},deltaNeff_cell{',num2str(count),'},'];
    plotstring2 = [plotstring2,'T{',num2str(count),'},xi_cell{',num2str(count),'}(1,:),'...
        'T{',num2str(count),'},xi_cell{',num2str(count),'}(2,:),'];
    legendcell{count} = ['(sin^2(2\theta), \delta m^2) = (',...
        num2str(sinsq_theta(count,1),3),',',num2str(delta_m2(count,1)*1.e18,3),')'];
    legendcell2{2*count-1} = legendcell{count};
    legendcell2{2*count} = legendcell{count};
    count = count +1;
    end
end
figure


deltaNeff=reshape(deltaNeff_final,4,4);
logsinsq_theta=reshape(log10(sinsq_theta),4,4);
logdelta_m2=reshape(sign(delta_m2).*log10(abs(delta_m2)*1.e18),4,4);
contourf(logsinsq_theta,logdelta_m2,deltaNeff)

figure
plotstring(end) = ')';
eval(plotstring)
legend(legendcell)

figure
plotstring2(end)=')';
eval(plotstring2)
legend(legendcell2)