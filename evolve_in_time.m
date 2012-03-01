clc;
clear;
close all;
%filename = 'd:\Shared\lasagna_svn\thermal_IH_1\dump_3_3.mat';
%filename = 'd:\Shared\lasagna_svn\thermal_NH_2\dump_3_3.mat';
%filename = 'd:\Shared\lasagna_svn\output\dump_7_7.mat';
%filename = 'd:\Shared\lasagna_svn\output\dump_L_nh.mat';
%filename = 'd:\Shared\lasagna_svn\output\too_late_for_repopulation.mat';
%filename = 'd:\Shared\lasagna_svn\NH_1em4_500_3x3\dump_000_000.mat';
filename = 'd:\Shared\lasagna_svn\output\dump2.mat'

S = load(filename,'T','L','x_grid','Ps_plus','Ps_minus','Pa_plus','Pa_minus',...
    'Py_plus','Py_minus','xi','alpha_rs');

mask = S.T~=0;
last_idx = sum(mask);

follow_index = 24;
speed = 5;
start_at = 1000;%1598;
count = 1;
for i=start_at:speed:last_idx
    x_grid = sqrt(S.x_grid(:,i)-1e-12);
    xi = sqrt([S.xi(:,i),S.xi(:,i)]);
    subplot(2,3,1)
    plot(x_grid,S.Ps_plus(:,i),'LineWidth',2)
    hold on; 
        plot(xi(1,:),ylim, xi(2,:),ylim);
        plot(x_grid(follow_index),S.Ps_plus(follow_index,i),'r*');
    hold off
    title({['Temperature: ',num2str(S.T(i)*1e3,4),'MeV.'];...
        ['L = ',num2str(S.L(i))];...
        'Ps^+'})
    subplot(2,3,4)
    plot(x_grid,S.Ps_minus(:,i),'LineWidth',2)
    hold on; plot(xi(1,:),ylim, xi(2,:),ylim);hold off
    title('Ps^-')
    subplot(2,3,2)
    plot(x_grid,S.Py_plus(:,i),'LineWidth',2)
    hold on; plot(xi(1,:),ylim, xi(2,:),ylim);hold off
    title('Py^+')
    subplot(2,3,5)
    plot(x_grid,S.Py_minus(:,i),'LineWidth',2)
    hold on; plot(xi(1,:),ylim, xi(2,:),ylim);hold off
    title('Py^-')
    subplot(2,3,3)
    plot(x_grid,S.Pa_plus(:,i),'LineWidth',2)
    hold on; 
        plot(xi(1,:),ylim, xi(2,:),ylim);
    hold off
    title('Pa^+')
    subplot(2,3,6)
    plot(x_grid,S.Pa_minus(:,i),'LineWidth',2)
    hold on; 
        plot(xi(1,:),ylim, xi(2,:),ylim);
    hold off
    title('Pa^-')
    drawnow
    tt = 1;
    xvec = S.x_grid(tt:end,i);
    I = trapz(xvec,0.5./(1+exp(xvec)).*S.Ps_plus(tt:end,i).*xvec.^3);
    J = trapz(xvec,0.5./(1+exp(xvec)).*S.Pa_plus(tt:end,i).*xvec.^3);
    K = trapz(xvec,0.5./(1+exp(xvec)).*4.*ones(size(S.Pa_plus(tt:end,i))).*xvec.^3);
    M = trapz(xvec,xvec.^2./(1+exp(xvec)).*S.Pa_minus(tt:end,i));
    N = 4*trapz(xvec,xvec.^2./(exp(xvec)-1).*ones(size(S.Pa_minus(tt:end,i))));
    
    I1 = 3.046*I/(7/20*pi^4);
    I2 = I1+(J/K-1);
    Neff(count) = I1;
    Neff2(count) = I2;
    L1(count) = S.L(i);
    L2(count) = M/N; %/(8*zeta(3));
    Tvec(count) = S.T(i)*1e3;
    count = count + 1;
end
figure
plot(Tvec,Neff,Tvec,Neff2)
set(gca,'Xdir','reverse')
figure
plot(Tvec,L1,Tvec,L2)
set(gca,'Xdir','reverse')

