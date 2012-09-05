clc;
clear;
close all;
%filename = 'd:\Shared\lasagna_svn\thermal_IH_1\dump_3_3.mat';
%filename = 'd:\Shared\lasagna_svn\thermal_NH_2\dump_3_3.mat';
%filename = 'd:\Shared\lasagna_svn\output\dump_7_7.mat';
%filename = 'd:\Shared\lasagna_svn\output\dump_L_nh.mat';
%filename = 'd:\Shared\lasagna_svn\output\too_late_for_repopulation.mat';
%filename = 'd:\Shared\lasagna_svn\NH_1em4_500_3x3\dump_000_000.mat';
%%filename = 'd:\Shared\lasagna_svn\trig_1e5_NH\dump_006_001.mat'
%filename = 'd:\Shared\lasagna_svn\output\dump2.mat'

%filename = 'd:\Shared\lasagna_svn\te_NH_1em2_ndf\dump_015_015.mat'
%filename = 'd:\Shared\lasagna_svn\te_NH_sup_ndf\dump_015_015.mat'
%filename = 'd:\Shared\lasagna_svn\te_NH_1em2_ndf\dump_015_015.mat'
%filename = 'd:\Shared\lasagna_svn\te_IH_0_ndf\dump_004_024.mat'
%filename = 'd:\Shared\lasagna_svn\te_IH_0_ndf\dump_010_024.mat'
%filename = 'd:\Shared\lasagna_svn\te_IH_0_ndf\dump_013_024.mat'
%filename = 'd:\Shared\lasagna_svn\te_IH_0_ndf\dump_004_016.mat'

%filename = 'd:\Shared\lasagna_svn\te_IH_0_ndf\dump_004_008.mat'
%filename = 'd:\Shared\lasagna_svn\te_IH_0_ndf\dump_008_004.mat'

%filename = 'd:\Shared\chaos\dump_004_008.mat'
%filename = 'd:\Shared\chaos\test_004_008.mat'
filename = 'output/no_chaos2_rs0p1_Ti_14p5.mat';
filename = 'output/dump.mat';

S = load(filename,'T_vec','L_vec','x_grid','Ps_plus','Ps_minus','Pa_plus','Pa_minus',...
    'Py_plus','Py_minus','xi_vec','alpha_rs','delta_m2_theta_zero','xi_vec','v_grid');

mask = S.T_vec~=0;
last_idx = sum(mask);

follow_index = 24;
speed = 1;
start_at = 1;%1598;
count = 1;
xlimits = [0 10];

sqrtstuff = sqrt(cos(2*S.delta_m2_theta_zero(2))*abs(S.delta_m2_theta_zero(1))*1e18/2);
for i=start_at:speed:last_idx
    x0 = 1.812e4*(S.T_vec(i)*1e3)^(-3)*sqrtstuff;
    
    x_grid = sqrt(S.x_grid(:,i)-1e-12);
    %x_grid = S.v_grid(:,i);
    %xi = sqrt([S.xi(:,i),S.xi(:,i)]);
    xi = sqrt([x0,x0;x0,x0]);
    subplot(3,3,1)
    plot(x_grid,S.Ps_plus(:,i),'LineWidth',2)
    xlim(xlimits);
    hold on; 
        plot(xi(1,:),ylim, xi(2,:),ylim);
        plot(x_grid(follow_index),S.Ps_plus(follow_index,i),'r*');
    hold off
    title({['Temperature: ',num2str(S.T_vec(i)*1e3,4),'MeV.'];...
        ['L = ',num2str(S.L_vec(i))];...
        'Ps^+'})
    subplot(3,3,4)
    plot(x_grid,S.Ps_minus(:,i),'LineWidth',2)
    xlim(xlimits);
    hold on; plot(xi(1,:),ylim, xi(2,:),ylim);hold off
    title('Ps^-')
    subplot(3,3,2)
    plot(x_grid,S.Py_plus(:,i),'LineWidth',2)
    xlim(xlimits);
    hold on; plot(xi(1,:),ylim, xi(2,:),ylim);hold off
    title('Py^+')
    subplot(3,3,5)
    plot(x_grid,S.Py_minus(:,i),'LineWidth',2)
    xlim(xlimits);
    hold on; plot(xi(1,:),ylim, xi(2,:),ylim);hold off
    title('Py^-')
    subplot(3,3,3)
    plot(x_grid,S.Pa_plus(:,i),'LineWidth',2)
    xlim(xlimits);
    hold on; 
        plot(xi(1,:),ylim, xi(2,:),ylim);
    hold off
    title('Pa^+')
    subplot(3,3,6)
    plot(x_grid,S.Pa_minus(:,i),'LineWidth',2)
    xlim(xlimits);
    hold on; 
        plot(xi(1,:),ylim, xi(2,:),ylim);
    hold off
    title('Pa^-')
    subplot(3,3,7)
    plot(x_grid,0.5*(S.Pa_plus(:,i)+S.Pa_minus(:,i)),'LineWidth',2)
    xlim(xlimits);
    hold on; 
        plot(xi(1,:),ylim, xi(2,:),ylim);
    hold off
    ylim([-1,5])
    title('Pa')
    %
    subplot(3,3,8)
    plot(x_grid,0.5*(S.Pa_plus(:,i)-S.Pa_minus(:,i)),'LineWidth',2)
    hold on; 
        plot(xi(1,:),ylim, xi(2,:),ylim);
    hold off
    xlim(xlimits);
    ylim([-1,5])
    title('Pa bar')
    %
    subplot(3,3,9)
    plot(x_grid,0.5*(S.Ps_minus(:,i)+S.Ps_plus(:,i)),'LineWidth',2)
    xlim(xlimits);
    hold on; plot(xi(1,:),ylim, xi(2,:),ylim);hold off
    xlim(xlimits);
    ylim([-1,5])
    title('Ps')
   
    %
    drawnow
    tt = 1;
    xvec = S.x_grid(tt:end,i);
    I = trapz(xvec,0.5./(1+exp(xvec)).*S.Ps_plus(tt:end,i).*xvec.^3);
    J = trapz(xvec,0.5./(1+exp(xvec)).*S.Pa_plus(tt:end,i).*xvec.^3);
    K = trapz(xvec,0.5./(1+exp(xvec)).*4.*ones(size(S.Pa_plus(tt:end,i))).*xvec.^3);
    M = trapz(xvec,xvec.^2./(1+exp(xvec)).*S.Pa_minus(tt:end,i));
    N = 4*trapz(xvec,xvec.^2./(exp(xvec)-1).*ones(size(S.Pa_minus(tt:end,i))));
    rho_ss = 0.25./(exp(xvec)+1).*(S.Ps_plus(tt:end,i)+S.Ps_minus(tt:end,i));
    rho_ss_bar = 0.25./(exp(xvec)+1).*(S.Ps_plus(tt:end,i)-S.Ps_minus(tt:end,i));
    Q = trapz(xvec,xvec.^2.*rho_ss);
    P = trapz(xvec,xvec.^2.*rho_ss_bar);
    %R = trapz(xvec,1.0./(1+exp(xvec)).*ones(size(S.Pa_plus(tt:end,i))).*xvec.^2);
    R = 0.25*N;
    
    I1 = I/K;
    %I2 = I1+(J/K-1);
    I2 = I/K+J/K-1;
    n_s(count) = Q/R;
    n_s_bar(count) = P/R;
    Neff(count) = I1;
    Neff2(count) = I2;
    L1(count) = S.L_vec(i);
    L2(count) = M/N; %/(8*zeta(3));
    Tvec(count) = S.T_vec(i)*1e3;
    count = count + 1;
    %pause
end

figure
plot(Tvec,Neff,Tvec,Neff2)
set(gca,'Xdir','reverse')
figure
plot(Tvec,L1,Tvec,L2)
set(gca,'Xdir','reverse')
figure
plot(Tvec,n_s,Tvec,n_s_bar)
set(gca,'Xdir','reverse')

