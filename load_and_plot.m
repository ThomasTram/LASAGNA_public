function S = ouput_matrix
%Read output matrix:
clear;clc;
close all;
%filename = 'output/sim_alpha_1.mat'
%momentum bin of special interest: in [1; vres]
mbin = 51;
analytic_static = true;
Ti_lowtemp = 6;%MeV
filename = 'output/dump.mat';

load(filename)
S = load(filename);
%Use T in MeV:
T = T*1000;
last_idx = length(T)-sum(T==0);

%Construct linear combinations of (more) physical quantities:
P0_plus = 0.5*(Pa_plus+Ps_plus);
P0_minus = 0.5*(Pa_minus+Ps_minus);
Pz_plus = 0.5*(Pa_plus-Ps_plus);
Pz_minus = 0.5*(Pa_minus-Ps_minus);

P0 = 0.5*(P0_plus+P0_minus);
P0_bar = 0.5*(P0_plus-P0_minus);
Pz = 0.5*(Pz_plus+Pz_minus);
Pz_bar = 0.5*(Pz_plus-Pz_minus);
Px = 0.5*(Px_plus+Px_minus);
Px_bar = 0.5*(Px_plus-Px_minus);
Py = 0.5*(Py_plus+Py_minus);
Py_bar = 0.5*(Py_plus-Py_minus);
%clear Pa_plus Pa_minus P0_plus P0_minus Pz_plus Pz_minus Pz_plus Pz_minus Py_plus P0_bar Px_plus Px_minus
%Auxillary stuff:
if is_electron
    C_alpha = 1.27;
else
    C_alpha = 0.92;
end
mask = ones(1,Tres_vres(2));

Vx_grid = Vx(mask,:)./x_grid;
V0_grid = V0(mask,:)./x_grid;
V1_grid = V1(mask,:).*x_grid;
Vz_grid = V0_grid+V1_grid;%+VL(mask,:);
Vz_bar_grid = V0_grid+V1_grid;%-VL(mask,:);

Gamma_grid = C_alpha*(1.16637e-5)^2*x_grid.*(1e-3*T(mask,:)).^5;
D_grid = 0.5*Gamma_grid;

%Static approximation:
if (analytic_static)
    Px_static_grid = Pz.*Vx_grid.*Vz_grid./(D_grid.^2+Vz_grid.^2);
    Px_bar_static_grid = Pz_bar.*Vx_grid.*Vz_bar_grid./(D_grid.^2+Vz_bar_grid.^2);
    Py_static_grid = Pz.*(-Vx_grid.*D_grid)./(D_grid.^2+Vz_grid.^2);
    Py_bar_static_grid = Pz_bar.*(-Vx_grid.*D_grid)./(D_grid.^2+Vz_bar_grid.^2);
else
    Px_static_grid = Vx_grid.*Vz_grid./(D_grid.^2+Vz_grid.^2);
    Px_bar_static_grid = Vx_grid.*Vz_bar_grid./(D_grid.^2+Vz_bar_grid.^2);
    Py_static_grid = (-Vx_grid.*D_grid)./(D_grid.^2+Vz_grid.^2);
    Py_bar_static_grid = (-Vx_grid.*D_grid)./(D_grid.^2+Vz_bar_grid.^2);
end
Py_static_plus = Py_static_grid + Py_bar_static_grid;
Py_static_minus = Py_static_grid - Py_bar_static_grid;
Px_static_plus = Px_static_grid + Px_bar_static_grid;
Px_static_minus = Px_static_grid - Px_bar_static_grid;


%Low temperature approximation
idx_Ti=max(1,sum(T>Ti_lowtemp));
% mask_integral = idx_Ti:last_idx;
% HT = sqrt(4*pi^3*3.4/45)*T(mask_integral).^3/(1.22e22);
% HT_grid = HT(mask,:);
% 
% phase = trapz(T(mask_integral),-sqrt(Vx_grid(:,mask_integral).^2+...
%     Vz_grid(:,mask_integral).^2)./HT_grid,2);
% phase_bar = trapz(T(mask_integral),-sqrt(Vx_grid(:,mask_integral).^2+...
%     Vz_bar_grid(:,mask_integral).^2)./HT_grid,2);

alpha_grid = asin(Vz_grid./sqrt(Vx_grid.^2+Vz_grid.^2));
alpha_bar_grid = asin(Vz_bar_grid./sqrt(Vx_grid.^2+Vz_bar_grid.^2));

Pz_lastidx = sin(alpha_grid(:,idx_Ti)).*sin(alpha_grid(:,last_idx)).*Pz(:,idx_Ti);
Pz_bar_lastidx = sin(alpha_bar_grid(:,idx_Ti)).*sin(alpha_bar_grid(:,last_idx)).*Pz_bar(:,idx_Ti);

Pz_lastidx3 = Vz_grid(:,idx_Ti).*Vz_grid(:,last_idx)...
    ./sqrt(Vx_grid(:,last_idx).^2+Vz_grid(:,last_idx).^2)...
    ./sqrt(Vx_grid(:,idx_Ti).^2+Vz_grid(:,idx_Ti).^2)...
    .*Pz(:,idx_Ti);

Pz_lastidx2 = Pz_lastidx + cos(alpha_grid(:,idx_Ti)).^2.*Pz(:,idx_Ti);
Pz_bar_lastidx2 = Pz_bar_lastidx + cos(alpha_bar_grid(:,idx_Ti)).^2.*Pz_bar(:,idx_Ti);

%Plot stuff
scrsz = get(0,'ScreenSize');

%Debug stuff:
Py_static_mean = mean(Py_static_grid(:,1:last_idx),2);
Py_mean = mean(Py(:,1:last_idx),2);
difPy = Py_mean-Py_static_mean;


figure('OuterPosition',[1 scrsz(4)/10 scrsz(3) 0.9*scrsz(4)])
subplot_rows = 2;
subplot_cols = 3;
nice_plot(subplot_rows,subplot_cols,1,T,Py_minus(1,:),last_idx,...
    'Title','Py_{minus}');
idx = last_idx;
nice_plot(subplot_rows,subplot_cols,2,x_grid(:,idx),Py_minus(:,idx)',Tres_vres(2),...
    'Title','Py_{minus}','xdir','normal','xlabel','x = p/T');
nice_plot(subplot_rows,subplot_cols,3,T,Ps_plus(mbin,:),last_idx,...
    'Title','Ps_{plus}');
nice_plot(subplot_rows,subplot_cols,4,T,Ps_minus(mbin,:),last_idx,...
    'Title','Ps_{minus}');
nice_plot(subplot_rows,subplot_cols,5,x_grid(:,idx),Ps_plus(:,idx)',Tres_vres(2),...
    'Title','Ps_{plus}','xdir','normal','xlabel','x = p/T');
nice_plot(subplot_rows,subplot_cols,6,x_grid(:,idx),Ps_minus(:,idx)',Tres_vres(2),...
    'Title','Ps_{minus}','xdir','normal','xlabel','x = p/T');

figure('OuterPosition',[1 scrsz(4)/10 scrsz(3) 0.9*scrsz(4)])
%figure('OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)])
subplot_rows = 2;
subplot_cols = 4;
nice_plot(subplot_rows,subplot_cols,1,T,P0(mbin,:),last_idx,...
    'Title','P0');
%nice_plot(subplot_rows,subplot_cols,2,T,P0_bar(1,:),last_idx,...
%    'Title','P0_{bar}');
nice_plot(subplot_rows,subplot_cols,2,T,xi(:,:),last_idx,...
    'Title','Resonances:');
nice_plot(subplot_rows,subplot_cols,3,T,Pz(mbin,:),last_idx,...
    'Title','Pz');

nice_plot(subplot_rows,subplot_cols,4,T,Px(mbin,:),last_idx,...
    'Title','Px');
hold on;plot(T(1:last_idx),Px_static_grid(mbin,1:last_idx),'-r','LineWidth',3); hold off;

nice_plot(subplot_rows,subplot_cols,5,T,Py(mbin,:),last_idx,...
    'Title','Py');
hold on;plot(T(1:last_idx),Py_static_grid(mbin,1:last_idx),'-r','LineWidth',3); hold off;

nice_plot(subplot_rows,subplot_cols,6,T,I_conserved(:,:),last_idx,...
    'Title','Conserved quantities');

if L(last_idx)<0
    sign_string = '(-)';
else
    sign_string = '(+)';
end
nice_plot(subplot_rows,subplot_cols,7,T,abs(L),last_idx,...
    'Title',['Leptonic asymmetry L. Sign: ',sign_string],'yscale','log');

subplot(subplot_rows,subplot_cols,8);
plot(v_grid(:,last_idx),x_grid(:,last_idx));
title('x(v)');

figure('OuterPosition',[1 scrsz(4)/10 scrsz(3) 0.9*scrsz(4)]);
linecell = {'-','--','-.',':'};
colorcell = {'k','b','b','b','b'};
subplot(2,2,1)
hold on
for i=1:4
    plot(x_grid(:,i),abs(Py(:,i)),[colorcell{i},linecell{i}],...
        x_grid(:,i),abs(Py_static_grid(:,i)),['r',linecell{i}]);
    legendcell{2*i-1} = ['T = ',num2str(T(i)),'MeV'];
    legendcell{2*i} = ['Stat. approx: T = ',num2str(T(i)),'MeV'];
end
set(gca,'xscale','log');set(gca,'yscale','log');
legend(legendcell)
clear legendcell
title('P_y(x)')

subplot(2,2,2)
hold on
for i=1:4
    plot(x_grid(:,i),abs(Px(:,i)),[colorcell{i},linecell{i}],...
        x_grid(:,i),abs(Px_static_grid(:,i)),['r',linecell{i}]);
    legendcell{2*i-1} = ['T = ',num2str(T(i)),'MeV'];
    legendcell{2*i} = ['Stat. approx: T = ',num2str(T(i)),'MeV'];
end
set(gca,'xscale','log');set(gca,'yscale','log');
legend(legendcell)
clear legendcell
title('P_x(x)')

subplot(2,2,3)
hold on
for i=1:4
    plot(x_grid(:,i),abs(Pz(:,i)),[colorcell{i},linecell{i}]);
    legendcell{i} = ['T = ',num2str(T(i)),'MeV'];
end
set(gca,'xscale','log');
legend(legendcell)
title('P_z(x)')

subplot(2,2,4)
hold on
for i=1:4
    delta = Vx(:,i)./sqrt(D_grid(:,i).^2+Vz_grid(:,i).^2);
    plot(x_grid(:,i),abs(delta),[colorcell{i},linecell{i}]);
    legendcell{i} = ['T = ',num2str(T(i)),'MeV'];
end
set(gca,'xscale','log','yscale','log');
legend(legendcell)
title('$|V_x/\sqrt{D^2+V_z^2}|$','Interpreter','latex')


figure
xvec = sqrt(x_grid(:,last_idx));
subplot(2,2,1)
plot(xvec,Pz(:,last_idx),xvec,Pz(:,idx_Ti),xvec,Pz_lastidx,xvec,Pz_lastidx2,xvec,Pz_lastidx3)
title(['Pz at T=',num2str(T(last_idx)),'MeV. Low temperature evolution imposed at T=',...
    num2str(T(idx_Ti)),'MeV']);
subplot(2,2,2)
plot(xvec,Pz_bar(:,last_idx),xvec,Pz_bar_lastidx,xvec,Pz_bar_lastidx2)
title(['Pz bar at T=',num2str(T(last_idx)),'MeV.']);
%figure
%nice_plot(1,1,1,T,vi,last_idx)
subplot(2,2,3)
semilogy(xvec,abs(Vx_grid(:,last_idx)./Vz_grid(:,last_idx)),...
    xvec,abs(Vx_grid(:,idx_Ti)./Vz_grid(:,idx_Ti)),...
    xvec,abs(Vx_grid(:,last_idx)./Vz_bar_grid(:,last_idx)),...
    xvec,abs(Vx_grid(:,idx_Ti)./Vz_bar_grid(:,idx_Ti)))
title('|Vx/Vz| and |Vx/Vz_{bar}|')
legend(['At T=',num2str(T(last_idx)),'MeV'],['At T=',num2str(T(idx_Ti)),'MeV'],...
    ['At T=',num2str(T(last_idx)),'MeV'],['At T=',num2str(T(idx_Ti)),'MeV'])
subplot(2,2,4)
%semilogy(xvec,Vx_grid(:,last_idx),xvec,abs(Vz_grid(:,last_idx)),...
%    xvec,abs(Vz_bar_grid(:,last_idx)))
plot(xvec,Vx_grid(:,last_idx),xvec,Vz_grid(:,last_idx),...
    xvec,Vz_bar_grid(:,last_idx))

figure('OuterPosition',[1 scrsz(4)/10 scrsz(3) 0.9*scrsz(4)])
%figure('OuterPosition',[1 scrsz(4) scrsz(3) scrsz(4)])
subplot_rows = 2;
subplot_cols = 3;
nice_plot(subplot_rows,subplot_cols,1,T,Px_plus(mbin,:),last_idx,...
    'Title','Px_{plus}');
hold on;plot(T(1:last_idx),Px_static_plus(mbin,1:last_idx),'-r','LineWidth',2); hold off;

nice_plot(subplot_rows,subplot_cols,2,T,Px_minus(mbin,:),last_idx,...
    'Title','Px_{minus}');
hold on;plot(T(1:last_idx),Px_static_minus(mbin,1:last_idx),'-r','LineWidth',2); hold off;

nice_plot(subplot_rows,subplot_cols,4,T,Py_plus(mbin,:),last_idx,...
    'Title','Py_{plus}');
hold on;plot(T(1:last_idx),Py_static_plus(mbin,1:last_idx),'-r','LineWidth',2); hold off;

nice_plot(subplot_rows,subplot_cols,5,T,Py_minus(mbin,:),last_idx,...
    'Title','Py_{minus}');
hold on;plot(T(1:last_idx),Py_static_minus(mbin,1:last_idx),'-r','LineWidth',2); hold off;

nice_plot(subplot_rows,subplot_cols,3,T,Px(mbin,:),last_idx,...
    'Title','Px');
hold on;plot(T(1:last_idx),Px_static_grid(mbin,1:last_idx),'-r','LineWidth',2); hold off;

nice_plot(subplot_rows,subplot_cols,6,T,Py(mbin,:),last_idx,...
    'Title','Py');
hold on;plot(T(1:last_idx),Py_static_grid(mbin,1:last_idx),'-r','LineWidth',2); hold off;
function nice_plot(rows, cols, number, x, y, lastidx, varargin)
subplot(rows, cols, number)
plot(x(1:lastidx),y(:,1:lastidx))
set(gca,'xdir','reverse');
xlabel('Temperature [MeV]');
for n=1:2:(nargin-6)
    switch lower(varargin{n})
        case 'title'
            title(varargin{n+1});
        case 'xlabel'
            xlabel(varargin{n+1});
        case 'ylabel'
            ylabel(varargin{n+1});
        otherwise
        set(gca,varargin{n},varargin{n+1})
    end
end
