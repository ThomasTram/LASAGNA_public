function Do_thermalisation_cartoon
clc;clear;
close all
%filename = 'd:\Shared\lasagna_svn\te_IH_0_ndf\dump_004_008.mat'
filename = 'd:\Shared\lasagna_svn\dump_004_008.mat'


dump_at_x = [0.1,1,2,3,5,10];

load(filename,'T','Ps_plus','Pa_plus','Tres_vres','delta_m2_theta_zero','x_grid',...
    'xmin_xext_xmax','V0','V1');

%Plot options
global axfont txfont axangle txangle axesfs legfs labelfs titlefs renderer clever_interp
axfont = 'Times';%'Helvetica';%'Courier';%'AvantGarde';%'Helvetica';%'cmr10';
txfont = 'Times';
axangle = 'italic';
txangle = 'normal';
axesfs = 18;
labelfs = 22;
legfs = 10;
titlefs = 22;
renderer = 'painters';


max_x = 15;
max_y = 140;
arrow_pos = 75;
arrow_siz = max_x/7.5;

maskT = T~=0;
lastidx = sum(maskT);
T = T(maskT)*1e3;
c2theta = cos(2*delta_m2_theta_zero(2));
dm2 = delta_m2_theta_zero(1)*1e18;
xmax = xmin_xext_xmax(3);

x0 = 1.812*10^4./T.^3*sqrt(c2theta*abs(dm2)/2);
x0p = sqrt(abs(V0(maskT)./V1(maskT)));
%plot(T,x0,T,x0p)
%set(gca,'xdir','reverse')




for i=1:length(dump_at_x)
    idx=find(x0>dump_at_x(i),1);
    if isempty(idx)
        break
    end
    
    if abs(x0(idx)-dump_at_x(i))> abs(x0(idx-1)-dump_at_x(i))
        idx_x(i) = idx-1;
    else
        idx_x(i) = idx;
    end
end


for i=1:length(idx_x)
    close all
    idx = idx_x(i);
    my_figure;
    
    plot(x_grid(:,idx),25*Pa_plus(:,idx),...
         x_grid(:,idx),25*Ps_plus(:,idx),...
        'LineWidth',2)
    axis([0 max_x 0 max_y])
    set_gca;
    my_xlabel('x=p/T')
    my_ylabel('degree of thermalisation [%]')
    
    hold on
    line([x0(idx) x0(idx)],[0 max_y],'Color','k','LineStyle','--','LineWidth',2)
    [arr_x, arr_y] = ax2fig([x0(idx), x0(idx)+arrow_siz, arrow_pos, arrow_pos]);
    annotation('arrow',arr_x,arr_y,...
        'LineWidth',2);
    text(x0(idx)+0.15, arrow_pos-10,...
        {'Resonance',['@',num2str(T(idx),3),'MeV']},...
        'FontSize',axesfs,'FontName',txfont,'FontAngle',txangle,...
        'FontWeight','bold');
    hold off
    legend({'Active Species','Sterile Species'},...
        'FontName',axfont,'FontWeight','normal','FontSize',axesfs);
    
    drawnow
    save_eps_and_pdf('resonant_thermalisation_',...
        'd:\Shared\lasagna_svn\plots\cartoon',...
        [num2str(i),'.eps'])
    end

function [posfig_x, posfig_y] = ax2fig(posax)
axfig = get(gca,'Position');
axax = axis;

alpha_x = axfig(3)/(axax(2)-axax(1));
beta_x = axfig(1)-alpha_x*axax(1);
alpha_y = axfig(4)/(axax(4)-axax(3));
beta_y = axfig(2)-alpha_y*axax(3);

posfig_x = alpha_x*posax(1:2)+beta_x;
posfig_y = alpha_y*posax(3:4)+beta_y;

function my_figure
global renderer
figure('Renderer',renderer,'PaperType','A4',...
    'PaperOrientation','landscape',...
    'Color',[1 1 1],'PaperPositionMode','auto');

function my_xlabel(xtit)
global labelfs txfont txangle
% Create xlabel
xlabel(xtit,'FontWeight','bold','FontSize',labelfs,...
    'FontName',txfont,'FontAngle',txangle);

function my_ylabel(ytit)
global labelfs txfont txangle
% Create xlabel
ylabel(ytit,'FontWeight','bold','FontSize',labelfs,...
    'FontName',txfont,'FontAngle',txangle);

function my_title(tit)
global titlefs txfont txangle
title(tit,...
    'FontWeight','bold',...
    'FontSize',titlefs,...
    'FontName',txfont,...
    'FontAngle',txangle);

function set_gca
global axfont axesfs axangle
set(gca,'FontName',axfont,'FontWeight','bold','FontSize',axesfs,...
    'FontAngle',axangle,'Ydir','normal','layer','top')

function save_eps_and_pdf(prefix_filename,imagepath,fname)
filename = [imagepath,'\',prefix_filename,fname];
saveas(gcf,filename,'epsc2');
useflateencode(filename);
eps2pdf(filename,...
    'C:\Program Files\gs\gs9.01\bin\gswin32c.exe',0);

function my_colorbar
global axfont axesfs
cmap = cmrmap(256,5);
colormap(cmap)
colorbar('FontName',axfont,'FontWeight','bold','FontSize',axesfs)

function my_legend(leg_cell,varargin)
global axfont legfs
lh = legend(leg_cell,'FontName',axfont,'FontWeight','normal','FontSize',legfs);
for i=1:2:(nargin-1)
    set(lh,varargin{i},varargin{i+1})
end