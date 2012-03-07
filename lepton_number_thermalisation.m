function lepton_number_thermalisation(varargin)
close all
if nargin==0
    clear; clc;
    prefix_filename = 'stardust23';
    %dirname = 'd:\Shared\lasagna_svn\thermal_NH_3\';
    %dirname = 'd:\Shared\lasagna_radau\output\run1\';
    %dirname = 'd:\Shared\lasagna_radau\output\';
    %dirname = '/users/ire/Desktop/RADAUS5/version_0.4/lasagnaLcheck/output/zero_nh/';
    %dirname = 'd:\Shared\lasagna_svn\thermal_NH_2\'
    %dirname = 'd:\Shared\lasagna_svn\thermal_NH_1616\';
    %dirname = 'd:\Shared\lasagna_svn\thermal_old\thermal_NH_8\';
    dirname  = 'd:\Shared\lasagna_svn\stardust23\'
    
    dm_res = 8;
    sin_res = 8;
    
    %Plot options interactive:
    Show_final_sterile_spectrum = false;
    Show_evolution_of_resonances = true;
    Show_evolution_of_Neff = true;
    Show_sweep_plot = true;
    use_old_calculation = false;
    use_new_calculation = true;
    save_plots = true;
    view_plots = false;
else
    prefix_filename = varargin{2};
    dirname = varargin{1};
    dm_res = varargin{3};
    sin_res = varargin{4};
    
    %Plot options non-interactive:
    Show_final_sterile_spectrum = false;
    Show_evolution_of_resonances = true;
    Show_evolution_of_Neff = true;
    Show_sweep_plot = true;
    use_old_calculation = false;
    use_new_calculation = true;
    save_plots = true;
    view_plots = false;
end

refine_mesh = 10;
%Aux plot options
imagepath = 'plots';
if ~exist(imagepath,'dir')
    mkdir(imagepath)
end
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
clever_interp = true;

S = dir([dirname,'*.mat']);
for i=1:length(S)
    disp([num2str(i),': ',S(i).name])
end
%mask_S = (length(S)-3):length(S);
mask_S = 8:8:64;
lenS = floor(sqrt(length(S)));

mask_S = floor(lenS/2):lenS:length(S);
plotstring = 'plot(';
plotstring2 = 'plot(';
plotstring3 = 'plot(';

count = 1;
deltaNeff_cell =cell(length(S),1);
for i=1:length(S)
    dat = load([dirname,S(i).name],'T','L','x_grid','Ps_plus','Pa_plus',...
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
    sinsq_theta(i,1)=sin(2*dat.delta_m2_theta_zero(2))^2;
    if (delta_m2(i,1)>0)
        hierarchy = 'NH';
    else
        hierarchy = 'IH';
    end
    
    dvec = zeros(1,idx_final);
    dvec2 = dvec;
    for j=1:idx_final
        xvec = dat.x_grid(:,j);
        Ps_plus_vec = dat.Ps_plus(:,j);
        I = trapz(xvec,1.0./(1+exp(xvec)).*Ps_plus_vec.*xvec.^3);
        J = trapz(xvec,1.0./(1+exp(xvec)).*dat.Pa_plus(:,j).*xvec.^3);
        K = trapz(xvec,1.0./(1+exp(xvec)).*4.*xvec.^3);
        dvec(j) = I/K;
        dvec2(j) = I/K+(J/K-1);
    end
    deltaNeff_cell{i} = dvec;
    deltaNeff_final(i)=dvec(end);
    deltaNeff_cell2{i} = dvec2;
    deltaNeff_final2(i)=dvec2(end);
    deltaNeff_notconverged(i) = sum(mask)<(numel(mask)-3);
    [xi_unique, xi_idx] = unique(xi_cell{i}(1,:));
    start_of_sweep(i) = max(1,spline(xi_unique,T{i}(xi_idx),0.1));
    if Show_final_sterile_spectrum
        loglog(xvec,Ps_plus_vec)
        hold on
    end
    if any(i==mask_S)
        plotstring = [plotstring,'T{',num2str(i),'},deltaNeff_cell{',num2str(i),'},'];
        plotstring3 = [plotstring3,'T{',num2str(i),'},deltaNeff_cell2{',num2str(i),'},'];
        plotstring2 = [plotstring2,'T{',num2str(i),'},xi_cell{',num2str(i),'}(1,:),'...
            'T{',num2str(i),'},xi_cell{',num2str(i),'}(2,:),'];
        legendcell{count} = ['(sin^2(2\theta), \delta m^2) = (',...
            num2str(sinsq_theta(i,1),3),',',num2str(delta_m2(i,1)*1.e18,3),')'];
        legendcell2{2*count-1} = legendcell{count};
        legendcell2{2*count} = legendcell{count};
        count = count +1;
    end
end


deltaNeff=reshape(deltaNeff_final,dm_res,sin_res)';
deltaNeff2=reshape(deltaNeff_final2,dm_res,sin_res)';
logsinsq_theta_org=reshape(log10(sinsq_theta),dm_res,sin_res)';
logdelta_m2_org=reshape(log10(abs(delta_m2)*1.e18),dm_res,sin_res)';
deltaNeff_notconverged = reshape(deltaNeff_notconverged,dm_res,sin_res)';
start_of_sweep = reshape(start_of_sweep,dm_res,sin_res)';

logdelta_m2_vec = linspace(min(min(logdelta_m2_org)),max(max(logdelta_m2_org)),refine_mesh*dm_res);
logsinsq_theta_vec = linspace(min(min(logsinsq_theta_org)),max(max(logsinsq_theta_org)),refine_mesh*sin_res);
[logsinsq_theta, logdelta_m2] = meshgrid(logsinsq_theta_vec,logdelta_m2_vec);

deltaNeff = prepare_contour(logsinsq_theta_org, logdelta_m2_org,deltaNeff,...
    deltaNeff_notconverged,logsinsq_theta, logdelta_m2);
deltaNeff2 = prepare_contour(logsinsq_theta_org, logdelta_m2_org,deltaNeff2,...
    deltaNeff_notconverged,logsinsq_theta, logdelta_m2);
start_of_sweep = prepare_contour(logsinsq_theta_org, logdelta_m2_org,start_of_sweep,...
    deltaNeff_notconverged,logsinsq_theta, logdelta_m2);

if use_old_calculation
    my_figure
    contourf(logsinsq_theta,logdelta_m2,deltaNeff)
    my_title(['\delta N_{eff}, L=',num2str(L{1}(1)), ' [',hierarchy,']'])
    my_xlabel('log_{10}(sin^2(2\theta_0))')
    my_ylabel('log_{10}(|\delta m^2|)')
    set_gca
    my_colorbar
    if save_plots
        filename = 'deltaNeff_oldcalc.eps';
        save_eps_and_pdf(prefix_filename,imagepath,filename);
    end
    if ~view_plots
        close(gcf)
    end
    if Show_evolution_of_Neff
        my_figure
        plotstring(end) = ')';
        eval(plotstring)
        set(gca,'xdir','reverse');
        my_title(['L = ',num2str(L{1}(1)), ' [',hierarchy,']'])
        my_xlabel('T (MeV)')
        my_ylabel('\delta N_{eff}')
        my_legend(legendcell,'Location','NorthWest')
        set_gca
        if save_plots
            filename = 'deltaNeff_evolution_oldcalc.eps';
            save_eps_and_pdf(prefix_filename,imagepath,filename);
        end
        if ~view_plots
            close(gcf)
        end
    end
end
if use_new_calculation
    my_figure
    contourf(logsinsq_theta,logdelta_m2,deltaNeff2)
    my_title(['\delta N_{eff}, L=',num2str(L{1}(1)),' (2)', ' [',hierarchy,']'])
    my_xlabel('log_{10}(sin^2(2\theta_0))')
    my_ylabel('log_{10}(|\delta m^2|)')
    set_gca
    my_colorbar
    if save_plots
        filename = 'deltaNeff.eps';
        save_eps_and_pdf(prefix_filename,imagepath,filename);
    end
    if ~view_plots
        close(gcf)
    end
    
    if Show_evolution_of_Neff
        my_figure
        plotstring3(end) = ')';
        eval(plotstring3)
        set(gca,'xdir','reverse');
        my_title(['L = ',num2str(L{1}(1)),' (2)', ' [',hierarchy,']'])
        my_xlabel('T (MeV)')
        my_ylabel('\delta N_{eff}')
        my_legend(legendcell,'Location','NorthWest')
        set_gca
        if save_plots
            filename = 'deltaNeff_evolution.eps';
            save_eps_and_pdf(prefix_filename,imagepath,filename);
        end
        if ~view_plots
            close(gcf)
        end
    end
end

if Show_sweep_plot
    my_figure
    contourf(logsinsq_theta,logdelta_m2,start_of_sweep)
    my_title(['Temperature in MeV when resonance passes x=0.1, L=',num2str(L{1}(1)), ' [',hierarchy,']'])
    my_xlabel('log_{10}(sin^2(2\theta_0))')
    my_ylabel('log_{10}(|\delta m^2|)')
    set_gca
    my_colorbar
    set(gca,'clim',[min(min(start_of_sweep)) max(max(start_of_sweep))])
    if save_plots
        filename = 'start_of_resonance.eps';
        save_eps_and_pdf(prefix_filename,imagepath,filename);
    end
    if ~view_plots
        close(gcf)
    end
end
    
if Show_evolution_of_resonances
    my_figure
    plotstring2(end)=')';
    eval(plotstring2)
    set(gca,'xdir','reverse');
    my_legend(legendcell2,'Location','NorthWest')
    set_gca
    if save_plots
        filename = 'resonance_evolution.eps';
        save_eps_and_pdf(prefix_filename,imagepath,filename);
    end
    if ~view_plots
        close(gcf)
    end
end


function my_figure
global renderer
figure('Renderer',renderer,'PaperType','A4',...
    'PaperOrientation','landscape',...
    'Color',[1 1 1]);

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
    'FontAngle',axangle,'Ydir','normal','layer','top','clim',[0 1])

function save_eps_and_pdf(prefix_filename,imagepath,fname)
filename = [imagepath,'\',prefix_filename,fname];
saveas(gcf,filename,'epsc2');
useflateencode(filename);
eps2pdf(filename,...
    'C:\Program Files\gs\gs9.01\bin\gswin32c.exe',0);

function my_colorbar
global axfont axesfs
colorbar('FontName',axfont,'FontWeight','bold','FontSize',axesfs)

function my_legend(leg_cell,varargin)
global axfont legfs
lh = legend(leg_cell,'FontName',axfont,'FontWeight','normal','FontSize',legfs);
for i=1:2:(nargin-1)
    set(lh,varargin{i},varargin{i+1})
end

function ZI = prepare_contour(X,Y,Z,notconverged,XI,YI)
global clever_interp
Zn = Z;
Zn(notconverged) = nan;
for i=1:size(notconverged,1)
    for j=1:size(notconverged,2)
        if notconverged(j,i)
            if clever_interp
                if (sum(~isnan(Zn(j,:)))>=2)
                    try
                        guess1 = spline(X(j,:),Zn(j,:),X(j,i));
                    catch
                        guess1 = 0;
                    end
                else
                    guess1 = nan;
                end
                
                
                if (sum(~isnan(Zn(:,j)))>=2)
                    
                    try
                        guess2 = spline(Y(:,j),Zn(:,j),Y(j,i));
                    catch
                        guess2 = 0;
                    end
                else
                    guess2 = nan;
                end
                if (~isnan(guess1))&&(~isnan(guess2))
                    Z(j,i) = 0.5*(guess1+guess2);
                elseif (~isnan(guess1))
                    Z(j,i) = guess1;
                elseif (~isnan(guess2))
                    Z(j,i) = guess2;
                else
                    Z(j,i) = nan;
                end
            else
                Z(j,i) = nan;
            end
        end
    end
end
try
    sizy1 = size(Y,1);
    sizy2 = size(Y,2);
    for i=1:sizy1
        for j=1:sizy2
            if (Y(i,j) ~= Y(i,floor(sizy2/2)))
                Y(i,j) = Y(i,floor(sizy2/2));
                disp('Y corrected!')
            end
            if (X(i,j) ~= X(floor(sizy1/2),j))
                X(i,j) = X(floor(sizy1/2),j);
                disp('X corrected!')
            end
        end
    end
    ZI = interp2(X,Y,Z, XI, YI);
    disp('success!');
catch
    ZI = 0;
end
