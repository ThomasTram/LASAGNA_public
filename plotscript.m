%Generate thermalisation plots
clc;
mypath = 'd:\Shared\lasagna_svn\';

plotcell = ...
    {'te_NH_sup_radau','radau_tmp'
    'te_NH_0_radau','radau_NH_0_'
    'te_NH_1em2_radau','radau_NH_1em2_'
    'te_IH_0_radau','radau_IH_0_'
    'te_IH_1em2_radau','radau_IH_1em2_'
    'te_NH_sup_ndf','ndf_tmp'
    'te_NH_0_ndf','ndf_NH_0_'
    'te_NH_1em2_ndf','ndf_NH_1em2_'
    'te_IH_0_ndf','ndf_IH_0_'
    'te_IH_1em2_ndf','ndf_IH_1em2_'
    'tm_NH_sup_ndf','ndf_tmp'
    'tm_NH_0_ndf','ndf_NHe_0_'
    'tm_NH_1em2_ndf','ndf_NHe_1em2_'
    'tm_IH_0_ndf','ndf_IHe_0_'
    'tm_IH_1em2_ndf','ndf_IHe_1em2_'}
    loadorsavesup = [-1,0,1,0,0,-1,0,1,0,0,-1,0,1,0,0];
sizeofside = [16,32,16,32,16,16,32,16,32,16,16,32,16,32,16];

% plotcell = ...
%     {'te_NH_sup_ndf','ndf_tmp'
%         'te_NH_1em2_ndf','ndf_NH_1em2_'}
% sizeofside = [16 16];
% loadorsavesup = [-1,1];


for i=1:size(plotcell,1)
    dm_res = sizeofside(i);
    sin_res = sizeofside(i);
    lepton_number_thermalisation([mypath,plotcell{i,1},'\'],...
        plotcell{i,2},dm_res,sin_res,loadorsavesup(i));
end
