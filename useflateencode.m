function useflateencode(filename)
fid = fopen(filename,'r');
newl = int32([13 10]);

insertafter{1}='systemdict /setdistillerparams known {';
insertafter{2}=['<< /AutoFilterColorImages false ',...
    '/ColorImageFilter /FlateEncode >> setdistillerparams'];
insertafter{3}='} if';


for K=1:1000
    linecell{K} = fgetl(fid);
    if strncmpi(linecell{K},'%%EndComments',12)
        break
    end
end
restoffile = fread(fid,inf);
fclose(fid);

fidw = fopen('tempeps.eps','w');
for l=1:K
    fwrite(fidw,linecell{l},'char');
    fwrite(fidw,newl,'char');
end
for l=1:numel(insertafter)
    fwrite(fidw,insertafter{l},'char');
    fwrite(fidw,newl,'char');
end
fwrite(fidw,restoffile);
fclose(fidw)
copyfile('tempeps.eps',filename);
%delete('tempeps.eps');
