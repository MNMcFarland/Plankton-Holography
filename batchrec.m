% do reconstruction in batches to avoid slow-down over time (due to memory fragmentation?)
% requires the following variables stored in C:\Users\mmcfarland\Documents\Malcbox\On\holography\temp.mat
%
% curdir - current working directory string (full path)
% recinc - reconstruction depth interval in um
% recdpt - total reconstruction depth in um
% res - resolution in um/pixel
% wvl - laser wavelength in um
% fnm - hologram file names (cell array)
% ind - indices of filenames in fnm to reconstruct
% lmag - (optional) image processing flag
% t - (optional) threshold for image processing

if ~exist('lmag','var') || isempty(lmag); lmag=1; end
if ~exist('t','var'); t=[]; end

load C:\Users\mmcfarland\Documents\Malcbox\On\holography\temp.mat;

% curdir=pwd;
if ~strcmp(curdir,pwd); cd(curdir); end
% if ~strcmp(curdir(end-2:end),'sub'); cd sub; curdir=pwd; end
% if ~strcmp(curdir,['E:\ES2015\Station ' num2str(st) '\Imperx\sub']);
%     cd(['E:\ES2015\Station ' num2str(st) '\Imperx\sub']); 
% end

if ~exist('strt','var'); strt=1; end
if ~exist('tm','var'); tm=[]; end
incr=50; % increment, # of holograms to process in each batch
stp=min(length(ind),strt+incr-1); % stop index
disp(tm)
tic;
% recsplatgpu(fnm(ind(strt:stp)),0:recinc:recdpt,res,wvl,lmag,t);
recsplatgpu2(fnm(ind(strt:stp)),0:recinc:recdpt,res,wvl,[0.5 1024],16);
tm1=toc;
tm=[tm; tm1];
strt=strt+incr; % increment start index
% save('C:\Users\mmcfarland\Documents\Malcbox\On\ES2015\temp.mat','strt','incr','tm','-append')
save('C:\Users\mmcfarland\Documents\Malcbox\On\holography\temp.mat','strt','incr','tm','-append')
if stp<length(ind)
%     system('mtlbrestart')
%     system('C:\Users\mmcfarland\Documents\Malcbox\On\ES2015\mtlbrestart')
    system('C:\Users\mmcfarland\Documents\Malcbox\On\holography\mtlbrestart')
end
