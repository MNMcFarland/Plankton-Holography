function A = holosub_nc(fnm,bkg,num)

% compute and subtract time averaged background
% when called with fnm argument: writes background subtracted images to /sub1 directory
% when called with empty fnm argument: returns average image
%
% USAGE: A = holosub1(fnm,bkg,num)
%
% EXAMPLES:
%   A = holosub1([],bkg,num); compute background from files in bkg and return average A
%   holosub1(fnm); compute background, subtract from each file in fnm, and save to /sub1 directory
%   holosub1(fnm,bkg,num); specify files to use to generate background
%
% INPUT:
%   fnm - (optional) name(s) of files to background subtract (cell or dir structure)
%   bkg - (optional) names of files used to compute background (cell or dir structure)
%   num - (optional) maximum number of frames to average (default=50 if bkg not specified)
%
% OUTPUT:
%   with fnm input argument: background subtracted images written to subfolder 'sub'
%   otherwise:
%   A - average image

if isstruct(fnm) && isfield(fnm,'name')
    fnm = {fnm.name};
elseif ~iscell(fnm) && ~isempty(fnm)
    fnm = {fnm}; 
end

if ~exist('bkg','var') || isempty(bkg)
    bkg = fnm;
    num = 50;
elseif isstruct(bkg) && isfield(bkg,'name')
    bkg = {bkg.name};
elseif ~iscell(bkg)
    bkg = {bkg}; 
end

if ~exist('num','var') %|| num>length(bkg)
    num = max(length(bkg),50); 
end

% set up a progress bar with cancel button
wb = waitbar(0,'computing background','CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
    'name','processing holograms');
hnd = findobj(wb,'type','axes');
set(get(hnd,'title'),'interpreter','none')
setappdata(wb,'canceling',0)
wbcntr = 0; % counter for progress bar
tot = num; % total # of images to process

%% compute average background
nci = ncinfo(bkg{1});
% group = nci.Groups.Name;
% nim = nci.Groups.Dimensions(1).Length;
dim = [nci.Groups.Dimensions.Length];
num1 = ceil(num/length(bkg)); % # images from each .nc file
A = zeros([dim(3) dim(2)]);
imcntr = 0;
for m = 1:length(bkg)
    nci = ncinfo(bkg{m});
    group = nci.Groups.Name;
    dim = [nci.Groups.Dimensions.Length];
    nim = nci.Groups.Dimensions(1).Length;
    if nim<num1
        num2 = nim-1;
    else
        num2 = num1;
    end
    for n = round(linspace(2,nim,num2))
        if getappdata(wb,'canceling') % check if cancelled
            delete(wb)
            return
        end
        H = ncread(bkg{m},[group '/Images'],[1 1 n],[dim(2:3) 1])';
        A = A + double(H);
        imcntr = imcntr+1;

        % update progress bar
        wbcntr = wbcntr+1;
        waitbar(wbcntr/tot,wb,['computing background - ' num2str(imcntr)]) 
    end
end
A = A/imcntr; % average image

%% 3. subtract average and write image files to sub directory
if ~isempty(fnm)
    waitbar(0,wb,'subtracting background')
    wbcntr = 0; % counter for progress bar
    tot = length(fnm); % total # of .nc files to process

    Az = A - mean(A,'all'); % zero centered average image
    for m = 1:length(fnm)
        mkdir([fnm{m}(1:end-3)]);
        if getappdata(wb,'canceling') % check if cancelled
            delete(wb)
            return
        end

        nci = ncinfo(fnm{m});
        group = nci.Groups.Name;
        dim = [nci.Groups.Dimensions.Length];
        for n = 2:dim(1)
            I = ncread(fnm{m},[group '/Images'],[1 1 n],[dim(2:3) 1])';
            I = double(I);
            I = I - Az;
            I = uint16(I);
            imwrite(I,[fnm{m}(1:end-3) '\' fnm{m}(1:end-3) '_' num2str(n,'%03d') '.tif'],'tif',...
                'Compression','none');
%             I = (I-min(I,[],'all'))/(max(I,[],'all')-min(I,[],'all'))*255;
%             I = uint8(I);
%             imwrite(I,[fnm{m}(1:end-3) '\sub\' fnm{m}(1:end-3) '_' num2str(n,'%03d') '.bmp'],'bmp');
        end
        % update progress bar
        wbcntr = wbcntr+1;
        waitbar(wbcntr/tot,wb,'subtracting background') 
    end
end
delete(wb)
