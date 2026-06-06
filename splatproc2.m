function [dat,L] = splatproc2(fnm,res,t,f,g,dst,crp)

% process and analyze splat images
% segments regions <= t (light background)
% with optional image cropping
% Uses regionprops.m
%
% USAGE: [dat,L] = splatproc2(fnm,res,t,f,g,dst,crp)
%
% INPUT:
%   fnm - splat image file name(s) (string, cell array, or dir stucture)
%   res - image resolution in um/pixel
%   t - threshold level, default = 10000
%   f - 6 element vector of wavelength cutoff values for frequency filters
%       [hi_hi hi_lo md_hi md_lo lo_hi lo_lo], optional
%       default = [3 10 10 30 30 150]
%   g - 3 element vector of gain levels for each frequency range
%       default = [.8 1.1 1.3]
%   dst - distance (in pixels) within which to merge regions
%       default = 10, set to 0 to prevent region merging
%   crp - (optional) crop rectangle [left top width height]
%
% OUTPUT: 
%   dat - data table for all input files
%       column -  data type
%       ---------------------------
%          1   -  Unfilled area (um^2)
%          2   -  Centroid
%          3   -  Bounding Box
%          4   -  MajorAxisLength (um)
%          5   -  MinorAxisLength (um)
%          6   -  Eccentricity, ratio of distance between ellipse foci and major axis length
%          7   -  Orientation, degrees from vertical
%          8   -  Convex Area, area of convex hull (um^2)
%          9   -  filled area (um^2)
%          10  -  Equivalent spherical diameter (um)
%          11  -  Solidity, proportion of the pixels in the convex hull that are also in the region
%          12  -  Extent, ratio of pixels in the region to pixels in the total bounding box
%          13  -  Mean intensity / distance
%          14  -  Min intensity / distance
%          15  -  Max intensity / distance
%          16  -  Image number
%   L - label matrix

%          10  -  Euler number, number of objects in the region minus the number of holes in those objects
%          14  -  Perimeter length (um), not valid for discontiguous regions

% Malcolm McFarland
% mmcfarland@fau.edu

if isnumeric(fnm) && all(size(fnm)>1) 
    S = fnm;
    fnm = {'0000000000000000000000'};
elseif isstruct(fnm) && isfield(fnm,'name')
    fnm={fnm.name};
    S=imread(fnm{1});
elseif ~iscell(fnm)
    fnm={fnm};
    S=imread(fnm{1});
end
if ~exist('res','var') || isempty(res); res=1; end
% if ~exist('t','var') || isempty(t); t=10000; end
if ~exist('f','var') || isempty(f); f=[3 10 10 30 30 150]; end
if ~exist('g','var') || isempty(g); g=[.8 1.1 1.3]; end
if ~exist('dst','var') || isempty(dst); dst=10; end
% t=gpuArray(t);

% S=imread(fnm{1});
[Ny,Nx]=size(S);
if ~exist('t','var') || isempty(t); t=0.15259*double(intmax(class(S))); end % t=10000 for uint16 or 39 for uint8
if ~exist('crp','var') || isempty(crp); crp=[1 1 Nx Ny]; end
S=S(crp(2):crp(2)+crp(4)-1,crp(1):crp(1)+crp(3)-1);
[Ny,Nx]=size(S);

%% generate frequency filter
x=((1:Nx)-Nx/2)/Nx;
y=((1:Ny)-Ny/2)/Ny;
[x,y]=meshgrid(x,y);
arg=(x.^2+y.^2);
q=3; % filter edge sharpness
K=zeros([Ny Nx 3]);
for m=1:3
    lsd=1/(2*f(2*m))/(2*log(2))^(1/(2*q));
    L=exp(-.5*arg.^q/lsd^(2*q)); % low pass
    hsd=1/(2*f(2*m-1))/(2*log(2))^(1/(2*q));
    H=exp(-.5*arg.^q/hsd^(2*q));
    K(:,:,m)=g(m)*(H-L);
end
k=fftshift(sum(K,3)); % do fftshift here
% k=gpuArray(k);
clear K L H x y lsd hsd q arg Nx Ny

%% process and analyze splat images
% set up a progress bar with cancel button
wb=waitbar(0,'file:','CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
    'name','analyzing splat images');
setappdata(wb,'canceling',0)
cntr=0; % counter for progress bar
tot=length(fnm); % total # of images

dat=[];
% E={};
for n=1:tot % length(fnm)
    if getappdata(wb,'canceling') % check if cancelled
        delete(wb)
        return
    end
    
    if n>1
        S=imread(fnm{n});
        S=S(crp(2):crp(2)+crp(4)-1,crp(1):crp(1)+crp(3)-1); % crop the image
    end
    S=double(S);
    znm=['splat2D' fnm{n}(7:end)];
    if exist(znm,'file')==2 
        Z=double(imread(znm));
        Z=Z(crp(2):crp(2)+crp(4)-1,crp(1):crp(1)+crp(3)-1);
    else
        Z=S;
%         warning([znm ' not found'])
    end
    
%     S=gpuArray(S);
    F=real(ifft2(k.*fft2(S))); % apply frequency filter
    T=F<=t;
    if dst>1
        D=bwdist(T); % distance map
        L=bwlabel(D<dst);
        L(~T)=0;
    else
        L=bwlabel(T);
    end
    dat1=regionprops('table',L,Z,...
        'Area',... 
        'Centroid',...
        'Eccentricity',... % ratio of distance between ellipse foci and major axis length
        'Extent',... % ratio of pixels in the region to pixels in the total bounding box
        'Solidity',... % proportion of the pixels in the convex hull that are also in the region
        'ConvexArea',... % area of convex hull
        'MajorAxisLength',... 
        'MinorAxisLength',... 
        'Orientation',... % angle from horizontal -90 to 90
        'FilledArea',... 
        'EquivDiameter',...
        'MinIntensity',... % min depth
        'MeanIntensity',... % mean depth
        'MaxIntensity',... % max depth
        'BoundingBox'); % bounding box
%         'Perimeter',... % perimeter length (not valid for discontiguous regions)
%         'EulerNumber',... % number of objects in the region minus the number of holes in those objects
    dat1.ImageNumber=repmat(n,size(dat1,1),1);
    dat=[dat; dat1]; % concatenate data
    cntr=cntr+1;
    waitbar(cntr/tot,wb,['file: ' num2str(n) ' of ' num2str(length(fnm))]) 
end
delete(wb)

dat.Area = dat.Area*res^2;
dat.FilledArea = dat.FilledArea*res^2;
dat.ConvexArea = dat.ConvexArea*res^2;
dat.MajorAxisLength = dat.MajorAxisLength*res;
dat.MinorAxisLength = dat.MinorAxisLength*res;
dat.EquivDiameter = dat.EquivDiameter*res;
% dat.Perimeter=dat.Perimeter*res;
