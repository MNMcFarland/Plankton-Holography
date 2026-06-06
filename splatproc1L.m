function [dat,L]=splatproc1L(fnm,res,t,f,g,dst,crp)

% process and analyze splat images
% with image cropping
% Uses Lanalyze.m
%
% USAGE: dat=splatproc2(fnm,res,t,f,g,dst,crp)
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
%   dat - data for all input files
%       column -  data type
%       -------------------------------
%          1   -  Unfilled area (um^2)
%          2   -  Min depth
%          3   -  Max depth
%          4   -  Mean depth
%          5   -  Min diameter (um)
%          6   -  Max diameter (um)
%          7   -  Mean diameter (um)
%          8   -  Orientation, degrees from vertical
%          9   -  Edge contact (logical)
%          10  -  Image number

%   L - label matrix
%   E - region outlines (cell)

if isstruct(fnm) && isfield(fnm,'name')
    fnm={fnm.name};
elseif ~iscell(fnm)
    fnm={fnm};
end
if ~exist('res','var') || isempty(res); res=1; end
% if ~exist('t','var') || isempty(t); t=10000; end
if ~exist('f','var') || isempty(f); f=[3 10 10 30 30 150]; end
if ~exist('g','var') || isempty(g); g=[.8 1.1 1.3]; end
if ~exist('dst','var') || isempty(dst); dst=10; end
% t=gpuArray(t);

S=imread(fnm{1});
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
k=ifftshift(sum(K,3)); % do ifftshift here
% k=gpuArray(k);
clear K L H x y

%% process and analyze splat images
% set up a progress bar with cancel button
wb=waitbar(0,'file:','CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
    'name','analyzing splat images');
setappdata(wb,'canceling',0)
cntr=0; % counter for progress bar
tot=length(fnm); % total # of planes to reconstruct

dat=table;
for n=1:length(fnm)
    if getappdata(wb,'canceling') % check if cancelled
        delete(wb)
        return
    end
    
    if n>1
        S=imread(fnm{n});
        S=S(crp(2):crp(2)+crp(4)-1,crp(1):crp(1)+crp(3)-1);
    end
    S=double(S);
    znm=['splatD' fnm{n}(7:end)];
    if exist(znm,'file')==2
        Z=double(imread(znm));
        Z=Z(crp(2):crp(2)+crp(4)-1,crp(1):crp(1)+crp(3)-1);
    else
        Z=S;
%         warning([znm ' not found'])
    end
    
%     S=gpuArray(S);
    F=real(ifft2(k.*fft2(S))); % apply frequency filter
    T=F>=t;
    if dst>0
        D=bwdist(T); % distance map
        L=bwlabel(D<dst); % label expanded regions
        L(~T)=0; % contract regions
    else
        L=bwlabel(T);
    end
    
    dat1=Lanalyze(L,Z,15); % skips some measurements for regions < 15 pixels = 20 um ESD
%     dat1=dat1(dat1.Edge,:); % avoid particles touching the image border
    dat1.ImageNumber=repmat(n,size(dat1,1),1);
    dat=[dat; dat1]; % concatenate data
    cntr=cntr+1;
    waitbar(cntr/tot,wb,['file: ' num2str(n) ' of ' num2str(length(fnm))]) 
end
delete(wb)

% dat.Orientation=90-dat.Orientation; % switch to angle from horizontal
dat.Area=dat.Area*res^2;
dat.FilledArea=dat.FilledArea*res^2;
dat.MajorAxisLength=dat.MajorAxisLength*res;
dat.MinorAxisLength=dat.MinorAxisLength*res;
dat.MaxFeret=dat.MaxFeret*res;
dat.MedianScan=dat.MedianScan*res;
dat.MeanScan=dat.MeanScan*res;
dat.Perimeter=dat.Perimeter*res;
dat.FilledArea=dat.FilledArea*res^2;
% dat.ConvexArea=dat.ConvexArea*res^2;

