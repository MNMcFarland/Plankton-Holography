function [S,D] = recsplat2(fin,nr,res,wvl,sd,w)
%
% Numerical reconstruction of in-line holograms
% Reconstructs images at depths indicated by nr
% Computes the splat image of reconstructed planes by
% using the local coefficient of variation of the gradient
% AKA the local Tamura coefficient of the gradient
%
% USAGE: [S,D]=recsplat2(fin,nr,res,wvl,sd,t,w)
%
% INPUT:
%   fin - hologram image or file name(s) (string, cell, or dir structure)
%   nr - row vector of actual (in air) plane depths in um
%   res - resolution in um/pixel, scalar or 2 element vector
%       of resolution along y (rows) and x (columns)
%   wvl - laser wavelength in um
%   sd - (optional) std devs for difference of gaussians bandpass filter in spatial
%       domain (pixels), default = [0.5 2048]
%   w - (optional) plane processing window size, default = 32
%
% OUTPUT: writes reconstructed splat image to subfolder rec
%   S - extended depth of focus splat image
%   D - depth map, plane depth for each pixel

if ~exist('sd','var') || isempty(sd); sd = [0.5 2048]; end
if ~exist('w','var') || isempty(w); w = 32; end % plane processing window size

if isnumeric(fin) && all(size(fin)>1) 
    I = fin;
    fnm = 0;
elseif isstruct(fin) && isfield(fin,'name')
    fnm = {fin.name};
elseif ~iscell(fin)
    fnm = {fin};
else
    fnm = fin;
end

% refractive index of medium alters apparent plane depth
n = 1.333; % refractive index of medium (i.e. water)
r = nr/n; % divide by refractive index of medium

if length(res)==2
    dy = res(1);  % Pixel resolution along Y (rows)
    dx = res(2);  % Pixel resolution along X (columns)
elseif length(res)==1
    dy = res;  % Pixel resolution along Y
    dx = res;  % Pixel resolution along X
else
    error('res must be scalar or a 2 element vector')
end

if nargout==0
    warning('off','MATLAB:MKDIR:DirectoryExists')
    mkdir('rec'); % create rec directory if it doesn't exist
end

% set up a progress bar with cancel button
wb=waitbar(0,'file:','CreateCancelBtn','setappdata(gcbf,''canceling'',1)',...
    'name','reconstructing hologram');
hnd = findobj(wb,'type','axes');
set(get(hnd,'title'),'interpreter','none')
setappdata(wb,'canceling',0)
cntr = 0; % counter for progress bar
tot = length(r)*length(fnm); % total # of planes to reconstruct

% read first image
if ~exist('I','var'); I = imread(fnm{1}); end

% generate the reconstruction filter kernel
[Ny,Nx] = size(I); % Ny = # rows; Nx = # columns
x = ((1:Nx)-Nx/2)/(Nx*dx);
y = ((1:Ny)-Ny/2)/(Ny*dy);
[x,y] = meshgrid(x,y);
arg = -1i*wvl*pi*(x.^2+y.^2); % for Kirchhoff-Fresnel kernel
arg = fftshift(arg);
% arg = single(arg);

% for Rayleigh-Sommerfeld kernel only
% k=2*pi/wvl;
% arg=sqrt(1-(wvl*x).^2-(wvl*y).^2);

% if isinteger(I)
%     imx=double(intmax(class(I)));
%     if (~exist('t','var') || isempty(t)); t=.1221*imx; end % threshold for reconstructed plane processing = 8000 for uint16
%     % if (~exist('t','var') || isempty(t)) && m; t=.091554*imx; end % threshold for reconstructed plane processing = 6000 for uint16
% else
%     if (~exist('t','var') || isempty(t)); t=30; end % threshold for reconstructed plane processing
% end

% % gaussian filter for smoothing
% x = ((1:Nx)-Nx/2)/Nx;
% y = ((1:Ny)-Ny/2)/Ny;
% [x,y] = meshgrid(x,y);
% fsd=1/(2*pi*sd); % std dev in frequency domain
% Gf = exp(-(x.^2+y.^2)/(2*fsd^2)); 
% Gf = fftshift(Gf);
% clear x y

% difference of gaussians bandpass filter for smoothing
x = ((1:Nx)-Nx/2)/Nx;
y = ((1:Ny)-Ny/2)/Ny;
[x,y] = meshgrid(x,y);
sd1 = 1/(2*pi*sd(1)); % standard deviation in frequency domain
G1 = exp(-(x.^2+y.^2)/(2*sd1^2));
% G1(floor(Nx/2)+1,floor(Ny/2)+1) = 0; % zero center
% Gf = gpuArray(fftshift(G1));
sd2 = 1/(2*pi*sd(2)); % standard deviation in frequency domain
G2 = exp(-(x.^2+y.^2)/(2*sd2^2));
Gf = fftshift(G1-G2);
clear x y

% Kc = ones(w); % kernel for counting pixels above t in gradient image
nbhd = true(w*2+1);
mnflt = nbhd/(w*2+1)^2; % kernel to calculate local mean with imfilter

for f=1:length(fnm)
    if f>1; I = imread(fnm{f}); end % first image already open
    
    D = zeros(Ny,Nx); % initialize depth map
%     S = zeros(Ny,Nx); % initialize splat image
    S = I; % initialize splat image
    G = zeros(Ny,Nx); % initialize array to store pixel counts
    
%     I=single(I);
    I = fft2(I);
    
    for z = r
        if getappdata(wb,'canceling') % check if cancelled
%             Si = .75*imx-S; % invert splat image
            S = uint16(S);
            delete(wb)
%             error('cancelled')
            return
        end

%         K = (2*pi*1i/k)*exp(-1i*2*pi^2*z*(x.^2 + y.^2)/k); % original
%         K = wvl*1i*exp(-1i*lambda*pi*z*(x.^2 + y.^2)); % Kirchhoff-Fresnel kernel for phase
%         K = exp(-1i*wvl*pi*z*(x.^2 + y.^2)); % Kirchhoff-Fresnel kernel
        K = exp(z*arg); % Kirchhoff-Fresnel kernel
%         K = exp(-1i*k*z*arg); % Rayleigh-Sommerfeld kernel
%         K = ifftshift(K);
%         K = single(K);
        M = I.*K;
        
        % process reconstructed image to determine correct plane
%         M1 = imgaussfilt(M,4,'Padding','symmetric'); % smooth
        M1 = M.*Gf; % smooth
        M1 = real(ifft2(M1));
        
        M1 = imgradient(M1,'central'); % gradient magnitude
%         M1 = imfilter(double(M1>t), Kc, 'symmetric'); % count pixels > threshold within local nhood (w x w)
%         M1(M1<t) = 0;
%         M1 = imfilter(M1, Kc, 'symmetric'); % sum pixels within local nhood (w x w)
%         M1 = imdilate(M1,Kc);
%         M1 = stdfilt(M1,nbhd) ./ imfilter(M1,mnflt,'symmetric'); % Tamura coefficient of the gradient (squared)
        Ex = imfilter(M1,mnflt,'symmetric'); % local mean of gradient, E[X]
        Ex2 = imfilter(M1.^2,mnflt,'symmetric'); % local mean of squared gradient, E[X^2]
        M1 = sqrt(Ex2 - Ex.^2)./Ex; % local coefficient of variation of gradient = std/mean
        
        ind = M1>=G; % find max pixel count
        G(ind) = M1(ind); % store max pixel count
        
        M = ifft2(M);
%         Ph = abs(atan(imag(M)./real(M))); % phase? (only seems to work with K-F kernel)
        M = real(M);
        
        S(ind) = M(ind); % update Splat image
        D(ind) = z*n; % update the depth map
        
        % update progress bar
        cntr = cntr+1;
        waitbar(cntr/tot,wb,['file: ' num2str(f) ' of ' num2str(length(fnm))...
            '; z = ' num2str(z*n)]) 
    end
%     if length(r)==2
%         assignin('base','M2',M2); % gradient image for optimizing sd and t 
%     end
    if nargout==0 && ~isnumeric(fin)
%         Si = .75 * imx - S; % invert splat image
%         Si = uint16(Si);
        S = S/4095*65535; % scale splat from 12 bit to 16 bit
        [~,fstr] = fileparts(fnm{f});
        fstr = [fstr '_' num2str(min(nr)/1000) '-' num2str(max(nr)/1000)]; % append depth range to file name in mm
%         imwrite(Si,['rec\splat2_' fstr '.tif'],'tif','compression','none')
%         imwrite(uint16(D),['rec\splat2D_' fstr '.tif'],'tif','compression','none')
        imwrite(uint16(S),['rec\splat2_' fstr '.png'],'png')
        imwrite(uint16(D),['rec\splat2D_' fstr '.png'],'png')
    end
end
delete(wb)
warning('on','MATLAB:MKDIR:DirectoryExists')
