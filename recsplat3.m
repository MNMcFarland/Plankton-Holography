function [S,D] = recsplat3(fin,nr,res,wvl,w)
%
% Numerical reconstruction of in-line holograms
% Reconstructs images at depths indicated by nr
% Computes the splat image of reconstructed planes by
% using the local coefficient of variation of the gradient
% AKA the local Tamura coefficient of the gradient
%
% USAGE: [S,D]=recsplat3(fin,nr,res,wvl,w)
%
% INPUT:
%   fin - hologram image or file name(s) (string, cell, or dir structure)
%   nr - row vector of actual (in air) plane depths in um
%   res - resolution in um/pixel, scalar or 2 element vector
%       of resolution along y (rows) and x (columns)
%   wvl - laser wavelength in um
%   w - (optional) plane processing window size, default = 32
%
% OUTPUT: writes reconstructed splat image to subfolder rec
%   S - extended depth of focus splat image
%   D - depth map, plane depth for each pixel

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
    mkdir('rec3'); % create rec directory if it doesn't exist
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

nbhd = true(w*2+1);
mnflt = nbhd/(w*2+1)^2; % kernel to calculate local mean with imfilter

for f=1:length(fnm)
    if f>1; I = imread(fnm{f}); end % first image already open
    
    D = zeros(Ny,Nx); % initialize depth map
%     S = zeros(Ny,Nx); % initialize splat image
    S = I; % initialize splat image
    F = zeros(Ny,Nx); % initialize array to store pixel counts
    
    I = fft2(I);
    
    for z = r
        if getappdata(wb,'canceling') % check if cancelled
            S = uint16(S);
            D = uint16(D);
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
        P = I.*K;
        
        P = ifft2(P);
%         Ph = abs(atan(imag(P)./real(P))); % phase? (only seems to work with K-F kernel)
        P = real(P);

        % process reconstructed image to determine correct plane
        G = imgradient(P,'central'); % gradient magnitude
%         G = stdfilt(G,nbhd) ./ imfilter(G,mnflt,'symmetric'); % Tamura coefficient of the gradient (squared)
        Eg = imfilter(G,mnflt,'symmetric'); % local mean of gradient, E[G]
        Eg2 = imfilter(G.^2,mnflt,'symmetric'); % local mean of squared gradient, E[G^2]
        G = sqrt(Eg2 - Eg.^2)./Eg; % local coefficient of variation of gradient = std/mean
        
        ind = G>=F; % find max pixel count
        F(ind) = G(ind); % store max pixel count
        
        S(ind) = P(ind); % update Splat image
        D(ind) = z*n; % update the depth map
        
        % update progress bar
        cntr = cntr+1;
        waitbar(cntr/tot,wb,['file: ' num2str(f) ' of ' num2str(length(fnm))...
            '; z = ' num2str(z*n)]) 
    end
    if nargout==0 && ~isnumeric(fin)
        S = S/4095*65535; % scale splat from 12 bit to 16 bit
        [~,fstr] = fileparts(fnm{f});
        fstr = [fstr '_' num2str(min(nr)/1000) '-' num2str(max(nr)/1000)]; % append depth range to file name in mm
        imwrite(uint16(S),['rec3\splat3_' fstr '.png'],'png')
        imwrite(uint16(D),['rec3\splat3D_' fstr '.png'],'png')
    end
end
delete(wb)
warning('on','MATLAB:MKDIR:DirectoryExists')
