function recgui_nc(varargin)

% interactive holographic reconstruction
% opens .nc files
% Requires: holosub_nc.m, recsplat2.m, and image processing toolbox
% 
% USAGE: recgui_nc
%   recgui_nc(I)
%   recgui_nc('res',4.584,'wvl',0.660,'maxz',20)
%
% INPUT:
%   I - (optional) hologram image
%   'res',resolution - (optional) name value pair um/pixel 
%   'wvl',wavelength - (optional) name value pair in microns 
%   'minz',minimum reconstruction depth - (optional) name value pair in millimeters
%   'maxz',maximum reconstruction depth - (optional) name value pair in millimeters
%   'zstep',reconstruction interval - (optional) name value pair in millimeters
%

% initialize data structure and set defaults
dat = struct;
dat.res = 2.5; % set default resolution um/pixel
dat.wvl = 0.532; % set default wavelength
dat.maxz = 12;
dat.minz = 0;
dat.zstep = 0.5;
dat.inum = 1;
dat.z = dat.minz;
dat.sd1 = .5; % stdev for smoothing
dat.sd2 = 2048; % stdev for smoothing
dat.w = 64; % window radius for splat

% parse input
for m=1:length(varargin)
    if all(size(varargin{m})>1)
        dat.H = varargin{m}; % put image in dat.H
        dat.name = [];
        if isinteger(dat.H)
            dat.imx = double(intmax(class(dat.H)));
        else
            dat.imx = max(dat.H,[],'all');
        end
    elseif ischar(varargin{m}) && strcmp(varargin{m},'res')
        m=m+1;
        dat.res=double(varargin{m});
    elseif ischar(varargin{m}) && strcmp(varargin{m},'wvl')
        m=m+1;
        dat.wvl=double(varargin{m});
    elseif ischar(varargin{m}) && strcmp(varargin{m},'minz')
        m=m+1;
        dat.minz=double(varargin{m});
    elseif ischar(varargin{m}) && strcmp(varargin{m},'maxz')
        m=m+1;
        dat.maxz=double(varargin{m});
    elseif ischar(varargin{m}) && strcmp(varargin{m},'zstep')
        m=m+1;
        dat.zstep=double(varargin{m});
    end
end

%% setup figure
dat.fh = figure('position',[40 45 1200 835],...
    'numbertitle','off',...
    'name','autoholo reconstruction',...
    'resize','on',...
    'menubar','none',...
    'color',[1 1 1]);

% setup image axes
dat.ax = axes('position',[0 .1 1 .9],'xtick',[],'ytick',[],'box','off','color','k');
if isfield(dat,'H'); dat = load_hologram(dat); end 

% Add menus with Accelerators
mymenu = uimenu('Parent',dat.fh,'Label','Tools');
uimenu('Parent',mymenu,'Label','Zoom','Accelerator','z','Callback',@(src,evt)zoom(dat.fh,'on'));
uimenu('Parent',mymenu,'Label','Pan','Accelerator','d','Callback',@(src,evt)pan(dat.fh,'on'));
% uimenu('Parent',mymenu,'Label','Rotate','Accelerator','r','Callback',@(src,evt)rotate3d(dat.fh,'on'));

%% setup toolbar
icon = 'baseline_add_photo_alternate_white_24dp.png';
[cdata,~,alpha] = imread(icon);
cdata = double(cdata)/255.*(double(alpha)/255);
cdata(cdata==0) = NaN;
cdata = 1-repmat(cdata,1,1,3);
uipushtool('tooltipstring','open .nc file','clickedcallback',@openfile_Callback,...
    'cdata',cdata);

icon = 'baseline_skip_next_white_24dp.png';
[cdata,~,alpha] = imread(icon);
cdata = fliplr(cdata);
alpha = fliplr(alpha);
cdata = double(cdata)/255.*(double(alpha)/255);
cdata(cdata==0) = NaN;
cdata = 1-repmat(cdata,1,1,3);
uipushtool('tooltipstring','previous image','clickedcallback',@imsel_Callback,...
    'cdata',cdata,'UserData',-1);

icon = 'baseline_skip_next_white_24dp.png';
[cdata,~,alpha] = imread(icon);
cdata = double(cdata)/255.*(double(alpha)/255);
cdata(cdata==0) = NaN;
cdata = 1-repmat(cdata,1,1,3);
uipushtool('tooltipstring','next image','clickedcallback',@imsel_Callback,...
    'cdata',cdata,'UserData',1);

icon = 'baseline_photo_library_white_24dp.png';
[cdata,~,alpha] = imread(icon);
cdata = double(cdata)/255.*(double(alpha)/255);
cdata(cdata==0) = NaN;
cdata = 1-repmat(cdata,1,1,3);
uipushtool('tooltipstring','subtract background','clickedcallback',@bkgrnd_Callback,...
    'cdata',cdata);

icon = 'baseline_zoom_in_white_24dp.png';
[cdata,~,alpha] = imread(icon);
cdata = double(cdata)/255.*(double(alpha)/255);
cdata(cdata==0) = NaN;
cdata = 1-repmat(cdata,1,1,3);
uipushtool('tooltipstring','toggle zoom','clickedcallback','zoom','cdata',cdata);

icon = 'baseline_pan_tool_white_18dp.png';
[cdata,~,alpha] = imread(icon);
cdata = double(cdata)/255.*(double(alpha)/255);
cdata(cdata==0) = NaN;
cdata = 1-repmat(cdata,1,1,3);
uipushtool('tooltipstring','toggle pan','clickedcallback','pan','cdata',cdata);

icon = 'baseline_save_white_24dp.png';
[cdata,~,alpha] = imread(icon);
cdata = double(cdata)/255.*(double(alpha)/255);
cdata(cdata==0) = NaN;
cdata = 1-repmat(cdata,1,1,3);
uipushtool('tooltipstring','save image','clickedcallback',@put_Callback,...
    'cdata',cdata);

icon = 'baseline_play_circle_filled_white_24dp.png';
[cdata,~,alpha] = imread(icon);
cdata = double(cdata)/255.*(double(alpha)/255);
cdata(cdata==0) = NaN;
cdata = 1-repmat(cdata,1,1,3);
uipushtool('tooltipstring','process current z plane','clickedcallback',@plnproc_Callback,...
    'cdata',cdata);

icon = 'baseline_play_circle_outline_white_24dp.png';
[cdata,~,alpha] = imread(icon);
cdata = double(cdata)/255.*(double(alpha)/255);
cdata(cdata==0) = NaN;
cdata = 1-repmat(cdata,1,1,3);
uipushtool('tooltipstring','splat image','clickedcallback',@splat_Callback,...
    'cdata',cdata);

icon = 'baseline_settings_applications_white_24dp.png';
[cdata,~,alpha] = imread(icon);
cdata = double(cdata)/255.*(double(alpha)/255);
cdata(cdata==0) = NaN;
cdata = 1-repmat(cdata,1,1,3);
uipushtool('tooltipstring','settings','clickedcallback',@settings_Callback,...
    'cdata',cdata);

icon = 'baseline_bathtub_white_24dp.png';
[cdata,~,alpha] = imread(icon);
cdata = double(cdata)/255.*(double(alpha)/255);
cdata(cdata==0) = NaN;
cdata = 1-repmat(cdata,1,1,3);
uipushtool('tooltipstring','save configuration','clickedcallback',@save_config_Callback,...
    'cdata',cdata);

icon = 'baseline_folder_white_24dp.png';
[cdata,~,alpha] = imread(icon);
cdata = double(cdata)/255.*(double(alpha)/255);
cdata(cdata==0) = NaN;
cdata = 1-cdata;
uipushtool('tooltipstring','load configuration','clickedcallback',@load_config_Callback,...
    'cdata',cdata);

%% setup controls
fs = 10; % fontsize

dat.dslider = uicontrol('style','slider',...
    'Units','normalized',...
    'Position',[.05 .05 .9 .025],...
    'Callback',@dslider_Callback,...
    'min',dat.minz,'max',dat.maxz,...
    'sliderstep',[dat.zstep/(dat.maxz-dat.minz) dat.zstep/(dat.maxz-dat.minz)],...
    'value',dat.minz);

dat.minzfld = uicontrol('style','edit',...
    'Units','normalized',...
    'Position',[.05 .01 .05 .025],...
    'Callback',@minz_Callback,...
    'BackgroundColor',[1 1 1],...
    'String',num2str(dat.minz),...
    'fontsize',fs);
annotation('textbox',[0 .01 .05 .025],'String','min z =',...
    'linestyle','none','fontsize',fs,'horizontalalignment','right');
annotation('textbox',[.1 .01 .04 .025],'String','mm',...
    'linestyle','none','fontsize',fs,'horizontalalignment','left');

dat.stpfld = uicontrol('style','edit',...
    'Units','normalized',...
    'Position',[.225 .01 .05 .025],...
    'Callback',@stp_Callback,...
    'BackgroundColor',[1 1 1],...
    'String',num2str(dat.zstep),...
    'fontsize',fs);
annotation('textbox',[0.175 .01 .055 .025],'String','z step =',...
    'linestyle','none','fontsize',fs,'horizontalalignment','right');
annotation('textbox',[.275 .01 .04 .025],'String','mm',...
    'linestyle','none','fontsize',fs,'horizontalalignment','left');

dat.depthfld = uicontrol('style','edit',...
    'Units','normalized',...
    'Position',[.4 .01 .05 .025],...
    'Callback',@depth_Callback,...
    'BackgroundColor',[1 1 1],...
    'String',num2str(dat.dslider.Value),...
    'fontsize',fs);
annotation('textbox',[.34 .01 .06 .025],'String','depth z =',...
    'linestyle','none','fontsize',fs,'horizontalalignment','right');
annotation('textbox',[.45 .01 .04 .025],'String','mm',...
    'linestyle','none','fontsize',fs,'horizontalalignment','left');

dat.bkgbtn = uicontrol('style','togglebutton',...
    'Units','normalized',...
    'Position',[.515 .01 .15 .025],...
    'Callback',@bkgbtn_Callback,...
    'BackgroundColor',[1 1 1],...
    'String','background disabled',...
    'fontsize',fs,...
    'Max',true,...
    'Min',false,...
    'Value',false);

dat.inumfld = uicontrol('style','edit',...
    'Units','normalized',...
    'Position',[.75 .01 .05 .025],...
    'Callback',@inum_Callback,...
    'BackgroundColor',[1 1 1],...
    'String',num2str(dat.inum),...
    'fontsize',fs);
annotation('textbox',[0.695 .01 .055 .025],'String','image =',...
    'linestyle','none','fontsize',fs,'horizontalalignment','right');
dat.nimtxt = annotation('textbox',[.8 .01 .05 .025],'String','of 0',...
    'linestyle','none','fontsize',fs,'horizontalalignment','left');

dat.maxzfld = uicontrol('style','edit',...
    'Units','normalized',...
    'Position',[.9 .01 .05 .025],...
    'Callback',@maxz_Callback,...
    'BackgroundColor',[1 1 1],...
    'String',num2str(dat.maxz),...
    'fontsize',fs);
annotation('textbox',[.85 .01 .055 .025],'String','max z =',...
    'linestyle','none','fontsize',fs,'horizontalalignment','right');
annotation('textbox',[.95 .01 .04 .025],'String','mm',...
    'linestyle','none','fontsize',fs,'horizontalalignment','left');

guidata(dat.fh,dat);

%% toolbar callbacks --------------------------------------------------------------------------------------
function openfile_Callback(src,~)
    dat = guidata(src);
    [fname,fpath] = uigetfile({'*.nc'},'multiselect','off');
    if ischar(fname)
        nci = ncinfo([fpath fname]);
        dat.group = nci.Groups.Name;
        dat.name = fname;
        dat.path = fpath;
        dim = [nci.Groups.Dimensions.Length];
        dat.nim = dim(1);
        dat.dim = [dim(2) dim(3) 1]; % [dim(2) dim(3) dim(1)];
        dat.H = ncread([fpath fname],[dat.group '/Images'],[1 1 dat.inum],dat.dim)';
        dat.H = double(dat.H);
        dat.imx = 4095;
        dat.nimtxt.String = ['of ' num2str(dat.nim)];
        
        dat = load_hologram(dat);
        set(dat.fh,'name',['Dataset - ' fname ' - ' nci.Attributes.Value])
    end
    guidata(src,dat);
    
function imsel_Callback(src,~)
    dat = guidata(src);
    inum = dat.inum + src.UserData;
    if inum == round(inum) && inum>0 && inum<=dat.nim
        dat.inum = inum;
        dat.H = ncread([dat.path dat.name],[dat.group '/Images'],[1 1 dat.inum],dat.dim)';
        dat.H = double(dat.H);
        dat = show_hologram(dat);
        set(dat.inumfld,'String',num2str(dat.inum));
    else
        warndlg('image number out of range','Warning','modal')
    end
    guidata(src,dat);

function bkgrnd_Callback(src,~)
    dat = guidata(src);
    [fname,fpath] = uigetfile({'*.nc'},'multiselect','on');
    if ischar(fname) || iscell(fname)
        if iscell(fname)
            fname = cellfun(@(x) [fpath x],fname,'uniformoutput',false); % prepend path
        elseif ischar(fname)
            fname = [fpath fname];
        end
            dat.Av = holosub_nc([],fname,50); % generate background image
            dat.Av = dat.Av - mean(dat.Av,'all'); % zero center
            dat.bkgbtn.Value = true;
            dat.bkgbtn.String = 'background enabled';
            if isfield(dat,'H')
                dat = show_hologram(dat);
            end
    end
    guidata(src,dat);

function put_Callback(src,~)
    dat=guidata(src);
    if isfield(dat,'M')
        [fnm,pth] = uiputfile({'*.tif';'*.jpg';'*.bmp';'*.png';'*.gif'});
        I = get(dat.im,'cdata');
        assignin('base','I',I);
        if ~(fnm==0)
            I = uint8(I/dat.imx*255);
            imwrite(I,[pth fnm]);
        end
    end
    guidata(src,dat);

function plnproc_Callback(src,~)
    dat=guidata(src);
    if isfield(dat,'I')
        P = planeproc(dat);
        set(dat.im,'cdata',P);
%         dat.ax.CLim = [0 max(P,[],'all')];
        dat.ax.CLim = [min(P,[],'all') max(P,[],'all')];
%         dat.cb = colorbar;
    end
    guidata(src,dat)

function splat_Callback(src,~)
    dat=guidata(src);
    if isfield(dat,'H')
        if isfield(dat,'Av') && dat.bkgbtn.Value
            A = dat.Av;
        else
            A = 0;
        end
        nr = (dat.minz:dat.zstep:dat.maxz)*1000;
%         S = recsplatgpu3(dat.H-A,nr,dat.res,dat.wvl,dat.w);
        [S,D] = recsplatgpu2(dat.H-A,nr,dat.res,dat.wvl,[dat.sd1 dat.sd2],dat.w);
%         [S,D] = recsplat2(dat.H-A,nr,dat.res,dat.wvl,[dat.sd1 dat.sd2],dat.w);
        set(dat.im,'cdata',S);
    else
        errordlg('load hologram','image error');
    end
    guidata(src,dat)

function settings_Callback(src,~)
    dat = guidata(src);
    definpt = {num2str(dat.res), num2str(dat.wvl), num2str(dat.sd1), num2str(dat.sd2), num2str(dat.w)};
    stg = inputdlg({'resolution (um/pixel):','wavelength (um):','high freq. cutoff (pixels):',...
        'low freq. cutoff (pixels):','window radius (pixels):'},'Settings',[1 35],definpt);
    if numel(stg)>0
        dat.res = str2double(stg{1});
        dat.wvl = str2double(stg{2});
        dat.sd1 = str2double(stg{3});
        dat.sd2 = str2double(stg{4});
        dat.w = str2double(stg{5});
        if isfield(dat,'H')
            dat = load_hologram(dat);
        end
    end
    guidata(src,dat)
    
function save_config_Callback(src,~)
    dat = guidata(src);
    cfg.res = dat.res;
    cfg.wvl = dat.wvl;
    cfg.sd1 = dat.sd1;
    cfg.sd2 = dat.sd2;
    cfg.w = dat.w;
    cfg.minz = dat.minz;
    cfg.maxz = dat.maxz;
    cfg.zstep = dat.zstep;
    cfg.date = datestr(now);
    [fnm,pth] = uiputfile('recgui.rgc','Save Settings');
    if fnm~=0
        save([pth fnm],'cfg','-mat')
    end
    guidata(src,dat)

function load_config_Callback(src,~)
    dat = guidata(src);
    [fnm,pth] = uigetfile('*.rgc','Load Settings');
    if fnm~=0
        load([pth fnm],'-mat','cfg')
        dat.res = cfg.res;
        dat.wvl = cfg.wvl;
        dat.sd1 = cfg.sd1;
        dat.sd2 = cfg.sd2;
        dat.w = cfg.w;
        dat.minz = cfg.minz;
        dat.maxz = cfg.maxz;
        dat.zstep = cfg.zstep;
        dat.z = dat.minz;
        guidata(src,dat)
        
        dat.stpfld.String = num2str(dat.zstep);
        stp_Callback(dat.stpfld);
        dat.minzfld.String = num2str(dat.minz);
        minz_Callback(dat.minzfld)
        dat.maxzfld.String = num2str(dat.maxz);
        maxz_Callback(dat.maxzfld)
        dat.depthfld.String = num2str(dat.z);
        depth_Callback(dat.depthfld)
        
        if isfield(dat,'H')
            dat = load_hologram(dat);
        end
    end
    guidata(src,dat)

%% control callbacks -----------------------------------------------------------------------------------
function dslider_Callback(src,~)
    dat = guidata(src);
    dat.z = get(src,'Value');
    set(dat.depthfld,'String',num2str(dat.z));
    if isfield(dat,'H')
        dat.M = recfun(dat);
        set(dat.im,'cdata',dat.M);
        dat.ax.CLim = [0 dat.imx];
    end
    guidata(src,dat)

function minz_Callback(src,~)
    dat = guidata(src);
    dat.minz = abs(str2double(src.String));
    if dat.minz>=dat.maxz
        dat.maxz = dat.minz+dat.zstep;
        dat.dslider.Max = dat.maxz;
        dat.maxzfld.String = num2str(dat.maxz);
    end
    if dat.minz>=dat.dslider.Value
        dat.dslider.Value = dat.minz;
        dslider_Callback(dat.dslider)
    end
    dat.dslider.Min = dat.minz;
    dat.dslider.SliderStep = [dat.zstep/(dat.maxz-dat.minz) dat.zstep/(dat.maxz-dat.minz)];
    guidata(src,dat)

function stp_Callback(src,~)
    dat = guidata(src);
    stp = get(src,'String');
    dat.zstep = abs(str2double(stp));
    if dat.zstep>dat.maxz-dat.minz
        dat.maxz = dat.minz + dat.zstep;
        dat.dslider.Max = dat.maxz;
        dat.maxzfld.String = num2str(dat.maxz);
    end
    dat.dslider.SliderStep = [dat.zstep/(dat.maxz-dat.minz) dat.zstep/(dat.maxz-dat.minz)];
    guidata(src,dat)

function depth_Callback(src,~)
    dat = guidata(src);
    dat.z = str2double(get(src,'String'));
    if dat.z>dat.maxz
        dat.dslider.Max = dat.z;
        dat.maxz = dat.z;
        dat.maxzfld.String = num2str(dat.z);
    elseif dat.z<dat.minz
        dat.minz = dat.z;
        dat.dslider.Min = dat.z;
        dat.minzfld.String = num2str(dat.z);
    end
    dat.dslider.Value = dat.z;
    if isfield(dat,'H')
        dat.M=recfun(dat);
        set(dat.im,'cdata',dat.M);
        dat.ax.CLim = [0 dat.imx];
    end
    guidata(src,dat);

function bkgbtn_Callback(src,~)    
    dat = guidata(src);
    if src.Value
        src.String = 'background enabled';
    else
        src.String = 'background disabled';
    end
    if isfield(dat,'H')
        dat = show_hologram(dat);
    end
    guidata(src,dat);

function inum_Callback(src,~)
    dat = guidata(src);
    inum = str2double(src.String);
    if inum==round(inum) && inum>0 && inum<=dat.nim
        dat.inum = inum;
        dat.H = ncread([dat.path dat.name],[dat.group '/Images'],[1 1 dat.inum],dat.dim)';
        dat.H = double(dat.H);
%         dat.H = double(dat.H) - dat.Av;
%         dat.H = double(dat.H)./dat.Av;
        dat = show_hologram(dat);
    else
        warndlg('image number out of range','Warning','modal')
        src.String = dat.inum;
    end
    guidata(src,dat)

function maxz_Callback(src,~)
    dat = guidata(src);
    dat.maxz = abs(str2double(src.String));
    if dat.maxz<=dat.minz
        dat.minz = dat.maxz-dat.zstep;
        dat.dslider.Min = dat.minz;
        dat.minfld.String = num2str(dat.minz);
    end
    if dat.dslider.Value>=dat.maxz
        dat.dslider.Value = dat.maxz;
        dslider_Callback(dat.dslider);
    end
    dat.dslider.Max = dat.maxz;
    dat.dslider.SliderStep = [dat.zstep/(dat.maxz-dat.minz) dat.zstep/(dat.maxz-dat.minz)];
    guidata(src,dat)

%% functions ---------------------------------------------------------------------------------
function dat = load_hologram(dat)
%     dat.H = double(dat.H);
    
    [Ny,Nx] = size(dat.H); % Image Size
    x = ((1:Nx)-Nx/2)/(Nx*dat.res);
    y = ((1:Ny)-Ny/2)/(Ny*dat.res);
    [x,y] = meshgrid(x,y);
    dat.arg = -1i*dat.wvl*pi*(x.^2+y.^2); % Kirchhoff-Fresnel kernel only
%     dat.k = 2*pi/dat.wvl; % for Rayleigh-Sommerfeld Kernel only
%     dat.arg = sqrt(1-(dat.wvl*dat.x).^2-(dat.wvl*dat.y).^2); % for Rayleigh-Sommerfeld kernel only
    
    % difference of gaussians filter for smoothing
    x = ((1:Nx)-Nx/2)/Nx;
    y = ((1:Ny)-Ny/2)/Ny;
    [x,y] = meshgrid(x,y);
    sd1 = 1/(2*pi*dat.sd1); % standard deviation in frequency domain
    G1 = exp(-(x.^2+y.^2)/(2*sd1^2));
    sd2 = 1/(2*pi*dat.sd2); % standard deviation in frequency domain
    G2 = exp(-(x.^2+y.^2)/(2*sd2^2));
    dat.Gf = fftshift(G1-G2); % band pass filter
    
    dat = show_hologram(dat);

function dat = show_hologram(dat)
    if dat.bkgbtn.Value && (~isfield(dat,'Av') || any(size(dat.Av)~=size(dat.H)))
        dat.bkgbtn.Value = false;
        dat.bkgbtn.String = 'background disabled';
        warning('background image not defined or size does not match hologram')
    end
    
    if dat.bkgbtn.Value
        dat.I = fft2(dat.H - dat.Av);
    else
        dat.I = fft2(dat.H); % pre-compute FFT of hologram
    end
    dat.M = recfun(dat);

    if isfield(dat,'im')
        set(dat.im,'cdata',dat.M);
        dat.ax.CLim = [0 dat.imx];
    else
        dat.im = imagesc(dat.ax,dat.M,[0 dat.imx]);
        axis image
        set(dat.ax,'visible','off')
        colormap(dat.ax,gray)
    end
    
function M = recfun(dat)
    z = dat.z*1000/1.333; % account for refractive index of medium (water = 1.333)
%     z = dat.z*1000; % don't account for refractive index of medium (water = 1.333)
%     n = 1i*dat.wvl*exp(-1i*dat.wvl*pi*dat.z*(dat.x.^2 + dat.y.^2)); % full Kirchhoff-Fresnel kernel
    n = exp(z*dat.arg); % Kirchhoff-Fresnel kernel
%     n = exp(-1i*dat.k*dat.z*dat.arg); % Rayleigh-Sommerfeld kernel
    n = fftshift(n);
    M = ifft2(dat.I.*n);
    M = real(M);

function P = planeproc(dat)
    z = dat.z*1000/1.333; % account for refractive index of medium (water = 1.333)
    n = exp(z*dat.arg); % Kirchhoff-Fresnel kernel
    n = fftshift(n);
    P = dat.I.*n; % reconstruct
    P = P.*dat.Gf; % smooth and zero center
    P = real(ifft2(P));
%     P(P>0) = 0;
    P = imgradient(P,'central');
    w = dat.w*2+1;
    SD = stdfilt(P,true(w)); % local standard deviation
    M = imfilter(P,ones(w)/w^2,'symmetric'); % local mean
    P = SD./M; % coefficient of variation of the gradient
%     threshold(P);
    
