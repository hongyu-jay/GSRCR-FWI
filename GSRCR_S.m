
clear ;
clc
close all;
if isempty(gcp('nocreate')), parpool; end
%% Key parameters (designed to be edited quickly)
save_data = 1;											% determine whether we save info
data_dir = '~D:\matlab\wavetools-master\SEGmodel2Dsalt/';% directory to save data in (if save_data=1)
fmin = 3;												% minimum frequency
fmax = 10;												% maximum frequency
nf = 12;												% number of frequencies
ctr = 1;												% contrast

%% Remaining parameters
tic
dx=20;
v_true = dlmread('SEGmodel2Dsalt.dat');
v_true = imresize(v_true, [200 600], 'bilinear')/1e3;
[ny nx]=size(v_true);                                   % number of gridpoints in y direction
                                                        % number of gridpoints in y direction
dom = domain([0 12 0 4],[nx ny]);						% rectangular domain
h=dom.hx;
wpml = 4*h;													% width of PML
freqs = linspace(fmin,fmax,nf);								% frequencies
ns = ceil(sqrt(nx*ny/(3*nf)));								% number of sources, in 1:3 ratio to receivers
nr = 3*ns;													% number of receivers
c_vec = [1 ctr+1];											% vector for colorbar
SNR = 0*1e-3;												% noise level 5%
sigma = SNR;												% noise level 0%
maxit = 10;													% maximum number of LBFGS iterations
smoothness = 20;											% intensity of smoothing filter
source_info.type = 'hline';									% create a horizontal line of sources
source_info.bounds = [6*h (nx-6)*h];								% left and right endpoints of sources
%source_info.bounds = [.1 2.9];								% left and right endpoints of sources
source_info.height = 3.8;									% height of sources
sources = sources_and_receivers(ns,source_info);			% x and y locations of sources
receiver_info.type = 'hline';								% create a horizontal line of receivers
receiver_info.bounds = [6*h (nx-6)*h];							% left and right endpoints of receivers
%receiver_info.bounds = [.1 2.9];							% left and right endpoints of receivers
receiver_info.height =3.8;									% height of receivers
receivers = sources_and_receivers(nr,receiver_info);		% x and y locations of receivers
window_type = 'rectangle_outer';						% change to get rectangular window around receivers
%window_type='all';
if strcmp(window_type,'all')								% decide whether or not to window
	window_info.type = 'all';								% take all points
else
	window_info.type = 'rectangle_outer';					% create a rectangular window
	window_info.bounds = [6*h (nx-6)*h source_info.height-15*h source_info.height+15*h];						% boundary of this window
end
[win_inds,W] = dom.window(window_info);						% indices of elements outside of window
pml_info.type = 'pml';										% pml info, for plotting purposes
pml_info.width = wpml;										% width of pml
[~,PML] = dom.window(pml_info);								% indicates whether a pixel is inside PML
% [v_true,ps]=mapminmax(v_true(:)',1,2);
% mapminmax('reverse',guiyiy,ps);
% v_true=1e-3*reshape(v_true,[ny nx]);
c_true = dom.mat2vec(v_true);					% true background velocity
vmin=min(c_true);
vmax=max(c_true);

v0=mean(v_true,2);
v0=repmat(v0,[1 nx]) ;
v0=smooth(v0,50);

% 
fsim =  cal_fsim(v_true,v0,0,0);
ssim =  cal_ssim(v_true,v0,0,0);
psnr =  PSNR(v_true,v0);
rms  =  RMS(v_true,v0);

fprintf('Initial FSIM= %0.2f\n',fsim);
fprintf('Initial SSIM= %0.2f\n',ssim);
fprintf('Initial PSNR = %0.2f\n',psnr);
fprintf('Initial RMS = %0.2f\n',rms);

%%
% Format and save to disk
set(gcf,'Position',[11 340 1362 462])
fprintf('Spatial dof:%d, inverse prob. dof:%d\n',dom.N,nf*ns*nr)
m_true = 1./c_true.^2;
lambda = 1.2; c1 = 8; c2 = 0.2;
alpha = 0.05;  beta = 0.05;  mu = 0.08;    err_or = 3.6E-05;
mtrue_f = gradient(m_true);
name =2;
[Opts]  =  FWI_Set(   I, name, lambda, alpha, beta, mu,  c1, c2) ;           
nw = length(win_inds);	
mm_true=dom.vec2mat(m_true);

tic
for i=1:10
            
     if(i==1)
       c0 = dom.mat2vec(v0);
       u_l=0*c0;
       q_l=0*c0;
       dm=0*c0;
       out.J=0;
       out.DJ=zeros(nw,1);
       lam=0.1;
       [m,out] = adjoint_state_2d(dom,freqs,sources,receivers,window_info,c_true,c0,SNR,maxit,u_l,q_l,v_true);
     else
         if mod(i,2)==0
       lam=lam*0.01;
         maxit = maxit + 10;
         end
       if lam <1e-5
           lam =0.5*lam;
       end
       c0=dom.mat2vec(reconstructed_image);  
       [m,out] = adjoint_state_admm(dom,freqs,sources,receivers,window_info,c_true,c0,SNR,maxit,u_l,q_l,lam,v_true);
     end

     if ~isempty(find(m<0))==1
      m(m<0.0498)=0.0498;
     end
%      

     reconstructed_image_sart = dom.vec2mat(sqrt(1./abs(m)));
     fsim =  cal_fsim(v_true,reconstructed_image_sart,0,0);
     ssim =  cal_ssim(v_true,reconstructed_image_sart,0,0);
     psnr =  PSNR(v_true,reconstructed_image_sart);
     rms  =  RMS(v_true,reconstructed_image_sart);
     fprintf('Initial fSIM= %0.3f\n',fsim);
     fprintf('Initial SSIM= %0.3f\n',ssim);
     fprintf('Initial PSNR = %0.3f\n',psnr);
     fprintf('Initial RMS = %0.3f\n',rms);

     reconstructed_image_s=dom.vec2mat(m);
     tic
     %%
     Sigma          =     sqrt(2);
     Opts.initial   =     reconstructed_image_s;
     Opts.org                   =     mm_true;
     Opts.blockSize              =  6;
     Opts.nSig                   =  Sigma; 
     E_Out_New                   =  Exter_NSS_Main  (reconstructed_image_s, Opts, Opts.blockSize);     
     Opts.nSig                   =  Sigma/255; 
     Inter_out                   =  Ixter_NSS_Main  (1*E_Out_New, Opts, Opts.nSig);
     reconstructed_m             =  ( Opts.alpha*  Opts.nSig^2* reconstructed_image_s -  Opts.alpha *  Opts.nSig^2* dom.vec2mat((q_l)) + mu * Inter_out)/ ( Opts.alpha *  Opts.nSig^2 + mu);
     toc
   
     q_l=dom.vec2mat(q_l)+(reconstructed_m)-dom.vec2mat(u_l);
     u_l=dom.mat2vec(reconstructed_m);
     
     
     reconstructed_image = sqrt(1./abs(reconstructed_m));
     
     fsim =  cal_fsim(v_true,reconstructed_image,0,0);
     ssim =  cal_ssim(v_true,reconstructed_image,0,0);
     psnr =  PSNR(v_true,reconstructed_image);
     rms  =  RMS(v_true,reconstructed_image);
     fprintf('Initial fSIM= %0.3f\n',fsim);
     fprintf('Initial SSIM= %0.3f\n',ssim);
     fprintf('Initial PSNR = %0.3f\n',psnr);
     fprintf('Initial RMS = %0.3f\n',rms);
    %%
    figure
    subplot(2,1,1)
    dom.imagesc(1e3*reconstructed_image_sart)
    caxis(1e3*[vmin,vmax])
    title('before')
    colorbar

    subplot(2,1,2)
    dom.imagesc(1e3*reconstructed_image)
    caxis(1e3*[vmin,vmax])
    title('after')
    colorbar
     
end
cpu_time = toc;
fprintf('Running time... %.2f (s).\n',cpu_time)
