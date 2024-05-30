
close all
clc
if isempty(gcp('nocreate')), parpool; end

%% Key parameters (designed to be edited quickly)
window_type = 'all';						           % change to get rectangular window around receivers
save_data = 1;											% determine whether we save info
data_dir = '~D:\matlab\wavetools-master\marm_data/';	% directory to save data in (if save_data=1)
fmin = 3;												% minimum frequency
fmax = 15;												% maximum frequency
nf = 12;												% number of frequencies
ctr = 1;												% contrast

%% Remaining parameters
nx = 640;									% number of gridpoints in x direction8
ny = 204;											% number of gridpoints in y direction8
% h=4.8*1e-3;
dom = domain([0 9.6 0 3.05],[nx ny]);				% rectangular domain
h=dom.hx;									        % width of PML
wpml = 5*h;	
freqs = linspace(fmin,fmax,nf);								% frequencies
ns = ceil(sqrt(nx*ny/(3*nf)));								% number of sources, in 1:3 ratio to receivers							% frequencies
nr = 3*ns;													% number of receivers
c_vec = [1 ctr+1];											% vector for colorbar
SNR =8*1e-3;    										    % noise level 5%
% SNR = 1/eps;												% noise level 0%
maxit = 15;													% maximum number of LBFGS iterations
smoothness =15;											    % intensity of smoothing filter
source_info.type = 'hline';									% create a horizontal line of sources
source_info.bounds = [2*h (nx-2)*h];						% left and right endpoints of sources
source_info.height = 2.9;                           		% height of sources
sources = sources_and_receivers(ns,source_info);			% x and y locations of sources
receiver_info.type = 'hline';								% create a horizontal line of receivers
receiver_info.bounds = [2*h (nx-2)*h];					    % left and right endpoints of receivers
receiver_info.height = 2.9;                                 % height of receivers
receivers = sources_and_receivers(nr,receiver_info);		% x and y locations of receivers
%window_type='all';

if strcmp(window_type,'all')								% decide whether or not to window
	window_info.type = 'all';								% take all points
else
	window_info.type = 'rectangle_outer';					% create a rectangular window
	window_info.bounds = [1*h (nx-1)*h 2.85 2.95];	        %boundary of this window
end
[win_inds,W] = dom.window(window_info);						% indices of elements outside of window
nw = length(win_inds);	
pml_info.type = 'pml';										% pml info, for plotting purposes
pml_info.width = wpml;										% width of pml
[~,PML] = dom.window(pml_info);								% indicates whether a pixel is inside PML
v_true = 1*marmousi(dom,ctr);					                % true background velocity
v_true=reshape(v_true,[ny nx]);
c_true = dom.mat2vec(v_true);				               	% true background velocity
v0=velsmooth(v_true,15,15,1200);
m_true = 1./c_true.^2;	
vmin=min(c_true);
vmax=max(c_true);

fsim =  cal_fsim(v_true,v0,0,0);
ssim =  cal_ssim(v_true,v0,0,0);
psnr = PSNR(v_true,v0);
rms  =   RMS(v_true,v0);
fprintf('Initial fSIM= %0.3f\n',fsim);
fprintf('Initial SSIM= %0.3f\n',ssim);
fprintf('Initial PSNR = %0.3f\n',psnr);
fprintf('Initial RMS = %0.3f\n',rms);


% Format and save to disk
set(gcf,'Position',[11 340 1362 462])
lambda = 1.2; c1 = 8; c2 = 0.2;
alpha = .8;  beta = 0.4;  mu =.5;    err_or = 3.6E-05;
mtrue_f = gradient(m_true);
mm_true=dom.vec2mat(mtrue_f);
name =1;
[Opts]  =  FWI_Set(   I, name, lambda, alpha, beta, mu,  c1, c2) ;  

% end
t=[];

%% Reconstruction via adjoint state method
fprintf('Spatial dof:%d, inverse prob. dof:%d\n',dom.N,nf*ns*nr)
tic
lam=0.1;
for i=1:10
            
     if(i==1)
       c0 = dom.mat2vec(v0);
       u_l=0*c0;
       q_l=0*c0;
        [m,out] = adjoint_state_2d(dom,freqs,sources,receivers,window_info,c_true,c0,SNR,maxit,u_l,q_l,v_true);
     else
         if mod(i,2)==0
         lam=0.05*lam;
        maxit=maxit+5;
         end
       if lam <1e-5
           lam = lam*0.5;
       end
       
        c0=dom.mat2vec(reconstructed_image);  
       [m,out] = adjoint_state_admm(dom,freqs,sources,receivers,window_info,c_true,c0,SNR,maxit,u_l,q_l,lam,v_true);
     end
       

     if ~isempty(find(m<0))==1
      m(m<0.0331)=0.0331;
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
     
     %%
     tic
     Sigma=sqrt(2);
     Opts.initial  =     reconstructed_image_s;
     Opts.org      =     mm_true;
     Opts.blockSize              =  7;
     Opts.nSig                   =  Sigma; 
     E_Out_New                   =  Exter_NSS_Main  (reconstructed_image_s, Opts, Sigma);      
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
