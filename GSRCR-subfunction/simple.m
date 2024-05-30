function M =marmousi(dom,ctr,type,smoothness)
% MARMOUSI   Generate the Marmousi model.
%   M = MARMOUSI(DOM,CTR,TYPE,SMOOTHNESS) return a matrix M which
%   represents the Marmousi velocity model on the grid associatd with
%   domain object DOM. See:
%      http://www.caam.rice.edu/~benamou/testproblem.html
%   The remaining parameters are optional.
%
%   The contrast is given by CTR, which is 1 by default.
%   Type can be 'hard' or 'smooth' for the corresponding hard or smooth
%   Marmousi models. SMOOTHNESS will apply a Gaussian smoothing filter of
%   50 pixels for a given smoothing intensity.
%
%   See also DOMAIN, SMOOTH.

if nargin<2, ctr=1; end
if nargin<3, type='hard'; end
switch type
	case 'hard'
		load overthrust_25.dat
		M=reshape(overthrust_25/1000,201,[]);
	case 'smooth'
		load marmsmooth.dat
		M=reshape(marmsmooth,122,[]);
end

% Original Marmousi model is 122 x 384 pixels.
[Xg,Yg]=ndgrid(1:201,1:801);
[Xq,Yq]=ndgrid(linspace(1,201,dom.ny),linspace(1,801,dom.nx));
% Use nearest neighbor interpolation so there is no artificial smoothing
Mgi=griddedInterpolant(Xg,Yg,M,'nearest');
M=Mgi(Xq,Yq);
if nargin==4 && smoothness>0
	M=smooth(M,smoothness);
end
M=1+ctr*(M-min(M(:)))/(max(M(:))-min(M(:)));