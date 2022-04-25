function [iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir, Mask]=generate_brain_scalar()
load xS
true_QSM = xS(:,:,:,6);
%example set of parameters
matrix_size = [256 256 98];
voxel_size = [.9375 .9375 1.5 ];
B0_dir = [ 0 0 1];    
delta_TE = 2.6e-3;
CF = 127e6;
SNR = 50;
numecho = 11;
TE =[1:numecho]*delta_TE;

%simulated the local field
D = dipole_kernel(matrix_size, voxel_size, B0_dir);
f = ifftn(D.*fftn(true_QSM))*CF*1e-6;

% simulate the magnitude image
mag = (sum(abs(xS),4)/3)+1;
mag = mag/max(mag(:)).*Mask;


iField =zeros([matrix_size numecho]);

for j = 1:numecho
    s0 = (mag.*exp(-2*pi*1i*f*TE(j) ));
    iField(:,:,:,j) = s0+1/SNR*randn(matrix_size)+1/SNR*1i*randn(matrix_size );
end



%%%%% finishing data simulation %%%%
save cheat.mat true_QSM 
save freq.mat f
clear f mag xS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







