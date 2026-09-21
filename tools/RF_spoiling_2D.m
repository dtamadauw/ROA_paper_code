function [Msig,Msig2,phase_motion]=RF_spoiling_2D(T1,T2,TR,TE,flip,inc,Nex,Nf_in, varargin)

if length(varargin) == 0
    phix = 2*pi;
    phiy = 2*pi;
else
    moment = varargin{1};
    phix = moment(1);
    phiy = moment(2);
end


Nf_sq = round(sqrt(Nf_in));
Nf = ((Nf_sq))^2;
mesh_temp = ((1.0:Nf_sq)/Nf_sq)-0.5;
%mesh_temp = ((1.0:Nf_sq)/Nf_sq);
[gx,gy] = meshgrid(mesh_temp, mesh_temp);
gx = reshape(gx, 1, Nf)*phix;
gy = reshape(gy, 1, Nf)*phiy;

M=zeros(3,Nf);
Msig = zeros(1,Nex);
Msig2 = zeros(1,Nex);

df = 0;

[Ate,Bte] = freeprecess(TE,T1,T2,df);
[Ate2,Bte2] = freeprecess(TR-2*TE,T1,T2,df);
[Atr,Btr] = freeprecess(TE,T1,T2,df);


M = [zeros(2,Nf);ones(1,Nf)];
on = ones(1,Nf);
phase_motion= zeros(1,Nf);
Rfph = zeros(1,Nf+1);	
Rfph(1) = 0;
Rfinc = inc;

dess_demod = 0;

fl_sign = varargin{2};


for n=1:Nex

	A = Ate * throt(flip*pi/180,Rfph(n)*pi/180);
	B = Bte;
    
	M = A*M+B*on;

	Msig(n) = mean( squeeze(M(1,:)+i*M(2,:)) ) * exp(-i*Rfph(n)*pi/180);
        
    M=Ate2*M+Bte2*on;
    
	for k=1:Nf
		M(:,k) = zrot((gx(k)*(fl_sign)^n)+(gy(k)*(1)^n))*M(:,k);
    end
    
    dess_dem = Rfph(n)+n*inc;
    Msig2(n) = mean( squeeze(M(1,:)+i*M(2,:)) ) * exp(-i*(dess_dem)*pi/180);
    

    M=Atr*M+Btr*on;
        
    Rfph(n+1) = Rfph(n)+Rfinc;
    Rfinc = (n+1)*inc;
    
end
