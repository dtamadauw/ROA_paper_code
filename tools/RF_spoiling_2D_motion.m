function [Msig,Msig2,phase_motion]=RF_spoiling_2D_motion(T1,T2,TR,TE,flip,inc,Nex,Nf_in, varargin)

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

fl_sign = 1;
r_int = ones(1,Nex);
v_int = ones(1,Nex);
%Calculate phase shift by motion
if length(varargin) > 1
    
    flip_on = varargin{3};
    
    if flip_on
        fl_sign = -1;
    end
    
    %disp('Motion Simulation')
    motion_str = varargin{2};
    fun = motion_str.fun;
    delta_t = motion_str.delta_t;
    pixScale = motion_str.pixScale;
    dTR = TR*1e-3;
    
    for n=1:Nex
        if flip_on
            G = motion_str.G * (-1)^n;
        else
            G = motion_str.G;           
        end
        [r0, v] = fun((n-1)*TR*1e-3);
        phase_motion(n) = pixScale*(G*r0*delta_t+0.5*v*G*delta_t*delta_t);
        %r_int(n) = r0;
        %v_int(n) = v;
    end
    
end



for n=1:Nex

	A = Ate * throt(flip*pi/180,Rfph(n)*pi/180);
	B = Bte;
    
	M = A*M+B*on;

	Msig(n) = mean( squeeze(M(1,:)+i*M(2,:)) ) * exp(-i*Rfph(n)*pi/180);
    
    M=Ate2*M+Bte2*on;
    
    M = zrot(phase_motion(n)) * M;
    
	for k=1:Nf
		M(:,k) = zrot((gx(k)*(fl_sign)^n)+(gy(k)*(1)^n))*M(:,k);
    end
    
    dess_dem = Rfph(n)+n*inc;
    Msig2(n) = mean( squeeze(M(1,:)+i*M(2,:)) ) * exp(-i*(dess_dem)*pi/180);
    
    M=Atr*M+Btr*on;
    
    addp = pi;
    
    Rfph(n+1) = Rfph(n)+Rfinc;
    Rfinc = (n+1)*inc;
    
    


end

