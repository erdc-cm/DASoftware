function kf = initializeSSM(obj,kf,param)
% initialize state space model(SSM) of obj(a forward model object) for given KF
% type (KF,HiKF,CSKF)

% map model properties to filter
kf.m = obj.m;
kf.n = obj.n;

if isprop(kf,'P')
	kf.P = zeros(kf.m,kf.m);
end

if isprop(kf,'Q')
	kf.Q = common.getQ(obj.loc,obj.kernel);
end

if isprop(kf,'R')
	kf.R = obj.obsvar*eye(kf.n,kf.n);
end

if isprop(kf,'H')
	kf.H = obj.H;
end

if isprop(kf,'K')
	kf.K = zeros(kf.m,kf.n);
end

if isprop(kf,'variance')
	kf.variance = zeros(kf.m,1);% TODO: change to VarOfState
end

if isprop(kf,'F')
	kf.F = eye(kf.m,kf.m);
end

if isprop(kf,'x')
	kf.x = obj.getx();
end

if isprop(kf,'PHT')
	kf.PHT = zeros(obj.m,obj.n);
end

if isprop(kf,'variance_Q')
	kf.variance_Q = ones(obj.m,1);
end

if isprop(kf,'QHT')
	kf.QHT = common.getKernelMatrixProd(obj.kernel,obj.loc,obj.H',param);
end

% CSKF
if isprop(kf,'A')
	kf.A = kf.getBasis(obj.loc,obj.kernel);
end

if isprop(kf,'C')
	kf.C = zeros(kf.N,kf.N);
end

if isprop(kf,'V')
	kf.V = kf.getCompCov(obj.loc,obj.kernel,kf.A);
end

end