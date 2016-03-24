function kf = initializeSSM(obj,kf,param)
% initialize state space model(SSM) of Saetrom for given kf
% type (KF,HiKF,CSKF)

% map model properties to filter
kf.m = obj.m;
kf.n = obj.n;

if isprop(kf,'P')
        %todo add different variances, correlation lengths?        
        kf.P = zeros(kf.m,kf.m);
        depth_ind = [1:1:obj.nxc]; vel_ind = [obj.nxc+1:1:obj.m];
        %%hh
	kf.P(1:obj.nxc,1:obj.nxc) = obj.xt_var(1)*common.getQ(obj.loc, ...
                                                          obj.kernel);
        %%hv
	kf.P(1:obj.nxc,obj.nxc+1:obj.m) = -sqrt(obj.xt_var(1))*sqrt(obj.xt_var(2))*common.getQ(obj.loc,obj.kernel);
        %%vh
	kf.P(obj.nxc+1:obj.m,1:obj.nxc) = kf.P(1:obj.nxc,obj.nxc+1:obj.m);
        %%vv 
        kf.P(obj.nxc+1:obj.m,obj.nxc+1:obj.m) = obj.xt_var(2)*common.getQ(obj.loc,obj.kernel);
end

if isprop(kf,'Q')
	kf.Q = zeros(kf.m,kf.m);
        depth_ind = [1:1:obj.nxc]; vel_ind = [obj.nxc+1:1:obj.m];
        %%hh
	kf.Q(1:obj.nxc,1:obj.nxc) = obj.xt_var(1)*common.getQ(obj.loc, ...
                                                          obj.kernel);
        %%hv
	kf.Q(1:obj.nxc,obj.nxc+1:obj.m) = -sqrt(obj.xt_var(1))*sqrt(obj.xt_var(2))*common.getQ(obj.loc,obj.kernel);
        %%vh
	kf.Q(obj.nxc+1:obj.m,1:obj.nxc) = kf.Q(1:obj.nxc,obj.nxc+1:obj.m);
        %%vv 
        kf.Q(obj.nxc+1:obj.m,obj.nxc+1:obj.m) = obj.xt_var(2)*common.getQ(obj.loc,obj.kernel);
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
	kf.variance = obj.xt_var(1)*ones(kf.m,1);% TODO: change to
                                                 % VarOfState
        kf.variance(obj.nxc+1:obj.m)=obj.xt_var(2);
end

if isprop(kf,'x')
    rng(obj.seed);
	kf.x = obj.getx();
end

if isprop(kf,'PHT')
	kf.PHT = common.getKernelMatrixProd(obj.kernel,obj.loc,obj.H');
end

if isprop(kf,'variance_Q')
	kf.variance_Q = zeros(obj.m,1);
end

if isprop(kf,'QHT')
	kf.QHT = zeros(obj.m,obj.n);
end

% CSKF
if isprop(kf,'A')
	kf.A = kf.getBasis(obj.loc,obj.kernel);
end

if isprop(kf,'C')
	kf.C = kf.getCompCov(obj.loc,obj.kernel,kf.A);
end

if isprop(kf,'V')
	kf.V = zeros(size(kf.C));
end
% EKF
end