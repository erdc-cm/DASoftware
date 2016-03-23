%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimate Bathymetry for Rigid Lid dynamics in 1d
%Bathymetry given by a Gaussian hump in the center of the physical
%domain
%with length $\ell_b$. Thais is
%%
% h^t = h_0 - h_b exp\left[-3\frac{\left(x - L_x/2\right)^2}{\ell_b^2}\right]
%
% for x \in [0,L_x]
% 
% Note h = z_0 - z_b
%% Rigid lid dynamics:
%Discharge (volumetric discharge) and free-surface elevation are
%constant.
%   Q = v*h = Q_b. 
%So, 
% v^t = Q_b/h^t = v_b h_b/h^t(x)
% 
%Assumes the extended state is \vec x = [\vec h^T,\vec v^T]^T where
%\vec h, \vec v \in \R^{m/2}

classdef RigidLid < handle
    % Created on 21/03/2016 by mwf
    % Modified on 21/03/2016 by mwf
    properties
        xt = [];    % true state, a structure (xt.vec, xt.t),
        % xt.t must be within range [1,10].
        % Because matlab index from 1 instead of 0 as in C++
        zt = [];    % true observation with noise added
    end
    
    properties(GetAccess = public,SetAccess = private)
        tspan = []; % maximum step
        m = [];     % number of unknowns, = 2*nxc
        n = [];     % number of observation
        H = [];     % measurement equation, fixed in time
        F = [];     % state transition equation, change with time
        xt_var = 1; % true state variance
        obsvar = 1; % true observation variance
        loc = [];   % m x nd matrix, location coordinates of each state variable
        method = ''; % algorithm name
        kernel;     % for initial covariance
        seed;       % seed for generating realizations
        Lx  = 500.0; % length of domain
        ell = 50.0;  % length of Gaussian hump
        nxc = 100;   % number of grid cells 
        Qb  = 5.0*0.5;   % inflow discharge per unit width, m^2/sec
        h0  = 5.0;       % background depth
        hb  = 1.0;       % height of bump
        dx  = [];    % cell width
        obs_x = [];  % location of the observations
        obs_dof = []; % indices of the observations
    end
    
    methods
        % Constructor
        function obj = RigidLid(param)
            % Check input
            % Initialize parameters by case1 = RigidLid(param);
        % enforce even
            obj.m = floor(param.m/2)*2;
            obj.nxc= floor(obj.m/2)
            obj.n = param.n;
            obj.xt_var = param.x_std;
            obj.obsvar = param.obsstd;
            obj.tspan = 9;
            obj.method = param.method;
            obj.kernel = param.kernel;
            obj.seed = param.seed;
            %%TODO add other input options
            obj.dx = obj.Lx/obj.nxc
            % Get geometry
            getLOC(obj);
            obj.obs_x=linspace(0,obj.Lx,obj.n);
            % Construct H that is fixed in time
            obj.H = getH(obj,0);
            % Initialize state structure
            rng(100)
            obj.xt = obj.getxt();%obj.getx();
        end
        
        function [xnew,z] = simulate(obj,x)
            % 1-step simulation
            % for generating true process
            % compute true state x
            xt_tmp = obj.getxt();
            Fmtx = obj.getF(xt_tmp);
            xnew.vec = xt_tmp.vec;
            xnew.t = x.t + 1;
            % compute true observation z
            Hmtx = obj.H;
            noise = sqrt(obj.obsvar)*randn(size(obj.n,1));
            z.noisefree = Hmtx*xnew.vec;
            z.vec = z.noisefree + noise; % todo: add noise
            obj.F = Fmtx;
            obj.xt = xnew;
            obj.zt = z;
        end
        function [xnew,z] = simulate_orig(obj,x)
            % 1-step simulation
            % for generating true process
            % compute true state x
            Fmtx = obj.getF(x);
            xnew.vec = Fmtx*x.vec;
            xnew.t = x.t + 1;
            % compute true observation z
            Hmtx = obj.H;
            noise = sqrt(obj.obsvar)*randn(size(obj.n,1));
            z.noisefree = Hmtx*xnew.vec;
            z.vec = z.noisefree + noise; % todo: add noise
            obj.F = Fmtx;
            obj.xt = xnew;
            obj.zt = z;
        end
        
        function y = h(obj,x)
            %% Observation equation
            % for filtering
            % obj and x must reflect the correct static and dynamic initial condition before calling this function
            Hmtx = obj.H; % H is initialized in constructor for static case, for dynamic case use getH
            y.vec = Hmtx*x.vec;
        end
        
        function x = f(obj,x)
            %% Forecast equation
            % used in filtering
            % x is the transformed unknowns used for data assimilation
            % obj and x is the static and dynamic initial condition that is not corrected by DA
            % overload property: obj.F changes at the end of call
            if (size(x.vec,1) ~= obj.m) || (size(x.vec,2) ~= 1)
                error('x is not the right size');
            else
                % update state x and transition matrix F
                x.vec(obj.nxc+1:obj.m)=obj.Qb./max(x.vec(1:obj.nxc),1.0e-8);
                x.t = x.t + 1;
                obj.F = obj.getF(x)
            end
        end
        
        function x = getx(obj)
            % called in constructor of FW and DA for initialization
            % generate a random realization from N(h0,Q0), where Q0 is
            % generated form kernelfun and loc. See getQ for details
            Q0 = common.getQ(obj.loc,obj.kernel);
            [A,C,~] = svd(Q0);
            x.vec = zeros(obj.m,1);
            x.vec(1:obj.nxc) = obj.h0*ones(obj.nxc,1);% + A*sqrt(C)*randn(obj.nxc,1);
            x.t=-1.;
            obj.f(x);
        end
        
        function visualize(obj,fw_list,da_list)
            % Plot initial condition
            figure;
            subplot(2,1,1)
            obj.plotstate(da_list{1},fw_list{1})
            title('Initial condition at step 0');
            
            %% plot final state
            % PLOT NEEDS TO BECOME EITHER A FUNCTION OR A CLASS
            subplot(2,1,2)
            obj.plotstate(da_list{end},fw_list{end})
            title([obj.method,' mean and STD at step 10']);
            
        end
        
        function plotstate(obj,da,fw)
            % plot the true and estimated state
            % input:
            % da: data asismialtion data
            % fw: true data
            if isfield(da,'P')
                pvar = diag(da.P);
            elseif isfield(da,'variance')
                pvar = da.variance;
            else
                pvar = zeros(size(fw.xt.vec));
            end
            xc = obj.loc(:,1);
            plot(xc,fw.xt.vec(1:obj.nxc),'r',...
                 xc,da.x.vec(1:obj.nxc),'b',...
                 xc,da.x.vec(1:obj.nxc)+ 1.96*sqrt(pvar(1:obj.nxc)),'g',...
                 xc,da.x.vec(1:obj.nxc) - 1.96*sqrt(pvar(1:obj.nxc)),'g')
            hold on;
            xi = obj.obs_dof;
            plot(obj.obs_x,fw.xt.vec(xi),'*')
            legend('True','Estimated','95% interval');
        end
        
        function xt = getxt(obj)
            % return true profile for bathymetry and the
            % corresponding velocity field
            xt.vec = zeros(obj.m,1);
            ht = obj.h0 - obj.hb*exp((-3.0*(obj.loc(:,1)-0.5* ...
                                            obj.Lx).^2)./(obj.ell*obj.ell));
            vt = obj.Qb./ht;
            xt.vec(1:obj.nxc)=ht; xt.vec(obj.nxc+1:obj.m)=vt;
            xt.t=0.;
        end
    end
    
    methods(Access = public)
        function F = getF(obj,input)
            % Get transition matrix F = [I 0; Q_b/h 0]  
            %                          
            % 
            t = input.t; 
            depth_ind = [1:1:obj.nxc]; vel_ind = [obj.nxc+1:1:obj.m];
            F = zeros(obj.m,obj.m);
            idx_hh = sub2ind(size(F), depth_ind,depth_ind);
            F(idx_hh) = 1.0;
            idx_vh = sub2ind(size(F), vel_ind,depth_ind); 
            F(idx_vh) = -obj.Qb./max(input.vec(depth_ind).^2,1.0e-8);

        end
        
        function H = getH(obj,x)
            % Get measurement matrix H
            % nearest neighbor cells?
            obs_i=min(max(floor(obj.obs_x/obj.dx),1),obj.nxc);
            % only measuring velocity, stored after depth
            I = [1:1:obj.n]; J=obj.nxc+obs_i;
            obj.obs_dof = J;
            % H is fixed in time
            Hmtx = zeros(obj.n,obj.m);
            idx = sub2ind(size(Hmtx), I, J);
            Hmtx(idx) = 1;
            H = Hmtx;
        end
        
        function getLOC(obj)
            xi = linspace(0.5*obj.dx,obj.Lx-obj.dx,obj.nxc)';
            yi = ones(obj.nxc,1);
            obj.loc = [xi yi];
        end
        %         function H = getHcoordinate(obj)
        %             % get H from user input
        %         end
        
    end
    
    methods
        da = initializeSSM(obj,da,param);
        % initialize state space model(SSM) of Rigid for a given filter
        % type (KF,HiKF,CSKF)
    end
end

%TODO: consider put these common functions (generate covariance or a realization
%form a kernel) outside, and use mexBBFMM or rSVD to generate realization
% function Q0 = getQ(loc,kernelfun)
% % Each column of loc saves the coordinates of each point
% % For 2D cases, loc = [x y]
% % For 3D cases, loc = [x y z]
% % np: number of points
% % nd: number of dimension
% [np,nd] = size(loc);
% h = zeros(np,np); % seperation between two points
% for i = 1:nd
%     xi = loc(:,i);
%     [xj,xl]=meshgrid(xi(:),xi(:));
%     h = h + ((xj-xl)).^2;
% end
% h = sqrt(h);
% Q0 = kernelfun(h);
% end