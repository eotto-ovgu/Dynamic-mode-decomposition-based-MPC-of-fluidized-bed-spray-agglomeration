classdef data
    properties
        xss             % steady state for linearization
        uss             % steady state input for linearization
        x               % state 
        u               % input 
        p               % parameter set 
        tspan           % time stamps of the data sequence
        ts              
        label           % label describing the data set
    end
    methods
        function obj = data() % constructor
        end
        function obj = createNewData(obj,label)
            obj.label = label;
            [obj.xss,obj.uss,obj.x,obj.u,obj.tspan,obj.ts,obj.p] = forwardSim();
        end        
        function obj = importData(obj,xss,uss,x,u,p,tspan,ts,label)
            obj.xss = xss;
            obj.uss = uss;
            obj.x = x;
            obj.u = u;
            obj.p = p;
            obj.tspan = tspan;            
            obj.ts = ts;            
            obj.label = label;            
        end
    end
end
function [xss,uss,x,u,tspan,ts,p] = forwardSim()
    % define model and simulation parameters
    p           =   parameters();
    opts        =   odeset('abstol',1e-3,'reltol',1e-3);

    % define time interval and sampling time of data
    T_end       =   6*3600;
    ts          =   1*60;
    tspan       =   linspace(0,T_end,T_end/ts+1);
    N_t         =   length(tspan)-1;

    % compute operating point
    figure(23); hold on; grid on
    count = 0;
    for i = [363]
        count = count + 1;
        uss         =   i; %593.1248e-6;
        xinit       =   p.n30';
    
        % determine steady state
        [~, x_act]  =   ode15s(@proces,tspan,xinit,opts,p,@(t) uss);
        xss         =   x_act(end,:)';
        
        plot(p.xgrid,xss)
        d32(count)  =   sauter(p.xgrid,p.dx,xss);
        d32out(count)=   sauterOutlet(p.T2,p.xgrid,p.dx,xss);
        m(count)    = bedmass(p.rho_p,p.xgrid,p.dx,xss);
    end

    
    % generate smooth input function
    inputFunction = generateInputFunction(uss,ts,N_t);

    % save inputs at discretization time steps
    u           =   inputFunction(tspan(1:end-1))';
          
    % compute states at discretization time steps
    [~, x]      =   ode15s(@proces,tspan,xss,opts,p,inputFunction);
    x           =   x';

    % add noise
    xn          =   noise(x);   
    x           =   xn;
end
