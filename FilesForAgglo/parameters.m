function p = parameters()

    % diameter grid, unit mm
    xmin          =       0.1;               % [mm]
    xmax          =       2.5;                  % [mm]
    
    vmin          =       pi/6*xmin^3;        % [mm^3]
    vmax          =       pi/6*xmax^3;        % [mm^3]
    
    N             =       100;

    % logspace
    vb            =       logspace(log10(vmin),log10(vmax),N+1);
    dv            =       diff(vb);
    v             =       vb(1:end-1) + dv/2;

    xbounds       =       (6/pi*vb).^(1/3);
    dx            =       diff(xbounds);
    xgrid         =       (6/pi*v).^(1/3);

    % 
    m_bed_0       =       8;        % [kg]
    rho_p         =       2400*1e-9;  % [kg/mm^3]
    
    % external nucleation
    m_dot_enuc    =       0.25/60;    % [kg/s]
    L_enuc        =       0.2;
    sig_enuc      =       0.02;
    n_dot_enuc    =       generateFeed(m_dot_enuc,L_enuc,sig_enuc,xgrid,dx,dv,rho_p);

    % screens
    L2            =       0.5;        % [mm]
    sig2          =       0.65;      % [mm] 
    T2            =       normcdf(xgrid,L2,sig2);
    load('FilesForAgglo\separationFunction.mat')
    T2            =       separationFunction(xgrid/1000)';

    % initial condition
    xinit         =       0.2;
    siginit       =       0.03;
    
    n0_x          =       normpdf(xgrid,xinit,siginit);
    n0            =       n0_x.*dx./dv;
    n0            =       m_bed_0/rho_p*n0/((v.*dv)*(n0'));

    % new finite volume scheme
%     beta_0        =       1e-11;
%     a             =       0;
%     b             =       0;
%     beta          =       beta_fun(v,v',2,a,b);
%                     
%     for i = 1:N
%         for j = 1:N
%             if v(i) + v(j) > vb(end)
%                 beta(i,j) = 0;
%             end
%         end
%     end
%     vnet            =       v + v';
% 
%     BSQsym = zeros(N*N,N);
%     BSQcol = zeros(N*N,N);
%     
%     for i = 1:N
%         Q = double(and(vnet>vb(i),vnet<vb(i+1)));
%         S = vnet/v(i);
%         BSQ = (beta.*S.*Q);
% 
%         BSQsym((i-1)*N+1:i*N,:) = BSQ + BSQ';
%         BSQcol((i-1)*N+1:i*N,:) = BSQ;
%     end
    
    % assemble in structure
    p.N = N;
    p.vb = vb;
    p.dv = dv;
    p.v = v;
    p.xbounds = xbounds;
    p.dx = dx;
    p.xgrid = xgrid;
    p.m_bed_0 = m_bed_0;
    p.rho_p = rho_p;
    p.n_dot_enuc = n_dot_enuc;
%     p.beta_0 = beta_0;
    p.n0 = n0;
    p.n30 = p.xgrid.^3.*p.n0.*p.dv./p.dx;
%     p.a = a;
%     p.b = b;
%     p.beta = beta;
%     p.BSQcol = sparse(BSQcol);
%     p.BSQsym = sparse(BSQsym);
    p.T2            = T2;

    p.vnet          = p.v + p.v';

    for i = 1:p.N
        for j = 1:p.N
            if p.v(i) + p.v(j) > p.vb(end)
                p.aggLim(i,j) = 0;
            else
                p.aggLim(i,j) = 1;
            end
        end
    end

    for i = 1:p.N
        for j = 1:p.N
            for k = 1:p.N
                if (p.vb(i) < p.v(j) + p.v(k)) && (p.v(j) + p.v(k) <= p.vb(i+1))
                    p.Q(i,j,k) = 1;
                else
                    p.Q(i,j,k) = 0;
                end

                p.S(i,j,k) = (p.v(j) + p.v(k))/p.v(i);
            end
        end
    end


end
    