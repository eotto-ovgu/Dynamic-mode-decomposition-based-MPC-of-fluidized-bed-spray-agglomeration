function dn3dt = proces(t,n3,p,inputFunction)
    n           = p.xgrid'.^(-3).*n3.*p.dx'./p.dv';
    
    % all PSD-related vectors are row vectors
    n           = n.';
    
    % external nucleation
    n_dot_enuc  = p.n_dot_enuc;   
    
    % aggregration
    T           = inputFunction(t);
    dn_agg      = agglomerationTerm(n,T,p);
    
    % outlet
    K = 8e-4*(1 + (T - 363)/50);
%     K = 8e-4;
    n_dot_out   = K*p.T2.*n;
        
    % balance
    dndt        = dn_agg - n_dot_out + n_dot_enuc;
    
    % compute volume balance
    dn3dt       = (p.xgrid.^3.*dndt.*p.dv./p.dx).';

end

