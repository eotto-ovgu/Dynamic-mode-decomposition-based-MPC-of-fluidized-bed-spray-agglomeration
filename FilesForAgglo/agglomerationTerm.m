function dn_agg = agglomerationTerm(n,pc,p)
    % n - number density distribution [1/mm]
    % pc - process condition variable (e.g. gas inlet temperatur)
    % p - parameter structure

    N           = n.*p.dv;

    % compute beta_0 and b from process condition
    T           = pc;

    [beta_0, b] = kernelParams(T,45,250);

    % compute kernel from kernel parameters
    beta        = beta_fun(1e-9*p.v,1e-9*p.v',2,0,b);
    beta        = beta.*p.aggLim;
    
    % compute birth and death term
    B = zeros(1,p.N);
    D = zeros(1,p.N);
    for i = 1:p.N
        B(i)    = 0.5*N*(beta.*squeeze(p.Q(i,:,:)).*squeeze(p.S(i,:,:)))*N';
        D(i)    = N(i)*beta(i,:)*N';
    end
    dn_agg      = beta_0*(B - D)./p.dv;
    
end