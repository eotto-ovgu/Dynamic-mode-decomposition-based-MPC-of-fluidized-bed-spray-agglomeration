function m = bedmass(rho_p,xgrid,dx,n3)
    m = rho_p*pi/6*Moments(n3,xgrid,dx,3);
end