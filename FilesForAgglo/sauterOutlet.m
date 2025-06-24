function d32 = sauterOutlet(T,xgrid,dx,n3)
    n3out = T'.*n3;
    d32 = Moments(n3out,xgrid,dx,3)./Moments(n3out,xgrid,dx,2);
end