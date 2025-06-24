function mu = Moments(n3,xgrid,dx,i)
    mu = sum(dx'.*xgrid'.^(i-3).*n3);
end