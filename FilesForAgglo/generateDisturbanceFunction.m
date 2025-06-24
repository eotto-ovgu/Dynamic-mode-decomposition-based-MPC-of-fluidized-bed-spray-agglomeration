function disturbanceFunction = generateDisturbanceFunction(d0,ts,N_t)
    for k = 1:1:N_t    
        if k < N_t/3
            d(:,k) = d0;
        elseif k < 2*N_t/3
            d(:,k) = d0 + 0.1;
        else
            d(:,k) = d0 - 0.05;
        end
    end

    disturbanceFunction = fit(ts*(1:1:N_t)',d','linearinterp');
end