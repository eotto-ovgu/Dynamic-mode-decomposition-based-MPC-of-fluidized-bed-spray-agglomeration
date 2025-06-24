function referenceFunction = generateReferenceFunction(r0,ts,N_t)
    for k = 1:1:N_t    
        if k < N_t/3
            r(:,k) = r0;
        elseif k < 2*N_t/3
            r(:,k) = r0 + 0.05;
        else
            r(:,k) = r0 - 0.05;
        end
    end

    referenceFunction = fit(ts*(1:1:N_t)',r','linearinterp');
end