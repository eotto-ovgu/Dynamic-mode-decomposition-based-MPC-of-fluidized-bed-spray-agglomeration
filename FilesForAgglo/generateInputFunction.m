function inputFunction = generateInputFunction(u0,ts,N_t)
    % define piecewise constant signal with random amplitudes and random
    % time intervals as well as noise
    for k = 1:1:N_t       
        if k == 1
            intLen = 0;
            count = 0;
        end
        if count == intLen
            intLen = randi(25);
            count = 1;
            u(:,k)      = u0 + 10*(2*rand()-1.0);
            u_safe      = u(:,k);
        else
            count = count + 1;            
            u(:,k)      = u_safe + 1*(2*rand()-1.0);
        end
    end

    % generate smooth function handle for improved simulation behavior
    inputFunction = fit(ts*(1:1:N_t)',u','smoothingspline');

    % plot generated inputs
%     figure()
%     plot(ts*(1:1:N_t)',u)
%     hold on
%     plot(linspace(0,ts*N_t,ts*N_t*10)',inputFunction(linspace(0,ts*N_t,ts*N_t*10)))
end