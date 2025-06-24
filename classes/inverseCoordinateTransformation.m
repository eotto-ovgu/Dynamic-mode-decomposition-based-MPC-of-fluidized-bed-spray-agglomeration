function state = inverseCoordinateTransformation(N_x,N_x_transformed,N_del,xop,stateScale,inputScale,ext,X)
    % 1) reverse time delaying
    N_x_extended = N_x_transformed/(N_del+1); % dimnesion of the extended state 
    state = [];
    for i = 1:N_del
        state = [state, X((i-1)*N_x_extended+1:i*N_x_extended,1)];
    end
    
    state = [state, X(N_del*N_x_extended+1:(N_del+1)*N_x_extended,:)];

    % 2) inverse nonlinear transformation
    if ext == true
        state = state(1:N_x,:);
    end

    % 3) rescaling
    [state, ~] = rescaleData(stateScale,inputScale,state,[]);

    % 4) Deviation variable
    state = state + xop;
end