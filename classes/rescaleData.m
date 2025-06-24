function [state, input] = rescaleData(stateScale,inputScale,scaledState,scaledInput)
    state = scaledState*stateScale;
    input = scaledInput*inputScale;
end
