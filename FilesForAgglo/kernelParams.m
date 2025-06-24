function [beta0, b] = kernelParams(Tg,Ms,Mf)
    % Tg - gas inlet temperature - unit: K - default 363.15 
    % Ms - spray rate - unit: g/min - default 45
    % Mf - Feed rate - unit: g/min - default: 250
    

    % model matrices
%     p = load('KernelParameter.mat');
%     Abeta = p.par.beta_q';
%     Ab = p.par.b_q';
    Ab = [-155.347571879434, -0.00925230557330011, 1.84215648984003, 0.644303231885698, 1.05767920903161e-06, -2.32859760886732e-06, -0.00483286365596097, 1.25165503574574e-05, -0.000803934807778313, -0.000598514388747507];
    Abeta = [2830.98375836467, 0.206032106168870, -34.5925832821791, -11.7068848567003, 0.000213432619726040, -9.18926476641933e-05, 0.0929025184258008, -0.000222515235568846, 0.00580700957767893, 0.0106150572243328];
    
    % compute kernel parameters
    input = [1; Mf; Ms; Tg; Mf.*Ms; Mf.*Tg; Ms.*Tg; Mf.^2; Ms.^2; Tg.^2];
    beta0 = 10^(Abeta*input);
    b = (Ab*input);
end


