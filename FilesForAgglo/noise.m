function xn = noise(x)
    % measurement noise acting on original n3 (no shifting and scaling and other coordinate transformations)

    % absolute noise
%     xn          =   x + max(xss)*randn(size(x))/20;
%     xn(xn<0)    =   0;

    % relative noise
    xn          =   x.*(1 + randn(size(x))/10);
    xn(xn<0)    =   0;
end

