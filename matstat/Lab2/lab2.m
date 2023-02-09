function lab2()
    X = [10.06,8.32,8.50,8.82,6.02,6.44,7.90,7.85,5.90,7.62,8.66,...
         6.38,7.24,8.21,6.82,7.43,6.06,8.21,9.07,5.85,6.72,8.17,...
         8.53,8.68,7.21,8.43,8.77,7.27,5.79,9.78,6.44,7.24,6.83,...
         6.61,7.58,10.15,8.82,7.87,7.35,9.60,5.82,6.65,10.15,6.92,...
         6.77,9.35,6.92,7.76,6.45,7.47,6.99,9.95,7.22,7.38,7.87,...
         6.24,8.00,8.47,7.25,7.03,7.45,6.75,7.37,7.98,9.58,8.91,...
         6.14,8.19,5.07,7.47,7.29,8.78,7.86,7.82,10.09,8.54,7.21,...
         8.57,6.67,9.82,9.26,9.69,8.39,8.26,7.44,6.58,8.45,7.49,...
         7.16,9.17,8.16,8.38,7.60,8.53,6.10,7.39,7.70,8.45,7.73,...
         9.21,8.02,7.62,6.90,9.55,5.73,7.21,6.14,7.54,9.87,8.14,...
         8.16,7.50,7.60,6.25,7.03,7.07,6.61,9.68,7.65,8.32]; 
         
    n = length(X);
    
    start = 1;

    gamma = 0.9;
    alpha = (1 - gamma) / 2;
    
    mu = mean(X);
    s2 = var(X);

    fprintf('mu^(MX) = %.2f\n', mu);
    fprintf('s2^(DX) = %.2f\n', s2);
    
    mu_up = mu + sqrt(s2 / n) * tinv((1 + gamma) / 2, n - 1);
    mu_down = mu + sqrt(s2 / n) * tinv((1 - gamma) / 2, n - 1);

    fprintf('mu up = %.2f\n', mu_up);
    fprintf('mu down = %.2f\n', mu_down);
    

    sigma2_up = (n - 1) * s2 / chi2inv((1-gamma)/2, n - 1);
    sigma2_down = s2 .* (n - 1) ./ chi2inv((1 + gamma) / 2, n - 1);

    fprintf('sigma2 up = %.2f\n', sigma2_up);
    fprintf('sigma2 down = %.2f\n', sigma2_down);
    
    N = start : n;

    M = zeros(start, length(N));
    S = zeros(start, length(N));
    
    for i=N
        M(i) = mean(X(1:i));
        S(i) = var(X(1:i));
    end

    M_up = M - sqrt(S ./ N) .* tinv(1 - alpha, N - 1);
    M_down = M + sqrt(S ./ N) .* tinv(1 - alpha, N - 1);

    S_up = S .* (N - 1) ./ chi2inv(alpha, N - 1);
    S_down = S .* (N - 1) ./ chi2inv(1 - alpha, N - 1);

    figure
    hold on;
    plot([N(10), N(end)], [mu, mu], 'm');
    plot(N(10:n), M(10:n), 'g');
    plot(N(10:n), M_up(10:n), 'b');
    plot(N(10:n), M_down(10:n), 'r');
    grid on;
    hold off;
    
    figure
    hold on;
    plot([N(1), N(end)], [s2, s2], 'm');
    plot(N(10:n), S(10:n), 'g');
    plot(N(10:n), S_down(10:n), 'b');
    plot(N(10:n), S_up(10:n), 'r');
    grid on;
    hold off;
    
    for i = 1:1
        fprintf('N1 = %.2f\n', S(i));

    
    
end


