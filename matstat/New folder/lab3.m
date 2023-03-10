function lab3()
    T = [-5.00,-4.80,-4.60,-4.40,-4.20,-4.00,-3.80,-3.60,-3.40,-3.20,...
        -3.00,-2.80,-2.60,-2.40,-2.20,-2.00,-1.80,-1.60,-1.40,-1.20,...
        -1.00,-0.80,-0.60,-0.40,-0.20,0.00,0.20,0.40,0.60,0.80,1.00,...
        1.20,1.40,1.60,1.80,2.00,2.20,2.40,2.60,2.80,3.00,3.20,3.40,...
        3.60,3.80,4.00,4.20,4.40,4.60,4.80,5.00,5.20,5.40,5.60,5.80,...
        6.00,6.20,6.40,6.60,6.80,7.00];
    Y = [430.98,374.58,332.96,298.53,272.76,244.71,251.34,204.86,...
        177.44,162.26,131.38,108.29,40.65,92.43,79.25,41.88,85.55,...
        15.80,87.68,99.36,50.52,-21.94,-46.30,8.08,58.75,-81.45,...
        78.98,43.31,84.70,11.70,17.70,28.37,46.03,44.19,98.39,...
        83.43,82.72,118.17,68.39,149.32,156.20,192.85,191.87,...
        284.01,246.53,241.26,269.22,301.30,375.81,381.93,426.70,...
        428.06,488.57,464.02,559.36,566.84,611.49,696.44,680.70,...
        695.45,795.42];
    
    [a] = polyfit(T, Y, 2)

    Yt = a(3) + a(2) * T + a(1) * T.^2;
    
    delta = sqrt(sum((Y - Yt).^2));

    figure
    hold on;
    plot(T, Y, '.b');
    plot(T, Yt, 'g');
    xlabel("T");
    ylabel("Y");
    hold off;
    grid on;

    fprintf('theta(1) = %.2f\n', a(1));
    fprintf('theta(2) = %.2f\n', a(2));
    fprintf('theta(3) = %.2f\n', a(3));
    fprintf('delta = %.2f\n', delta);
end