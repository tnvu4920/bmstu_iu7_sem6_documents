function lab1()
    clear all;
    clc
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
     X = [X zeros(1,200)+15];
    X = sort(X);
    
    Mmax = max(X);
    Mmin = min(X);
    
    fprintf('Mmin = %s\n', num2str(Mmin));
    fprintf('Mmax = %s\n', num2str(Mmax));
    
    R = Mmax - Mmin;
    fprintf('R = %s\n', num2str(R));
    
    MX = getMU(X);
    fprintf('MX = %s\n', num2str(MX));
    
    Ssqr = getSsqr(X);
    fprintf('S^2 = %s\n', num2str(Ssqr));
    
    %m = getNumberOfIntervals(X);
    m = floor(log2(length(X)) + 2);
    fprintf('m = %s\n', num2str(m))
    
    %createGroup(X);
    %hold on;
    %distributionDensity(X, MX, Ssqr, m);

    figure;
    empiricF(X);
    hold on;
    distribution(X, MX, Ssqr, m);
end

function mu = getMU(X)
    n = length(X);
    mu = sum(X)/n;
end

function Ssqr = getSsqr(X)
    n = length(X);
    MX = getMU(X);
    Ssqr = sum((X - MX).^2) / (n-1);
end

function m = getNumberOfIntervals(X)
    m = floor(log2(length(X)) + 2);
end

function createGroup(X)
    n = length(X);
    m = getNumberOfIntervals(X);
    
    m = floor(log2(length(X)) + 2);
    
    intervals = zeros(1, m+1);
    numCount = zeros(1, m+1);
   
    MinX = 0;%min(X);
    Delta = (25 - 0) / m;
    fprintf('Delta = %s\n', num2str(Delta));
    
    for i = 0: m
        intervals(i+1) = MinX + Delta * i;
    end
    
    j = 1;
    count = 0;
    for i = 1:n
        while (and( j < m, X(i) >= intervals(j+1))) 
            j = j + 1; 
        end
        numCount(j) = numCount(j) + 1;
        count = count + 1;
    end
    
    for i = 1:m-1
        fprintf('[%5.2f; %5.2f) ', intervals (i), intervals(i+1));
    end 
    fprintf('[%5.2f, %5.2f]\n', intervals(m), intervals(m+1));
    
    for i = 1:m 
        fprintf('%8d       ', numCount(i));
    end
    fprintf('\n\n');

	graphBuf = numCount(1:m+1);
    for i = 1:m+1
        graphBuf(i) = numCount(i) / (n*Delta); 
    end
    
    stairs(intervals, graphBuf),grid;
end

function distributionDensity(X, MX, DX, m)
    R = X(end) - X(1);
    delta = R/m;
    Sigma = sqrt(DX);
    
    %Xn = (MX - R): delta/50 :(MX + R);
    Xn = min(X):delta/50:max(X);
    Y = normpdf(Xn, MX, Sigma);
    plot(Xn, Y), grid;
end

function distribution(X, MX, DX, m)
    R = X(end) - X(1);
    delta = R/m;
    
    Xn = (MX - R): delta :(MX + R);
    %Xn = min(X):delta / 50:max(X)+0.1;
    Y = 1/2 * (1 + erf((Xn - MX) / sqrt(2*DX))); 
    plot(Xn, Y, 'r'), grid;
    
    [yy, xx] = ecdf(X);
    stairs(xx, yy); grid;
    Xmin = min(X);
    Xmax = max(X);
    xxx = 0:Xmin / 50:Xmin;
    yyy = zeros(1, 51);
    stairs(xxx, yyy, 'b'); grid;
    
    xxxx = Xmax:(25 - Xmax) / 50:25;
    yyyy = zeros(1, 51)+1;
    
    
    %stairs(xxxx, yyyy, 'b'); grid;
end

function empiricF(X)  
    [yy, xx] = ecdf(X);
    
    Xmin = min(X);
    %plot(xx, yy, 'g'), grid;
    %stairs(xx, yy); grid;
    
    xxx = 0:Xmin / 50:Xmin;
    yyy = zeros(1, 51);
  
    
    %plot(xxx, yyy, 'b'), grid;
    %stairs(xxx, yyy);
end
