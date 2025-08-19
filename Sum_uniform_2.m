function [sum_1,sum_2,sum_3,sum_4,sum_5,sum_6,sum_7,sum_8,sum_9,sum_10,sum_11,sum_12] = Sum_uniform_2(L,L_1,alpha,z,c,beta,t)
    % Refactored: use finite for-loops with relative & absolute tolerances
    Kmax    = 500;        % maximum number of terms
    tol_rel = 1e-6;       % relative tolerance for sum convergence
    tol_abs = 1e-12;      % absolute tolerance for term magnitude

    % SUM 1
    sum_1 = 0; prev_sum = 0;
    for k = 1:Kmax
        lambda_k = pi * k / L;
        mu_k     = ((alpha^2)*(cos(beta)^2)/4 + lambda_k^2) / c;
        term     = (lambda_k / mu_k) * exp(-mu_k * t);
        sum_1    = sum_1 + (-1)^k * term * sin(lambda_k * z);
        if abs(term) < tol_abs || abs(sum_1 - prev_sum) < tol_rel * abs(sum_1)
            break;
        end
        prev_sum = sum_1;
    end

    % SUM 2
    sum_2 = 0; prev_sum = 0; flip = +1;
    for k = 1:Kmax
        lambda_k = pi * k / L;
        mu_k     = ((alpha^2)*(cos(beta)^2)/4 + lambda_k^2) / c;
        mu_k2    = (alpha^2)*(cos(beta)^2) + 4*lambda_k^2;
        term     = (lambda_k / mu_k2) * exp(-mu_k * t);
        sum_2    = sum_2 + flip * term * sin(lambda_k * z);
        if abs(term) < tol_abs || abs(sum_2 - prev_sum) < tol_rel * abs(sum_2)
            break;
        end
        prev_sum = sum_2;
        flip     = -flip;
    end

    % SUM 3
    sum_3 = 0; prev_sum = 0;
    for k = 1:Kmax
        lambda_k = pi * k / L;
        mu_k     = ((alpha^2)*(cos(beta)^2)/4 + lambda_k^2) / c;
        mu_k2    = (alpha^2)*(cos(beta)^2) + 4*lambda_k^2;
        term     = (lambda_k / mu_k2) * exp(-mu_k * t);
        sum_3    = sum_3 + term * sin(lambda_k * z) * cos(lambda_k * L_1);
        if abs(term) < tol_abs || abs(sum_3 - prev_sum) < tol_rel * abs(sum_3)
            break;
        end
        prev_sum = sum_3;
    end

    % SUM 4
    sum_4 = 0; prev_sum = 0;
    for k = 1:Kmax
        lambda_k = pi * k / L;
        mu_k     = ((alpha^2)*(cos(beta)^2)/4 + lambda_k^2) / c;
        mu_k2    = (alpha^2)*(cos(beta)^2) + 4*lambda_k^2;
        term     = exp(-mu_k * t) / mu_k2;
        sum_4    = sum_4 + term * sin(lambda_k * z) * sin(lambda_k * L_1);
        if abs(term) < tol_abs || abs(sum_4 - prev_sum) < tol_rel * abs(sum_4)
            break;
        end
        prev_sum = sum_4;
    end

    % SUM 5 (identical to sum_4)
    sum_5 = sum_4;  % reused since same form

    % SUM 6
    sum_6 = 0; prev_sum = 0;
    for k = 1:Kmax
        lambda_k = pi * k / L;
        mu_k     = ((alpha^2)*(cos(beta)^2)/4 + lambda_k^2) / c;
        mu_k2    = (alpha^2)*(cos(beta)^2) + 4*lambda_k^2;
        term     = (lambda_k / mu_k2) * exp(-mu_k * t);
        sum_6    = sum_6 + term * sin(lambda_k * z) * cos(lambda_k * L_1);
        if abs(term) < tol_abs || abs(sum_6 - prev_sum) < tol_rel * abs(sum_6)
            break;
        end
        prev_sum = sum_6;
    end

    % SUM 7
    sum_7 = 0; prev_sum = 0; flip = +1;
    for k = 1:Kmax
        lambda_k = pi * k / L;
        mu_k     = ((alpha^2)*(cos(beta)^2)/4 + lambda_k^2) / c;
        mu_k2    = (alpha^2)*(cos(beta)^2) + 4*lambda_k^2;
        term     = (lambda_k / mu_k2) * exp(-mu_k * t);
        sum_7    = sum_7 + flip * term * sin(lambda_k * z);
        if abs(term) < tol_abs || abs(sum_7 - prev_sum) < tol_rel * abs(sum_7)
            break;
        end
        prev_sum = sum_7;
        flip     = -flip;
    end

    % SUM 8
    sum_8 = 0; prev_sum = 0; flip = -1;
    for k = 1:Kmax
        lambda_k = pi * k / L;
        mu_k     = ((alpha^2)*(cos(beta)^2)/4 + lambda_k^2) / c;
        mu_k2    = (alpha^2)*(cos(beta)^2) + 4*lambda_k^2;
        term     = (lambda_k / (mu_k2^2)) * exp(-mu_k * t);
        sum_8    = sum_8 + flip * term * sin(lambda_k * z);
        if abs(term) < tol_abs || abs(sum_8 - prev_sum) < tol_rel * abs(sum_8)
            break;
        end
        prev_sum = sum_8;
        flip     = -flip;
    end

    % SUM 9
    sum_9 = 0; prev_sum = 0; flip = -1;
    for k = 1:Kmax
        lambda_k = pi * k / L;
        mu_k     = ((alpha^2)*(cos(beta)^2)/4 + lambda_k^2) / c;
        mu_k2    = (alpha^2)*(cos(beta)^2) + 4*lambda_k^2;
        term     = (lambda_k^3 / (mu_k2^2)) * exp(-mu_k * t);
        sum_9    = sum_9 + flip * term * sin(lambda_k * z);
        if abs(term) < tol_abs || abs(sum_9 - prev_sum) < tol_rel * abs(sum_9)
            break;
        end
        prev_sum = sum_9;
        flip     = -flip;
    end

    % SUM 10
    sum_10 = 0; prev_sum = 0;
    for k = 1:Kmax
        lambda_k = pi * k / L;
        mu_k     = ((alpha^2)*(cos(beta)^2)/4 + lambda_k^2) / c;
        mu_k2    = (alpha^2)*(cos(beta)^2) + 4*lambda_k^2;
        term     = (lambda_k / (mu_k2^2)) * exp(-mu_k * t);
        sum_10   = sum_10 + term * sin(lambda_k * z) * cos(lambda_k * L_1);
        if abs(term) < tol_abs || abs(sum_10 - prev_sum) < tol_rel * abs(sum_10)
            break;
        end
        prev_sum = sum_10;
    end

    % SUM 11
    sum_11 = 0; prev_sum = 0;
    for k = 1:Kmax
        lambda_k = pi * k / L;
        mu_k     = ((alpha^2)*(cos(beta)^2)/4 + lambda_k^2) / c;
        mu_k2    = (alpha^2)*(cos(beta)^2) + 4*lambda_k^2;
        term     = exp(-mu_k * t) / (mu_k2^2);
        sum_11   = sum_11 + term * sin(lambda_k * z) * sin(lambda_k * L_1);
        if abs(term) < tol_abs || abs(sum_11 - prev_sum) < tol_rel * abs(sum_11)
            break;
        end
        prev_sum = sum_11;
    end

    % SUM 12
    sum_12 = 0; prev_sum = 0;
    for k = 1:Kmax
        lambda_k = pi * k / L;
        mu_k     = ((alpha^2)*(cos(beta)^2)/4 + lambda_k^2) / c;
        mu_k2    = (alpha^2)*(cos(beta)^2) + 4*lambda_k^2;
        term     = (lambda_k^2 / (mu_k2^2)) * exp(-mu_k * t);
        sum_12   = sum_12 + term * sin(lambda_k * z) * sin(lambda_k * L_1);
        if abs(term) < tol_abs || abs(sum_12 - prev_sum) < tol_rel * abs(sum_12)
            break;
        end
        prev_sum = sum_12;
    end
end



