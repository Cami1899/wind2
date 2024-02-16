function [X_f, Y_f] = stochasticRungeKutta2(alpha, omega, sigma, weight, T, dt, X0, Y0)
    % Parametri del metodo di Runge-Kutta stocastico di ordine 2
    b1 = 1/2;
    b2 = 1/2;
    c1 = 0;
    c2 = 1;

    % Numero totale di passi temporali
    numSteps = round(T / dt);

    % Inizializzazione delle soluzioni
    X = zeros(1, numSteps + 1);
    Y = zeros(1, numSteps + 1);

    % Assegnamento delle condizioni iniziali
    X = X0*ones(length(alpha),numSteps);
    Y = Y0*ones(length(alpha),numSteps);

for i=1:length(alpha)
    
    
    % Iterazione attraverso i passi temporali
    for n = 1:numSteps
        
     
        % Incrementi deterministici e stocastici
        K1 = -alpha(i) * X(i,n)- omega(i) * Y(i,n);
        L1 = sigma(i) * randn();
        K2 = -alpha(i) * (X(i,n) + dt * K1) - omega(i) * (Y(i,n) + dt * K1);
        L2 = sigma(i) * randn();

        % Calcolo della soluzione utilizzando il metodo di Runge-Kutta stocastico di ordine 2
        X(i,n+1) = X(i,n) + dt * (b1 * K1 + b2 * K2) + sqrt(dt) * (c1 * L1 + c2 * L2);
        Y(i,n+1) = Y(i,n) + dt * (b1 * K1 + b2 * K2);   % + sqrt(dt) * (c1 * L1 + c2 * L2)
    end
end
X_f = sum(X .* sqrt(weight), 1);
Y_f = sum(Y .* sqrt(weight), 1);
end