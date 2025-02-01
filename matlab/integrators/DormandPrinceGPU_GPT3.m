function [t, x] = DormandPrinceGPU_GPT3(f, t0, x0, t_final, dt_max, dt_min, tol)%#codegen


  % Numero di dimensioni dello stato
  n = numel(x0);

  % Alloca gli array per i coefficienti e gli stati intermedi
  k = gpuArray.zeros(7, n);
  x = gpuArray.zeros(7, n);

  % Inizializza il tempo e lo stato
  t = t0;
  x(1,:) = x0;

  % Fissa il passo iniziale
  dt = dt_max;

  % Loop fino al tempo finale
  while t < t_final
    % Calcola i coefficienti k1, k2, ..., k7
    k(1,:) = dt * f(t, x(1,:));
    k(2,:) = dt * f(t + dt / 5, x(1,:) + k(1,:) / 5);
    k(3,:) = dt * f(t + 3 * dt / 10, x(1,:) + 3 * k(1,:) / 40 + 9 * k(2,:) / 40);
    k(4,:) = dt * f(t + 4 * dt / 5, x(1,:) + 44 * k(1,:) / 45 - 56 * k(2,:) / 15 + 32 * k(3,:) / 9);
    k(5,:) = dt * f(t + 8 * dt / 9, x(1,:) + 19372 * k(1,:) / 6561 - 25360 * k(2,:) / 4371 + 64448 * k(3,:) / 4371 - 212 * k(4,:) / 143);
    k(6,:) = dt * f(t + dt, x(1,:) - 9017 * k(1,:) / 3168 + 355 * k(2,:) / 33 + 46732 * k(3,:) / 4371 - 19854 * k(4,:) / 4371 + 5292 * k(5,:) / 143);
    k(7,:) = dt * f(t + dt, x(1,:) + 35 * k(1,:) / 384 + 500 * k(3,:) / 1113 + 125 * k(4,:) / 192 - 2187 * k(5,:) / 6784 + 11 * k(6,:) / 84);

    % Calcola lo stato a t + dt utilizzando il metodo di RK4 di ordine 5 (k1, ..., k5) e il metodo di RK4 di ordine 4 (k6)
    x_new = x(1,:) + (5179 * k(1,:) + 7571 * k(3,:) + 393 * k(4,:) - 920 * k(5,:) + 187 * k(6,:) + k(7,:)) / 7168;

    % Calcola l'errore relativo
    error = norm(x_new - x(7,:)) / norm(x_new);

    % Se l'errore è troppo grande, riduci il passo e ripeti il calcolo
    if error > tol
      dt = dt / 2;
      continue;
    end

    % Se l'errore è accettabile, aggiorna il tempo e lo stato
    t = t + dt;
    x(1,:) = x_new;

    % Se l'errore è abbastanza piccolo, aumenta il passo
    if error < tol / 3
      dt = min(dt * 2, dt_max);
    end

    % Se il passo raggiunge il limite massimo o minimo, fissalo
    dt = max(dt, dt_min);
    dt = min(dt, dt_max);
  end

  % Copia il risultato finale in un array di output
  t = gather(t);
  x = gather(x(1,:));
end

% function dxdt = f(~, x, sigma, rho, beta)
% dxdt = [sigma * (x(2) - x(1)); x(1) * (rho - x(3)) - x(2); x(1) * x(2) - beta * x(3)];
% end