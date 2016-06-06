function valovna () 

  % velikost mreze
  n = 60;
  
  % pozicije
  u = zeros(n);
  
 
  % hitrosti
  v = zeros(n);
  
  % utezi
  w = zeros(n);
  
  % zacetni pogoji
  v(round(n/3):n/2 , round(n/3):n/2) = 1;
  %v(2*n/3, n/3) = -2;
  
  % hitrost valovanja
  c = 1;
  
  % faktor dusenja
  k = 0.0001;
  
  % razlike med tockami
  h = 0.2;
  
  % casovni korak
  dt = 0.001;
  
  % pogoj za rob
  edgeCondition = 0.0;

  i = 0;
  while (i < 100000)
    
    uNew = zeros(n);
    vNew = zeros(n);
    a = zeros(n);
    
    a = c^2/h^2 * ([u(2:end, :); zeros(1, n)] + [zeros(1, n); u(1:end-1, :)] + [u(:, 2:end) zeros(n, 1)] + [zeros(n, 1) u(:, 1:end-1)] - 4*u) - k*v;
    
   
    
    a(1, :) = a(end, :) = a (:, 1) = a(:, end) = edgeCondition;
    
    newV = v + a*dt;
    newU = u + v*dt;
    
    uNew(1, :) = uNew(end, :) = uNew(:, 1) = uNew(:, end) = edgeCondition;
    vNew(1, :) = vNew(end, :) = vNew(:, 1) = vNew(:, end) = edgeCondition;
  
    u = newU;
    v = newV;
    i += 1;
    
    if (mod(i, 10) == 0)
      surf(u);
      axis ([0 n 0 n -1 1]);
      caxis([-1 1]);
      pause(0.0001);
    
    endif
    
  endwhile
  
  

endfunction