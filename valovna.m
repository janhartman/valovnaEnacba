function valovna () 

  % velikost mreze
  n = 60;
  
  % pozicije
  u = zeros(n);
  
 
  % hitrosti
  v = zeros(n);
  
  % utezi
  w = ones(n);
  
  % hitrost valovanja
  c = 2;
  
  % faktor dusenja
  k = 0.001;
  
  % razlike med tockami
  h = 0.2;
  
  % casovni korak
  dt = 0.001;
  
  % nacin robov: 0 - fiksni, 1 - prosti
  edgeCondition = 1;

  % maksimalno stevilo iteracij
  stIteracij = 100000;
  
  % ali uporabimo runge-kutta metodo
  rungeKutta = 0;
  
  % zacetni pogoji
  v (n/3-3:n/3+3, 2*n/3-3:2*n/3+3) = 2;
  v (2*n/3-3:2*n/3+3, n/3-3:n/3+3) = -2;
  
  okolica = 0;
  zacetniPogoj = 999;
  
  if okolica == 1
   w(30:50, 30:50) = 0;
   w(11:15,5:15) = 0;
  elseif okolica == 2
   w(51:52, 1:n) = 0;
  elseif okolica == 3
   w(41:42, 1:round((n/2)-5)) = 0; 
   w(41:42, round((n/2)+6):n) = 0; 
  elseif okolica == 4
   w(round(n/2)-9:round(n/2)+10,round(n/2)-9:round(n/2)+10) = 0; 
  endif
  
  
  %STATICNI ELEMENTI
  
  %Otok 20X20
  if okolica == 1 ||  okolica == 4
      vOtok = zeros(20);
      vOtok(1:19, 1:19) = 0.05;
      vOtok(2:18, 2:18) = 0.11;
      vOtok(3:17, 3:17) = 0.15;
      vOtok(4:16, 4:16) = 0.19;
      vOtok(5:15, 5:15) = 0.22;
      vOtok(6:14, 6:14) = 0.23;
      vOtok(7:13, 7:13) = 0.25;
      vOtok(10:12, 10:12) = 0.28;
  endif    
  %Otok 5X11   
  if okolica == 1 || okolica == 4
      dOtok = zeros([5,11]);    
      dOtok(1:5, 1:11) = 0.05;
      dOtok(2:4, 2:9) = 0.11;
      dOtok(3:3, 4:7) = 0.12;
  endif
  %Zid
  if okolica == 2 
      zid = zeros([2,n]);    
      zid(1:2, 1:n) = 0.5; 
  endif    
  if okolica == 3 
      zid1 = zeros([2,round((n/2)-5)]);   
      zid1(1:2, 1:round((n/2)-5)) = 0.5;   
  endif 
  
 
  
  % zacetni pogoji
  if zacetniPogoj == 0
    v(round(n/3):n/2 , round(n/3):n/2) = 1;
  elseif zacetniPogoj == 1
    v(5:10 , 5:10) = 1;
    v((n-10):(n-5),(n-10):(n-5)) = 2;
    v((n-20):(n-15),20:30) = 2;
  elseif zacetniPogoj == 2
    
    v(10:20,20:40) = 1;
    
  elseif zacetniPogoj == 3
    
    v((n-7):(n-1) , (round(n/2)-5):(round(n/2)+5)) = 1;
    
  endif
  
 
    
  for i = 0:stIteracij
    
    uNew = zeros(n);
    vNew = zeros(n);
    a = zeros(n);
    
    if (rungeKutta == 0)
      a = c^2/h^2 * ([u(2:end, :); zeros(1, n)] + [zeros(1, n); u(1:end-1, :)] + [u(:, 2:end) zeros(n, 1)] + [zeros(n, 1) u(:, 1:end-1)] - 4*u) - k*v;
     
      if (edgeCondition == 0)        
        vNew = v + a*dt;
        uNew = u + v*dt;
        
        uNew(1, :) = uNew(end, :) = uNew(:, 1) = uNew(:, end) = edgeCondition;
        vNew(1, :) = vNew(end, :) = vNew(:, 1) = vNew(:, end) = edgeCondition;
      
      else
        a(1,:) = c^2/h^2 * (u(2, :) + [0, u(1,2:end)] + [u(1,1:end-1), 0] -3*u(1,:)) - k*v(1,:);
        a(end,:) = c^2/h^2 * (u(end-1, :) + [0, u(end,2:end)] + [u(end,1:end-1), 0] -3*u(end,:)) -k*v(end,:) ;
        a(:,1) = c^2/h^2 * (u(:,2) + [0;u(2:end,1)] + [u(1:end-1,1); 0] -3*u(:,1)) -k*v(:,1) ; 
        a(:,end) = c^2/h^2 * (u(:,end-1) + [0;u(2:end,end)] + [u(1:end-1,end); 0] -3*u(:,end)) -k*v(:,end) ;
        
        a(1,1) = c^2/h^2 * (u(2,1) + u(1,2) - 2*u(1,1)) -k*v(1,1);
        a(1, end) = c^2/h^2 * (u(1,end-1) + u(2,end) - 2*u(1,end)) -k*v(1,end);
        a(end, 1) = c^2/h^2 * (u(end-1,1) + u(end,2) - 2*u(end,1)) -k*v(end,1);
        a(end, end) = c^2/h^2 * (u(end,end-1) + u(end-1,end) - 2*u(end,end)) -k*v(end,end);
     
        vNew = v + a*dt;
        uNew = u + v*dt;
      endif
      
    else
      k1 = dt*c^2/h^2 * ([u(2:end, :); zeros(1, n)] + [zeros(1, n); u(1:end-1, :)] + [u(:, 2:end) zeros(n, 1)] + [zeros(n, 1) u(:, 1:end-1)] - 4*u) - k*v;
      k2 = dt*c^2/h^2 * ([u(2:end, :); zeros(1, n)] + [zeros(1, n); u(1:end-1, :)] + [u(:, 2:end) zeros(n, 1)] + [zeros(n, 1) u(:, 1:end-1)] - 4*u) - k*v;
      k3 = dt*c^2/h^2 * ([u(2:end, :); zeros(1, n)] + [zeros(1, n); u(1:end-1, :)] + [u(:, 2:end) zeros(n, 1)] + [zeros(n, 1) u(:, 1:end-1)] - 4*u) - k*v;
      k4 = dt*c^2/h^2 * ([u(2:end, :); zeros(1, n)] + [zeros(1, n); u(1:end-1, :)] + [u(:, 2:end) zeros(n, 1)] + [zeros(n, 1) u(:, 1:end-1)] - 4*u) - k*v;
      
      vNew = v + (k1 + 2*k2 + 2*k3 + k4)/6;
      uNew = u + dt*v;
      
    endif
    
    u = uNew;
    v = vNew;
    u = u.*w;
    i += 1;
    
    if (mod(i, 10) == 0)
      
      m = u;
      if okolica == 1 
        m(31:50,31:50) = vOtok(1:20,1:20);
        m(11:15,5:15) = dOtok(1:5,1:11);
        
      elseif okolica == 2 
        m(51:52,1:n) = zid(1:2,1:n);
        
      elseif okolica == 3   
        m(41:42, 1:round((n/2)-5)) = zid1(1:2,1:round((n/2)-5)); 
        m(41:42, round((n/2)+6):n) = zid1(1:2,1:round((n/2)-5));
        
      elseif okolica == 4
        m(round(n/2)-9:round(n/2)+10,round(n/2)-9:round(n/2)+10) = vOtok(1:20,1:20);
        
      endif
    
      %printf("energija: %d\n", calculateEnergy(u, v));
      surf(m);
      axis ([0 n 0 n -1.5 1.5]);
      caxis([-1 1]);
      pause(0.001);
    
    endif
    
    
  endfor
  
  

endfunction


function W = calculateEnergy (u, v)
  W = sum(sum(u.^2));
  W += sum(sum(v.^2))/2;

endfunction