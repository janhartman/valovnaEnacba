function valovna_special_edition (zacetniPogoj, okolica) 
% MM 2. projekt, 2015/16
% Avtorji: Manja Cafuta, Jan Hartman, Zan Jelen, Ziga Kokelj
% Funkcija simulira valovanje z resevanjem valovne enacbe.
% Parametra sta dva: zacetni valovi in okolica (staticni deli)
% Ostale parametre se zaradi njihove številnosti nastavlja kar v datoteki.


  % velikost mreze : do 100 je hitrost delovanja dovolj dobra
  n = 60;
  
  % pozicije
  u = zeros(n);
  
  % hitrosti
  v = zeros(n);
  
  % utezi
  w = ones(n);
  
  % hitrost valovanja
  c = 2;
  
  % faktor dusenja - mora biti majhen
  k = 0.00001;
  
  % razlika med tockami
  h = 0.2;
  
  % casovni korak
  dt = 0.001;
  
  % nacin robov: 0 - fiksni, 1 - prosti
  rob = 0;

  % maksimalno stevilo iteracij
  stIteracij = 10000;
  
  % uporaba metode za resevanje: 0 - Euler, 1 - RK4, 2 - Leapfrog
  metoda = 1;
  
  % zacetni pogoj: kaksni bodo valovi na zacetku : od 0 do 4
  
  % okolica: kaksen bo nepremicen del: od 0 do 4
  
  % ali bomo merili energijo
  energija = 0;
  
  if (energija == 1)
    % matriki za shranjevanje energij
    eN = zeros(1, stIteracij);
    eK = zeros(1, stIteracij);
  endif
 
  % nastavljanje utezi
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
  elseif okolica == 5
    w(n/2-1:n/2+3,19:25) = w(n/2:n/2+2,18) = w(n/2:n/2+2,26) = w(n/2+1,17) = w(n/2+1,27) = w(n/2-2, 20:24) = w(n/2-3,21:23) = w(n/2-4, 22) = 0;
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
  
  % zid2
  if okolica == 3 
      zid1 = zeros([2,round((n/2)-5)]);   
      zid1(1:2, 1:round((n/2)-5)) = 0.5;   
  endif 
  
  % zid z dvema luknjama - dopolni
  % if ...
  
  % endif
 
  
  % Razlicni zacetni pogoji:
  % kaksno valovanje induciramo z nastavljanjem zacetnih hitrosti
  if zacetniPogoj == 0
    v(round(n/3):n/2 , round(n/3):n/2) = 2;
    
  elseif zacetniPogoj == 1
    v(5:10 , 5:10) = 2;
    v((n-10):(n-5),(n-10):(n-5)) = 2.5;
    v((n-20):(n-15),20:30) = 2;
    
  elseif zacetniPogoj == 2
    v(10:20,20:40) = 1.5;
    
  elseif zacetniPogoj == 3
    v((n-7):(n-1) , (round(n/2)-5):(round(n/2)+5)) = 2;
  
  endif
    
  % glavna zanka 
  for i = 1:stIteracij
    
    uNew = zeros(n);
    vNew = zeros(n);
    
    % Eulerjeva metoda
    if (metoda == 0)
        a = pospesek(u, v, c, h, k, n, rob);
        vNew = v + a*dt;
        uNew = u + v*dt;
        
        % implicitna Eulerjeva metoda (korektor, da simulacija ne podivja)
        a = pospesek(uNew, vNew, c, h, k, n, rob);
        vNew = v + a*dt;
        uNew = u + v*dt;
           
    % Runge-Kutta metoda 4. reda
    elseif (metoda == 1)
      k1 = dt*c^2/h^2 * ([u(2:end, :); zeros(1, n)] + [zeros(1, n); u(1:end-1, :)] + [u(:, 2:end) zeros(n, 1)] + [zeros(n, 1) u(:, 1:end-1)] - 4*u) - k*v;
      k2 = dt*c^2/h^2 * ([u(2:end, :); zeros(1, n)] + [zeros(1, n); u(1:end-1, :)] + [u(:, 2:end) zeros(n, 1)] + [zeros(n, 1) u(:, 1:end-1)] - 4*u) - k*v;
      k3 = dt*c^2/h^2 * ([u(2:end, :); zeros(1, n)] + [zeros(1, n); u(1:end-1, :)] + [u(:, 2:end) zeros(n, 1)] + [zeros(n, 1) u(:, 1:end-1)] - 4*u) - k*v;
      k4 = dt*c^2/h^2 * ([u(2:end, :); zeros(1, n)] + [zeros(1, n); u(1:end-1, :)] + [u(:, 2:end) zeros(n, 1)] + [zeros(n, 1) u(:, 1:end-1)] - 4*u) - k*v;
      
      vNew = v + (k1 + 2*k2 + 2*k3 + k4)/6;
      uNew = u + v*dt;
      
    % Leapfrog
    elseif (metoda == 2)
      a = pospesek(u, v, c, h, k, n, rob);
      if (mod(i, 2) == 0)
        vNew = v + a*dt;
        uNew = u;
      else
        uNew = u + v*dt;
        vNew = v;
      endif
      
      
    endif
    
    % robovi so fiksirani
    if (rob == 0)         
      uNew(1, :) = uNew(end, :) = uNew(:, 1) = uNew(:, end) = rob;
      vNew(1, :) = vNew(end, :) = vNew(:, 1) = vNew(:, end) = rob;                
    endif
    
    % prepisovanje starih matrik v nove in upostevanje utezi
    u = uNew;
    v = vNew;
    u = u.*w;
    
    % ne izrisemo v vsakem koraku iteracije   
    if (mod(i, 15) == 0)
      
      % nastavljanje okolice za prikaz
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
      elseif okolica == 5
        mc = m(n/2+3, 22);
        m(n/2-1:n/2+3,19:25) = m(n/2:n/2+2,18) = m(n/2:n/2+2,26) = m(n/2+1,17) = m(n/2+1,27) = m(n/2-2, 20:24) = m(n/2-3,21:23) = m(n/2-4, 22) = 0.1;
        m(n/2:n/2+1,19:25) = m(n/2+1,18) = m(n/2+1,26) = m(n/2+2,19:21) = m(n/2+2,23:25) = m(n/2-1,20:24) = m(n/2-2,21:23) = m(n/2-3,22) = 0.2;
        m(n/2+3, 22) = mc;

      endif
    
      % sam izris: uporablja se ukaz surf, nato še nastavimo osi
      surf(m);
      axis ([0 n 0 n -1 1]);
      caxis([-1 1]);
      pause(0.0001);
    
    endif
    
    if (energija == 1)
      % zapis energij
      eN(i) = nihajnaEnergija(u);
      eK(i) = kineticnaEnergija(v);
    endif
    
  endfor
  
  if (energija == 1)
    % izris grafov energij
    eN = eN(eN < 1000);
    plot(eN);
    hold on
    plot(eK, 'r');
    hold off
  endif
  
endfunction

% funkcija izracuna pospesek
function a = pospesek (u, v, c, h, k, n, rob)
    a = zeros(size(u));
    a = c^2/h^2 * ([u(2:end, :); zeros(1, n)] + [zeros(1, n); u(1:end-1, :)] + [u(:, 2:end) zeros(n, 1)] + [zeros(n, 1) u(:, 1:end-1)] - 4*u) - k*v;
    
   if (rob == 1)
    
      a(1,:) = c^2/h^2 * (u(2, :) + [0, u(1,2:end)] + [u(1,1:end-1), 0] -3*u(1,:)) - k*v(1,:);
      a(end,:) = c^2/h^2 * (u(end-1, :) + [0, u(end,2:end)] + [u(end,1:end-1), 0] -3*u(end,:)) -k*v(end,:) ;
      a(:,1) = c^2/h^2 * (u(:,2) + [0;u(2:end,1)] + [u(1:end-1,1); 0] -3*u(:,1)) -k*v(:,1) ; 
      a(:,end) = c^2/h^2 * (u(:,end-1) + [0;u(2:end,end)] + [u(1:end-1,end); 0] -3*u(:,end)) -k*v(:,end) ;
      
      a(1,1) = c^2/h^2 * (u(2,1) + u(1,2) - 2*u(1,1)) -k*v(1,1);
      a(1, end) = c^2/h^2 * (u(1,end-1) + u(2,end) - 2*u(1,end)) -k*v(1,end);
      a(end, 1) = c^2/h^2 * (u(end-1,1) + u(end,2) - 2*u(end,1)) -k*v(end,1);
      a(end, end) = c^2/h^2 * (u(end,end-1) + u(end-1,end) - 2*u(end,end)) -k*v(end,end);
   endif
    
  endfunction

% funkcija izracuna energijo nihanja
function W = nihajnaEnergija (u)
  W = sum(sum(u.^2));
endfunction

% funkcija izracuna kineticno energijo
function W = kineticnaEnergija (v)
  W = sum(sum(v.^2))/2;
endfunction