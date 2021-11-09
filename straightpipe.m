function t = straightpipe(kL, r0c)
  arguments
    kL (1,1) = sym('k', {'positive', 'real'})*sym('L', {'positive', 'real'})
    r0c (1,1) = sym('rho_0', {'positive', 'real'})*sym('c', {'positive', 'real'})
  end
  %   t = [];
  %   t(1,1) = cos(kL);
  %   t(1,2) = -1i*r0c*sin(kL);
  %   t(2,1) = -1i.*sin(kL)/r0c;
  %   t(2,2) = cos(kL);
  
  t = [cos(kL), -1i*r0c*sin(kL); -1i.*sin(kL)/r0c, cos(kL)];
end