clear
close all
realpos = {'real', 'positive'};
k = sym('k', realpos);
l = [ ...
  sym('L_1', realpos); ...
  sym('L_2', realpos); ...
  sym('L_3', realpos) ...
  ];
phi21 = sym('phi_21', realpos);
phi31 = sym('phi_31', realpos);
r = sym('rho_0', realpos);
g1 = sym('g_1', realpos);
g2 = sym('g_2', realpos);
c = sym('c', realpos);
f = sym('f', realpos);
w = sym('omega', realpos);
%% TMs
L1 = straightpipe(k*l(1))
L2 = straightpipe(k*l(2))
L3 = straightpipe(k*l(3))

Sb = movingarea(g1)
Sf = movingarea(1/g1)
St = movingarea(g2)

Jb = junction(g1)  
Jf = junction(1/g1)
Jt = junction(g2)

T_back = (Sb+Jb*L1)
T = (L3*(S3+Jt*L2*(Sf+Jf*(Sb+Jb*L1))))
xoverp = pressuredeflection(T, T_back, k*c)
simplify(xoverp)
xovera = vibrationdeflection(T, T_back, k*c)
simplify(xovera)
%%

%%
povera = xovera/xoverp
% povera = simplify(povera)
povera = collect((xovera/xoverp), [k r g2])
collect(simplify(povera), g2)
%% cf smalldiaphragmtaylor with l3=0 


old = l(3);
new = 0;
xoverpsmalldiaphragm = simplify(subs(xoverp, old, new))
xoverasmalldiaphragm = simplify(subs(xovera, old, new))
poverasmalldiaphragm = simplify(subs(povera, old, new))
%% Set g2=1 to compare with threeelementtaylor


old = g2;
new = 1;
xoverpsmalldiaphragm = simplify(subs(xoverp, old, new))
xoverasmalldiaphragm = simplify(subs(xovera, old, new))
poverasmalldiaphragm = simplify(subs(povera, old, new))
%% substitute the taylor expansion

sub1 = [sin(l(1)*k), cos(l(1)*k), sin(l(2)*k), cos(l(2)*k), sin(l(3)*k), cos(l(3)*k)]
sub2 = [k*l(1), 1-(k*l(1)).^2/2, k*l(2), 1-(k*l(2)).^2/2, k*l(3), 1-(k*l(3)).^2/2]
poveraTaylor = subs(povera, sub1, sub2)
poveraTaylor = expand(poveraTaylor)
simplify(poveraTaylor)
%% prove its valid for reasonable values

fpovera = matlabFunction(povera)
fpoveraTaylor = matlabFunction(poveraTaylor)
loglog(s.f, fpovera(s.l1, s.l2, s.l3, s.g2, s.k, s.r), s.f, fpoveraTaylor(s.l1, s.l2, s.l3, s.g2, s.k, s.r), '--', 'linewidth', 2)
xlabel('Frequency [hz]')
ylabel('P_a [m/m/s^2]')
title(['Taylor approximation vs sin/cos. ' s.namestr])
grid on
semilogx(s.f, 100*abs(fpovera(s.l1, s.l2, s.l3, s.g2, s.k, s.r)-fpoveraTaylor(s.l1, s.l2, s.l3, s.g2, s.k, s.r))./fpovera(s.l1, s.l2, s.l3, s.g2, s.k, s.r), 'linewidth', 2)
grid on
xlabel('Frequency')
ylabel('percent error')
title(['Percent error is much less than 1. ' s.namestr])
%% 
%% get in terms of wavelength 

% check
subs(2*pi/k, k, 2*pi*f/c); % should be c/f
kwavelength = solve(wavelength == 2*pi/k, k)
poveraTaylor2 = subs(poveraTaylor, k, kwavelength)
poveraTaylor2 = collect(poveraTaylor2, g2)

%% nondimensionalize: divide by $\rho_{0\;} L_1$

div = r*l(1)
right = poveraTaylor2/div
right = expand(right)
poveraxTaylor = right;
old = [l(2), l(3)]
new = [phi21*l(1), phi31*l(1)]
poveraxTaylor = subs(poveraxTaylor, old, new)
old = [l(1)/wavelength]
new = [x1]
poveraxTaylor = subs(poveraxTaylor, old, new)
poveraxTaylor = collect(poveraxTaylor, g2)
%%

% a = (children(collect(collect(poveraxTaylor, g2), x1)));
poverax = povera/div
fpoverax = matlabFunction(poverax)
%% for reference, the wavelength

vpa(simplify(343*u.m/u.s/(20e3*u.Hz)),3)
%% 
%% plot it

% % fv = logspace(log10(10), log10(20e3), 30);
% fv = 10e3
% cv = 344;
% rv = 1.22
% x1v = logspace(-4, 0, 30);
% x2v = logspace(-5, 0, 35);
% mesh(g2, x1v, x2v, abs(fpoverax(g2, x1v, x2v.')))
% title('taylor expansion')
% xlabel('L_1 / \lambda')
% ylabel('L_2 / \lambda')
% zlabel('abs(P_a / L_1 \rho_0)')
% set(gca, 'xscale', 'log')
% set(gca, 'yscale', 'log')
% set(gca, 'zscale', 'log')
% set(gca, 'colorscale', 'log')
% hold on
% s.l1 = 0.003;
% s.l2 = 1e-3;
% plot3(s.l1./(s.c./s.f), s.l2./(s.c./s.f), abs(fpoverax(s.l1./(s.c./s.f), s.l2./(s.c./s.f))))
% s.l1 = 0.03;
% s.l2 = 1e-3;
% plot3(s.l1./(s.c./s.f), s.l2./(s.c./s.f), abs(fpoverax(s.l1./(s.c./s.f), s.l2./(s.c./s.f))))
% s.l1 = 0.003;
% s.l2 = 1e-2;
% plot3(s.l1./(s.c./s.f), s.l2./(s.c./s.f), abs(fpoverax(s.l1./(s.c./s.f), s.l2./(s.c./s.f))))
% hold off
% % view([-90 40])
% fname = ['two element taylor nondim']
% savefig(fname)
% exportgraphics(gcf, [fname '.pdf'])
%% do it again with the non-taylor function

% fv = logspace(log10(10), log10(20e3), 30);
fv = 10e3
cv = 344;
rv = 1.22
x1v = logspace(-5, 0, 100);

% phi21v = 1e-3
phi31v = logspace(-2, 1, 100);

l1v = s.l1;
lambdav = l1v./x1v;
kv = 2*pi./lambdav;
phi21power =3
phi21v = 1*10^(phi21power);
l2v = l1v.*phi21v ;
l3v = l1v.*phi31v;

thesurf = surf(x1v, phi31v.', abs(fpoverax(l1v, l2v, l3v.', 100, kv, rv)), 'edgecolor', 'none');
hold on
thesurf = surf(x1v, phi31v.', abs(fpoverax(l1v, l2v, l3v.', 10, kv, rv)), 'edgecolor', 'none');
% plot3(cases(1).x1, cases(1).phi.', abs(fpoverax(cases(1).l1, cases(1).l2, cases(1).k, cases(1).r)), 'linewidth', 4)
% plot3(cases(2).x1, cases(2).phi.', abs(fpoverax(cases(2).l1, cases(2).l2, cases(2).k, cases(2).r)), 'linewidth', 4)
% plot3(cases(3).x1, cases(3).phi.', abs(fpoverax(cases(3).l1, cases(3).l2, cases(3).k, cases(3).r)), 'linewidth', 4)
% plot3(cases(4).x1, cases(4).phi.', abs(fpoverax(cases(4).l1, cases(4).l2, cases(4).k, cases(4).r)), 'linewidth', 4)
hold off
zlim([1e-2, 2e4])
title(['not taylor expansion, L_2/L_1=' num2str(phi21v)])
xlabel('L_1 / \lambda')
ylabel('L_2 / L_1')
zlabel('abs(P_a / L_1 \rho_0)')
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'log')
set(gca, 'zscale', 'log')
set(gca, 'colorscale', 'log')
fname = ['seven element taylor nondim']
% savefig(fname)
% exportgraphics(gcf, [fname '.pdf'])
%% plot it like the two element taylor

l1vals = [1e-3 1e-3 2e-3 2e-3];
l2l3vals = [0.5e-3 1e-3 0.75e-3 1e-3];
for n = 1:4
  cases(n) = s;
  cases(n).l1 = l1vals(n);
  cases(n).l2 = l2l3vals(n)/2;
  cases(n).l3 = l2l3vals(n)/2;
  cases(n).x1 = cases(n).l1*cases(n).f./cases(n).c;
  cases(n).phi = cases(n).l2./cases(n).l1*ones(size(cases(n).x1));
  cases(n).lambda = cases(n).c./cases(n).f;
  cases(n).namestr = ['L_1=' num2str(1e3*cases(n).l1) 'mm L_2=' num2str(1e3*cases(n).l2) 'mm.'];
end

li1 = loglog( ...
  cases(1).f, abs(fpovera(cases(1).l1, cases(1).l2, cases(1).l3, 1, cases(1).k, cases(1).r)), ...
  cases(2).f, abs(fpovera(cases(2).l1, cases(2).l2, cases(2).l3, 1, cases(2).k, cases(2).r)), ...
  cases(3).f, abs(fpovera(cases(3).l1, cases(3).l2, cases(3).l3, 1, cases(3).k, cases(3).r)), ...
  cases(4).f, abs(fpovera(cases(4).l1, cases(4).l2, cases(4).l3, 1, cases(4).k, cases(4).r)), ...
  'linewidth', 2);
hold on
li20 = loglog( ...
  cases(1).f, abs(fpovera(cases(1).l1, cases(1).l2, cases(1).l3, 20, cases(1).k, cases(1).r)), '--', ...
  cases(2).f, abs(fpovera(cases(2).l1, cases(2).l2, cases(2).l3, 20, cases(2).k, cases(2).r)), '--', ...
  cases(3).f, abs(fpovera(cases(3).l1, cases(3).l2, cases(3).l3, 20, cases(3).k, cases(3).r)), '--', ...
  cases(4).f, abs(fpovera(cases(4).l1, cases(4).l2, cases(4).l3, 20, cases(4).k, cases(4).r)), '--', ...
  'linewidth', 2);
hold off
axis('padded')
for lino = 1:numel(li20)
  li20(lino).Color = li1(lino).Color;
end
title('solid lines: s_2=s_1.  dashed lines: s_2=s_1/20.')
fmt = '%05.2f'
legend( ...
  ['L_1= ' num2str(cases(1).l1*1e3, fmt) 'mm L_2+L_3= ' num2str((cases(1).l2+cases(1).l3)*1e3, fmt) 'mm'], ...
  ['L_1= ' num2str(cases(2).l1*1e3, fmt) 'mm L_2+L_3= ' num2str((cases(2).l2+cases(2).l3)*1e3, fmt) 'mm'], ...
  ['L_1= ' num2str(cases(3).l1*1e3, fmt) 'mm L_2+L_3= ' num2str((cases(3).l2+cases(3).l3)*1e3, fmt) 'mm'], ...
  ['L_1= ' num2str(cases(4).l1*1e3, fmt) 'mm L_2+L_3= ' num2str((cases(4).l2+cases(4).l3)*1e3, fmt) 'mm'], ...
  ['Low frequency approximation'], ...
  'location', 'best')
grid on

labels('Acceleration noise due to 1 m/s^2 acceleration', 'Frequency', 'Hz', 'Acceleration noise', 'Pa/(m/s^2)')


axis('padded')



savefig('seven element exact')
exportgraphics(gcf, 'seven element exact.pdf', 'contenttype', 'vector')
%%

clf
li1 = loglog( ...
  cases(1).f, abs(fpovera(cases(1).l1, cases(1).l2, cases(1).l3, 20, cases(1).k, cases(1).r)), ...
  cases(2).f, abs(fpovera(cases(2).l1, cases(2).l2, cases(2).l3, 20, cases(2).k, cases(2).r)), ...
  cases(3).f, abs(fpovera(cases(3).l1, cases(3).l2, cases(3).l3, 20, cases(3).k, cases(3).r)), ...
  cases(4).f, abs(fpovera(cases(4).l1, cases(4).l2, cases(4).l3, 20, cases(4).k, cases(4).r)), ...
  'linewidth', 2);
hold on
% li20 = loglog( ...
%   cases(1).f, abs(fpovera(cases(1).l1, cases(1).l2, cases(1).l3, 20, cases(1).k, cases(1).r)), '--', ...
%   cases(2).f, abs(fpovera(cases(2).l1, cases(2).l2, cases(2).l3, 20, cases(2).k, cases(2).r)), '--', ...
%   cases(3).f, abs(fpovera(cases(3).l1, cases(3).l2, cases(3).l3, 20, cases(3).k, cases(3).r)), '--', ...
%   cases(4).f, abs(fpovera(cases(4).l1, cases(4).l2, cases(4).l3, 20, cases(4).k, cases(4).r)), '--', ...
%   'linewidth', 2);
% hold off
for li20no = 1:4
  li20(li20no) = yline(anoiseestimate(cases(li20no).l1, cases(li20no).l2+cases(li20no).l2, cases(li20no).r));
end


axis('padded')
for lino = 1:numel(li20)
%   li20(lino).Color = li1(lino).Color;
  li20(lino).Color = [0.7 0.7 0.7];
  li20(lino).LineWidth = 2; 
  li20(lino).LineStyle = '--'; 
end
title('Acceleration noise when s_2/s_1=20')
fmt = '%05.2f';
leg = legend( ...
  ['L_1= ' num2str(cases(1).l1*1e3, fmt) 'mm L_2+L_3= ' num2str((cases(1).l2+cases(1).l3)*1e3, fmt) 'mm'], ...
  ['L_1= ' num2str(cases(2).l1*1e3, fmt) 'mm L_2+L_3= ' num2str((cases(2).l2+cases(2).l3)*1e3, fmt) 'mm'], ...
  ['L_1= ' num2str(cases(3).l1*1e3, fmt) 'mm L_2+L_3= ' num2str((cases(3).l2+cases(3).l3)*1e3, fmt) 'mm'], ...
  ['L_1= ' num2str(cases(4).l1*1e3, fmt) 'mm L_2+L_3= ' num2str((cases(4).l2+cases(4).l3)*1e3, fmt) 'mm'], ...
  ['Taylor approximation'], ...
  'location', 'best');

leg.Position = [0.146982142857143,0.297083337236018,0.39642856232822,0.276190468385106];
grid on



axis('padded')
labels('Pressure related acceleration sensitivity', 'Frequency', 'Hz', 'Magnitude', 'Pa/(m/s^2)')
% ylim([10e-5 10e-2])
savefig('seven element vs est')
exportgraphics(gcf, 'seven element vs est.pdf', 'contenttype', 'vector')