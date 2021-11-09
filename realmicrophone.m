clear
close all
k = sym('k', {'real', 'positive'}); % the wavenumber
rho0 = sym('rho_0', {'real', 'positive'}); % density of air 
c = sym('c', {'real', 'positive'}); % speed of sound
f = sym('f', {'real', 'positive'}); % frequency 
w = sym('omega', {'real', 'positive'}); % angular frequency 

l1 = sym('L_1', {'real', 'positive'}); % the length of the back volume
l2 = sym('L_2', {'real', 'positive'}); % the length of the first part of the front volume
l3 = sym('L_3', {'real', 'positive'}); % the length of the attached tubes

s1 = sym('s_1', {'real', 'positive'}); % the area of the back volume 
sd = sym('s_d', {'real', 'positive'}); % the area of the diaphragm
s3 = sym('s_3', {'real', 'positive'}); % the area of the tube

S_e = sym('S_e'); % the electronic sensitivity 
assume(S_e ~= 0); 
%% create elemental TFs 
L1 = straightpipe(k*l1) % the uniform back volume 
L2 = straightpipe(k*l2) % the uniform lower front volume
L3 = straightpipe(k*l3) % the uniform tube

Sb = movingarea(s1/sd) % the vibration of the area between the back volume and the diaphragm 
Sf = movingarea(sd/s1) % the vibration of the area between the diaphragm and the front volume
St = movingarea(s1/s3) % the vibration of the area at the bottom of the tube 

Jb = junction(s1/sd) % the junction between the back volume and the diaphragm 
Jf = junction(sd/s1) % the junction between the diaphragm and the front volume
Jt = junction(s1/s3) % the junction at the bottom of the tube 

B = Sb + Jb*L1;
T = L3*(St + Jt*L2*(Sf + Jf*(Jb*L1 + Sb)));

S_p = -S_e*B(2,1)/(T(1,1)*2*pi*f)
S_a = -S_e/(2*pi*f)^2*(-B(2,1)*T(1,2)/T(1,1) + B(2,2) - 1)
P_a = simplify(S_a/S_p)
%% analytic values
% length constants
l1constant = 1.0e-3; 
l2constant = 0.5e-3;
l3constant = 0.5e-3;
% other constants 
cconstant = 343; 
rho0constant = 1.22
P_a_fun = matlabFunction(P_a) 
fvec = logspace(log10(20), log10(20e3), 10) 
kvec = 2*pi*fvec./343

P_a_vec = P_a_fun(l1constant, l2constant, l3constant, cconstant, fvec, kvec, rho0constant, 1, 1/20)
%% approximate values 
P_a_approx = rho0constant * (l1constant/2 + l2constant + l3constant) .* ones(size(kvec));
%% plot the comparison
semilogx(fvec, abs(P_a_vec), fvec, P_a_approx, '--')
grid on 
xlabel('Frequency') 
ylabel('|P_a|') 
legend('Exact', 'Taylor Approximation')