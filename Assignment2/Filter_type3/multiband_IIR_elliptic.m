%Bandstop Filter
del = 0.15;
omega_p = 1;
omega_s = 1.108;
B = 1.832;
w0 = 0.964;

D1 = 1/(1-del)^2 - 1;
D2 = 1/del^2 - 1;
epsilon_p = sqrt(D1);
epsilon_s = sqrt(D2);


a1= epsilon_p/epsilon_s;
a = omega_p/omega_s;
%% 

syms phi k F(phi,k) K(k) Kprime(k) 
F(phi,k) = 1./sqrt(1-(k.*sin(phi)).^2);

K(k) = int(F(phi,k),phi,0,pi/2);
Kprime(k) = subs(K(k),k,sqrt(1-k^2));

N = eval((feval(K,a)*feval(Kprime,a1))/(feval(K,a1)*feval(Kprime,a)))
N = ceil(N);

a = double(ellipdeg(N,a1));

W_snew = 1/a

l = floor(N/2);
r = mod(N,2);
%% 

i = 1:1:l;
u = double(((2*i)-1)/N)
zeta = double(cde(u,a))
h_zeroes = double((1j)./(a.*zeta));
h_zeroes = [h_zeroes,conj(h_zeroes(1:l))]

syms x y darcsn(x,y) arcsn(x,y)
darcsn(x,y) = 1./(sqrt(1-x.^2)*sqrt(1-(x*y).^2));
arcsn(x,y) = int(darcsn(x,y),x,0,x);

v = double(-1j*asne(1i/epsilon_p,a1)/N);

h_poles = 1j*double(cde((u-1j*v),a));

if(r==1)
    p0 = 1j*double(cde((1-1j*v),a));
    h_poles = [h_poles,conj(h_poles(1:l)),p0]
else
    h_poles = [h_poles,conj(h_poles(1:l))];
end

h_poles
num = 0.17647*0.85*real(poly(h_zeroes))
den = real(poly(h_poles))

syms s

analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);
analog_bsf(s) = analog_lpf(((B*s)/(s*s+w0*w0)));
[ns_bsf,ds_bsf] = numden(analog_bsf(s));
ns_bsf = sym2poly(expand(ns_bsf));
ds_bsf = sym2poly(expand(ds_bsf));
k = ds_bsf(1);
ns_bsf = ns_bsf/k
ds_bsf = ds_bsf/k

syms z
discrete_bsf(z) = analog_bsf((z-1)/(z+1));
[nz_bsf,dz_bsf] = numden(discrete_bsf(z));
nz_bsf = sym2poly(expand(nz_bsf));
dz_bsf = sym2poly(expand(dz_bsf));
k = dz_bsf(1);
nz_bsf = nz_bsf/k
dz_bsf = dz_bsf/k

fvtool(nz_bsf,dz_bsf)
fvtool(num,den)

[H,f] = freqz(nz_bsf,dz_bsf,1024*1024,600e3);
plot(f,abs(H))
%% 

%Bandpass Filter
del = 0.15;
omega_p = 1;
omega_s = 1.131;
B = 3.492;
w0 = 0.946;

D1 = 1/(1-del)^2 - 1;
D2 = 1/del^2 - 1;
epsilon_p = sqrt(D1);
epsilon_s = sqrt(D2);

a1= epsilon_p/epsilon_s;
a = omega_p/omega_s;
%% 

syms phi k F(phi,k) K(k) Kprime(k) 
F(phi,k) = 1./sqrt(1-(k.*sin(phi)).^2);

K(k) = int(F(phi,k),phi,0,pi/2);
Kprime(k) = subs(K(k),k,sqrt(1-k^2));

N = eval((feval(K,a)*feval(Kprime,a1))/(feval(K,a1)*feval(Kprime,a)))
N = ceil(N);

a = double(ellipdeg(N,a1));

W_snew = 1/a

l = floor(N/2);
r = mod(N,2);
%% 

i = 1:1:l;
u = double(((2*i)-1)/N)
zeta = double(cde(u,a))
h_zeroes = double((1j)./(a.*zeta));
h_zeroes = [h_zeroes,conj(h_zeroes(1:l))]

syms x y darcsn(x,y) arcsn(x,y)
darcsn(x,y) = 1./(sqrt(1-x.^2)*sqrt(1-(x*y).^2));
arcsn(x,y) = int(darcsn(x,y),x,0,x);

v = double(-1j*asne(1i/epsilon_p,a1)/N);

h_poles = 1j*double(cde((u-1j*v),a));

if(r==1)
    p0 = 1j*double(cde((1-1j*v),a));
    h_poles = [h_poles,conj(h_poles(1:l)),p0]
else
    h_poles = [h_poles,conj(h_poles(1:l))];
end

h_poles
num = 0.17647*0.85*real(poly(h_zeroes))
den = real(poly(h_poles))
%% 

syms s

analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);
analog_bsf(s) = analog_lpf(((s*s+w0*w0)/(B*s)));
[ns_bpf,ds_bpf] = numden(analog_bsf(s));
ns_bpf = sym2poly(expand(ns_bpf));
ds_bpf = sym2poly(expand(ds_bpf));
k = ds_bpf(1);
ns_bpf = ns_bpf/k
ds_bpf = ds_bpf/k

syms z
discrete_bsf(z) = analog_bsf((z-1)/(z+1));
[nz_bpf,dz_bpf] = numden(discrete_bsf(z));
nz_bpf = sym2poly(expand(nz_bpf));
dz_bpf = sym2poly(expand(dz_bpf));
k = dz_bpf(1);
nz_bpf = nz_bpf/k
dz_bpf = dz_bpf/k

fvtool(nz_bpf,dz_bpf)
fvtool(num,den)

[H,f] = freqz(nz_bpf,dz_bpf,1024*1024,600e3);
plot(f,abs(H))

%% 

dz_final = conv(dz_bpf, dz_bsf);
nz_final = conv(nz_bpf,nz_bsf);
% k = dz_final(1);
% dz_final = double(dz_final/k);
% nz_final=double(nz_final/k);

fvtool(nz_final, dz_final);                                           %frequency response

%magnitude plot (not in log scale) 
[H,f] = freqz(nz_final,dz_final, 1024*1024, 600e3);
plot(f,abs(H))
% plot(f, angle(H))
grid