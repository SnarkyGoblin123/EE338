
tol = 0.07;
%Chebyshev LPF parameters
D1 = 1/(1 - tol)^2-1;     %since delta is 0.15
D2 = 1/tol^2 - 1;   
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
w_sL = 1.131;        %Stop band of Transformed LPF
N1 = ceil(acosh(sqrt(D2/D1))/acosh(w_sL));

% f = 1 + 0.1562*(256x^9 - 576x^7 + 432x^5 - 120x^3 + 9x);
c9 = [256 0 -576 0 432 0 -120 0 9 0];
c9_2 = conv(c9, c9);
f1 = D1*c9_2 + [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
sol = i*roots(f1);
sol = sol(sol<0);
r = poly(sol);
[num, den] = zp2tf([], sol, double(r(end)));
% sol
%Parameters for Bandpass Transformation
W0 = 0.946;
B = 3.492;

%Evaluating Frequency Response of Final Filter
syms s z;
% num = [0.0252];
% den = @(s) 1*s.^9 + 0.8108*s.^6 + 2.0787*s.^5 + 1.2296*s.^4 + 1.2427*s.^3 + 0.4618*s^2 + 0.1877*s + 0.0252
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bpf(s) = analog_lpf((s*s +W0*W0)/(B*s));     %bandpass transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));          %bilinear transformation

%coeffs of analog BPF
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;


%coeffs of discrete BPF
[nz, dz] = numden(discrete_bpf(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
k = dz(1);                                          %normalisation factor
dz = dz/k;
nz = nz/k;


fvtool(nz,dz)
%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 600e3);
plot(f,abs(H))
plot(f,angle(H))
grid