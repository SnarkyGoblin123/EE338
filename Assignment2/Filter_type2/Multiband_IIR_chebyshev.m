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

%Parameters for Bandpass Transformation
W0 = 0.946;
B = 3.492;

%Evaluating Frequency Response of Final Filter
syms s z;
% num = [0.0252];
% den = @(s) 1*s.^7 + 0.8108*s.^6 + 2.0787*s.^5 + 1.2296*s.^4 + 1.2427*s.^3 + 0.4618*s^2 + 0.1877*s + 0.0252
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bpf(s) = analog_lpf((s*s +W0*W0)/(B*s));     %bandpass transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));          %bilinear transformation

%coeffs of analog BPF
[ns_bpf, ds_bpf] = numden(analog_bpf(s));                   %numerical simplification
ns_bpf = sym2poly(expand(ns_bpf));                          
ds_bpf = sym2poly(expand(ds_bpf));                          %collect coeffs into matrix form
k = ds_bpf(1);    
ds_bpf = ds_bpf/k;
ns_bpf = ns_bpf/k;


%coeffs of discrete BPF
[nz_bpf, dz_bpf] = numden(discrete_bpf(z));                 %numerical simplification
nz_bpf = sym2poly(expand(nz_bpf));                          
dz_bpf = sym2poly(expand(dz_bpf));                          %collect coeffs into matrix form
k = dz_bpf(1);                                          %normalisation factor
dz_bpf = dz_bpf/k;
nz_bpf = nz_bpf/k;
dz_bpf
nz_bpf

fvtool(nz_bpf,dz_bpf)
%magnitude plot (not in log scale) 
% [H,f] = freqz(nz,dz,1024*1024, 600e3);
% plot(f,abs(H))
% plot(f,angle(H))
% grid


tol = 0.07;
%Chebyshev LPF parameters
D1 = 1/(1 - tol)^2-1;       %since delta is 0.07
D2 = 1/tol^2 - 1;   
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
w_sL = 1.108;        %Stop band of Transformed LPF
N2 = ceil(acosh(sqrt(D2/D1))/acosh(w_sL));

c10 = [512 0 -1280 0 1120 0 -400 0 50 0 -1];
c10_2 = conv(c10, c10);
f1 = D1*c10_2 + [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
sol = i*roots(f1);
sol = sol(sol<0);
r = poly(sol);
[num, den] = zp2tf([], sol, double(r(end)/sqrt(1+D1)));
%Parameters for Bandstop Transformation
W0 = 0.964;
B = 1.832;

%Evaluating Frequency Response of Final Filter
syms s z;
% num = [0.0252];
% den = @(s) 1*s.^7 + 0.8108*s.^6 + 2.0787*s.^5 + 1.2296*s.^4 + 1.2427*s.^3 + 0.4618*s^2 + 0.1877*s + 0.0252
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bsf(s) = analog_lpf((B*s)/(s*s +W0*W0));     %bandpass transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));          %bilinear transformation

%coeffs of analog BSF
[ns_bsf, ds_bsf] = numden(analog_bsf(s));                   %numerical simplification
ns_bsf = sym2poly(expand(ns_bsf));                          
ds_bsf = sym2poly(expand(ds_bsf));                          %collect coeffs into matrix form
k = ds_bsf(1);    
ds_bsf = ds_bsf/k;
ns_bsf = ns_bsf/k;

%coeffs of discrete BSF
[nz_bsf, dz_bsf] = numden(discrete_bsf(z));                 %numerical simplification
nz_bsf = sym2poly(expand(nz_bsf));                          
dz_bsf = sym2poly(expand(dz_bsf));                          %collect coeffs into matrix form
k = dz_bsf(1);                                          %normalisation factor
dz_bsf = dz_bsf/k;
nz_bsf = nz_bsf/k;
dz_bsf
nz_bsf

fvtool(nz_bsf,dz_bsf)
% %magnitude plot (not in log scale) 
% [H,f] = freqz(nz,dz,1024*1024, 600e3);
% plot(f,abs(H))
% % plot(f,angle(H))
% grid

dz_final = conv(dz_bpf, dz_bsf);
nz_final = conv(nz_bpf,nz_bsf);
% k = dz_final(1);
% dz_final = double(dz_final/k);
% nz_final=double(nz_final/k);

fvtool(nz_final, dz_final);                                           %frequency response

%magnitude plot (not in log scale) 
[H,f] = freqz(nz_final,dz_final, 1024*1024, 600e3);
plot(f,abs(H))
plot(f, angle(H))
grid