%Butterworth Analog LPF parameters
Wc = 1.027;              %cut-off frequency
N = 20;                  %order 

%poles of Butterworth polynomial of degree N in the open CLHP
p1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
for k=0:19
    if (k<10)
        pole = Wc*(cos(pi/2 + pi/(2*N) + k*pi/N) + 1i*sin(pi/2 + pi/(2*N) + k*pi/N));
        p1(k+1) = pole;
    else
        x = 19-k;
        pole = Wc*(cos(pi/2 + pi/(2*N) + x*pi/N) - 1i*sin(pi/2 + pi/(2*N) + x*pi/N));
        p1(k+1) = pole;
    end
end

%Band Pass speifications
fp1 = 45;
fs1 = 40;
fs2 = 255;
fp2 = 250;

%Transformed Band Pass specs using Bilinear Transformation
f_samp = 600;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([], p1,Wc^N);   %TF with poles p and numerator Wc^N and no zeroes
                                                %numerator chosen to make the DC Gain = 1
                                      

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf1(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bpf(s) = vpa(analog_lpf1((s*s + W0*W0)/(B*s)));        %bandpass transformation
discrete_bpf(z) = vpa(analog_bpf((z-1)/(z+1)));
%bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bpf(s)); %numerical simplification to collect coeffs
ns = coeffs(ns,'All');
ds = coeffs(ds, 'All');
ns = ns/ds(1);
ds = ds/ds(1);
%coeffs of discrete bsf

[nz, dz] = numden(discrete_bpf(z));         %Extracting Numerator and Denominator                  
nz_bpf = coeffs(nz, 'All');
dz_bpf = coeffs(dz, 'All');
k = dz_bpf(1);
nz_bpf = nz_bpf/k;              %Normalizing
dz_bpf = dz_bpf/k;              %Normalizing
nz_bpf = double(nz_bpf);
dz_bpf = double(dz_bpf);
% fvtool(nz_bpf,dz_bpf)                                           %frequency response
% 
% %magnitude plot (not in log scale) 
% [H,f] = freqz(nz_bpf,dz_bpf, 1024*1024, 600e3);
% plot(f,abs(H))
% grid

%Butterworth Analog LPF parameters
Wc = 1.022;              %cut-off frequency
N = 24;                  %order 

%poles of Butterworth polynomial of degree N in the open CLHP
p = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
for k=0:N-1
    if (k<N/2)
        pole = Wc*(cos(pi/2 + pi/(2*N) + k*pi/N) + 1i*sin(pi/2 + pi/(2*N) + k*pi/N));
        p(k+1) = pole;
    else
        x = N-1-k;
        pole = Wc*(cos(pi/2 + pi/(2*N) + x*pi/N) - 1i*sin(pi/2 + pi/(2*N) + x*pi/N));
        p(k+1) = pole;
    end
end

%Band Stop speifications
fp1 = 75;
fs1 = 80;
fs2 = 215;
fp2 = 220;

%Transformed Band Pass specs using Bilinear Transformation
f_samp = 600;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([], p,Wc^N);   %TF with poles p and numerator Wc^N and no zeroes
                                                %numerator chosen to make the DC Gain = 1
                                      

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bsf(s) = vpa(analog_lpf((B*s)/(s*s + W0*W0)));        %bandpass transformation
discrete_bsf(z) = vpa(analog_bsf((z-1)/(z+1)));
%bilinear transformation

%coeffs of analog bsf
[ns1, ds1] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns_bsf = coeffs(ns1,'All');
ds_bsf = coeffs(ds1, 'All');
ns_bsf = ns_bsf/ds_bsf(1);
ds_bsf= ds_bsf/ds_bsf(1);

%coeffs of discrete bsf

[nz1, dz1] = numden(discrete_bsf(z));         %Extracting Numerator and Denominator                  
nz_bsf = coeffs(nz1, 'All');
dz_bsf = coeffs(dz1, 'All');
k = dz_bsf(1);
nz_bsf = nz_bsf/k;              %Normalizing
dz_bsf = dz_bsf/k;              %Normalizing
nz_bsf = double(nz_bsf);
dz_bsf = double(dz_bsf);
%fvtool(nz_bsf,dz_bsf)                                           %frequency response

%magnitude plot (not in log scale) 
% [H,f] = freqz(nz_bsf,dz_bsf, 1024*1024, 600e3);
% plot(f,abs(H))
% grid

a1 = coeffs(dz1*dz);
a2 = coeffs(nz1*nz);
k = a1(1);
a1 = double(a1/k);
a2=double(a2/k);

fvtool(a2, a1);                                           %frequency response

%magnitude plot (not in log scale) 
[H,f] = freqz(a2,a1, 1024*1024, 600e3);
plot(f,abs(H))
grid