function [res, f, psdx] = get_spectrum (data,Fs)
%FS - sampling Frequency
data = data-mean(data);
L=length(data);
NFFT = 2^nextpow2(L);
Y = fft(data,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
% figure; plot(f,2*abs(Y(1:NFFT/2+1))) 
% title('Single-Sided Amplitude Spectrum of data')
% xlabel('Frequency (Hz)')
% ylabel('|data(f)|')
res=2*abs(((Y(1:NFFT/2+1))));

psdx = (1/(Fs)) * res.^2;
psdx(2:end-1) = psdx(2:end-1);