function wl_mt = pSequence(mt,l,N)
wl_mt = 0;
for k = 0:N-1
    wl_mt = wl_mt + exp(1j*2*pi/N*k*(mt/2*(k+1)+l));
end
wl_mt = wl_mt/sqrt(N);
end
