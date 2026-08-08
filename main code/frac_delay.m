function y=frac_delay(x,D)

N=length(x);

Dint=floor(D);
Dfrac=D-Dint;

L=81;
mid=(L-1)/2;

k=0:L-1;

h=sinc(k-mid-Dfrac);
h=h.*hamming(L)';
h=h/sum(h);

z=conv(x,h);

% remove FIR inherent integer group delay
z=z(mid+1:mid+N);

% add integer part of propagation delay
y=zeros(1,N);

if Dint<N
    y(Dint+1:end)=z(1:N-Dint);
end

end