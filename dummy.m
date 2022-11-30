function X = dft(x)
[N] = size(x);

n=0 : N-1;

for k =n
    for t = n
        X(k+1) = X(k+1) + x(t+1)*exp(-i*2*pi*k*t/N);
    end
end