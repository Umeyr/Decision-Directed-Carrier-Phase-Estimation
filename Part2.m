
u = 0;

for g = 1:1000

clear;
N0 = .5;
T = 1;
M = 100;
fc = .1 * 10^3;

%information bits
a = rand(1, M);
a = 1 - 2*round(a);

%discrete time
step = .01;
t = 0:step:T-step;


%pulse shape
k  = 50 ;
%p = exp(-k*((t-0.5*T).^2));  %only in the interval [0,T]
p = ones(1,100);

%noise generation
n = randn (1, M*size(t,2) );
n = sqrt(N0/2).*n;

%randomly chosen teta
teta = .3 * 10^-14;

%generation of r(t)
r = zeros(1, M*size(t,2) );
for i = 1:M

r((i-1)*size(t,2)+1:i*size(t,2)) =  cos(2*pi*fc*(t+i-1) + teta)*a(i).*p + n((i-1)*size(t,2)+1:i*size(t,2)) ;

end

%+ n((i-1)*size(t,2)+1:(i-1)*size(t,2)+size(t,2))


%downconversion
rI = zeros(1, M*size(t,2) );
rQ = zeros(1, M*size(t,2) );

for i = 1:M

rI((i-1)*size(t,2)+1:i*size(t,2)) = cos(2*pi*fc*t) .* r((i-1)*size(t,2)+1:i*size(t,2)) ;

end

for i = 1:M

rQ((i-1)*size(t,2)+1:i*size(t,2)) = sin(2*pi*fc*t) .*r((i-1)*size(t,2)+1:i*size(t,2)) ;

end

figure(1);
plot(rI(1:1000));
figure(3);
plot(rQ(1:1000));



rI = transpose(rI);
rQ = transpose(rQ);

fc = 15;
fs = 1000;

[b,f] = butter(6,fc/(fs/2));


rI = filter(b,f,rI);
rQ = filter(b,f,rQ);




figure(2);
plot(rI(1:1000));
figure(4);
plot(rQ(1:1000));

rl = rI + 1i*rQ;

rl = transpose(rl);

%outputs yi

y = zeros(1,M);



for i = 1:M
    
   y(i) =  sum( rl((i-1)*size(t,2)+1:i*size(t,2)) .* p);
    
end


tetaML = atan( sum( imag(y).*a ) / sum( real(y).*a ) );

u = u + (teta-tetaML)^2;
end

u = u / 1000;
