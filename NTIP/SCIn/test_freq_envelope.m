fs=96000;
length = 10;
f0=200;
df=10;
t = transpose((1:sum(length)*fs)/fs);

x = sin(2*pi*(f0+df)*t + 2*pi*(f0+df)/fs);
y = sin(2*pi*(f0)*t + 2*pi*f0/fs);
z = [x,y];

aud1=audioplayer(z,fs);
playblocking(aud1)

alpha = z(:,1) - z(:,2);
alphaz = [alpha,alpha];

aud2=audioplayer(alphaz,fs);
playblocking(aud2)

figure; plot(1/96000:1/96000:size(alpha,1)/96000,alpha')