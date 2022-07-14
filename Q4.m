clear; clc; close all;

% reading audio file
[music,Fs] = audioread('Music.wav'); 
music11 = music(:,1);
music22 = music(:,2);
music1 = [transpose(music11) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
music2 = [transpose(music22) 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];


% generate window function
hi = zeros(64,1);
for i = 1:64
    hi(i) = sqrt(2)*sin((i-0.5)*pi/64);
end

% define Essentials
block_size = 32;
b = 8;
L = 4;
q = 2*L/(2^b-1);
len_music = length(music1);

signal1 = [zeros(1,32),music1(1,1:block_size*floor(len_music/block_size)),zeros(1,32)];
signal2 = [zeros(1,32),music2(1,1:block_size*floor(len_music/block_size)),zeros(1,32)];

% generate MDCT matrix
MDCT = zeros(block_size,2*block_size);
for i = 1:block_size
   for j = 1:2*block_size
      MDCT(i,j) = cos(pi/block_size*(i-1+1/2)*(j-1+1/2+block_size/2)); 
   end
end
MDCT = sqrt(2/block_size)*MDCT;
I_MDCT = transpose(MDCT);

% reconstruction audio
recon_audio1 = [];
recon_audio2 = [];
for k = 1:floor(length(signal1)/block_size)-1
    block1 = transpose(signal1(1+(k-1)*block_size:2*block_size+(k-1)*block_size));
    block2 = transpose(signal2(1+(k-1)*block_size:2*block_size+(k-1)*block_size));
    
    % using window function before MDCT
    blockh1 = zeros(64,1);
    for i = 1:64
        blockh1(i) = hi(i).*block1(i);
    end
    
    blockh2 = zeros(64,1);
    for i = 1:64
        blockh2(i) = hi(i).*block2(i);
    end
    
    y1 = MDCT * blockh1;
    z1 = round(y1/q);
    y_bar1 = z1*q;
    
    y2 = MDCT * blockh2;
    z2 = round(y2/q);
    y_bar2 = z2*q;
    
    % using window function after inverse MDCT
    y_inv1 = I_MDCT * y_bar1;
    y_invh1 = zeros(64,1);
    for i = 1:64
        y_invh1(i) = hi(i).*y_inv1(i);
    end
    
    y_inv2 = I_MDCT * y_bar2;
    y_invh2 = zeros(64,1);
    for i = 1:64
        y_invh2(i) = hi(i).*y_inv2(i);
    end
    
    w1(:,k) = y_invh1;
    if (k>1)
       w21 = w1(block_size+1:2*block_size , k-1);
       w31 = w1(1:block_size , k);
       recon_audio1 = [recon_audio1 ; (w21+w31)/2];
    end
    
    w2(:,k) = y_invh2;
    if (k>1)
       w22 = w2(block_size+1:2*block_size , k-1);
       w32 = w2(1:block_size , k);
       recon_audio2 = [recon_audio2 ; (w22+w32)/2];
    end
end

recon_audio = [recon_audio1 recon_audio2];
recon_audio((884209:884224),:) = [];
% play reconstruction audio
pause(2);
sound(recon_audio,Fs);

% save audios
% audiowrite('reconstructed music_b12_using window function.wav',recon_audio,Fs);


% calculate RMSE
RMSE = sqrt(mean(mean((music(:) - recon_audio(:)).^2)));
display(RMSE)


len_music = length(music);
% compare main signal and reconstructed signal in Time Domain
TotalTime1 = len_music/Fs;
t1 = 0:TotalTime1/len_music:TotalTime1-TotalTime1/len_music;
TotalTime2 = length(recon_audio)/Fs;
t2 = 0:TotalTime2/length(recon_audio):TotalTime2-TotalTime2/length(recon_audio);
figure('Name','Time Domain');
subplot(2,1,1)
plot(t1,music,'r','linewidth',2);
title('Main Signal'); 
xlabel('$ t $','Interpreter','latex');
subplot(2,1,2)
plot(t2,recon_audio,'b','linewidth',2);
title('Reconstructed Signal'); 
xlabel('$ t $','Interpreter','latex');

% compare main signal and reconstructed signal in Frequency Domain
w1=-Fs/2:1/TotalTime1:Fs/2-1/TotalTime1;
w2=-Fs/2:1/TotalTime2:Fs/2-1/TotalTime2;

figure('Name','Frequency Domain');
subplot(2,1,1)
plot(w1,abs(fftshift(fft(music))),'m','linewidth',2);
title('|fft(Main Signal)|'); 
xlabel('$ \omega $','Interpreter','latex');
subplot(2,1,2)
plot(w2,abs(fftshift(fft(recon_audio))),'c','linewidth',2);
title('|fft(Reconstructed Signal)|'); 
xlabel('$ \omega $','Interpreter','latex');




