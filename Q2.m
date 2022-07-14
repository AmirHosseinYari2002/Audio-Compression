clear; clc; close all;

% generate pure tone
Fs = 8192;
t = 2/Fs:2/Fs:1;
music = cos(2*pi*128*t);

% play audio
sound(music,Fs);

% define Essentials
block_size = 32;
b = 8;
L = 4;
q = 2*L/(2^b-1);
len_music = length(music);

signal = [zeros(1,32),music(1,1:block_size*floor(len_music/block_size)),zeros(1,32)];

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
recon_audio = [];
for k = 1:floor(length(signal)/block_size)-1
    block = transpose(signal(1+(k-1)*block_size:2*block_size+(k-1)*block_size));
    y = MDCT * block;
    z = round(y/q);
    y_bar = z*q;
    w(:,k) = I_MDCT * y_bar;
    if (k>1)
       w2 = w(block_size+1:2*block_size , k-1);
       w3 = w(1:block_size , k);
       recon_audio = [recon_audio ; (w2+w3)/2];
    end
end

% play reconstruction audio
pause(2);
sound(recon_audio,Fs);

% save audios
% audiowrite('pure tone 4.wav',music,Fs);
% audiowrite('reconstructed pure tone 4.wav',recon_audio,Fs);


% calculate RMSE
RMSE = sqrt(mean((music(:) - recon_audio(:)).^2));
display(RMSE)


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

% plot RMSE in diffrent frequence
figure('Name','RMSE');
plot(w1,sqrt(mean((abs(fftshift(fft(music))) - abs(fftshift(fft(recon_audio)))).^2)),'g','linewidth',2);
title('RMSE'); 
xlabel('$ \omega $','Interpreter','latex');


%% plot RMSE according to frequence tone
clear; clc; close all;

RMSE = zeros(1,20);

for ft=1:20
    % generate pure tone
    Fs = 8192;
    t = 2/Fs:2/Fs:1;
    music = cos(2*pi*128*ft*t);

    % define Essentials
    block_size = 32;
    b = 8;
    L = 4;
    q = 2*L/(2^b-1);
    len_music = length(music);

    signal = [zeros(1,32),music(1,1:block_size*floor(len_music/block_size)),zeros(1,32)];

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
    recon_audio = [];
    for k = 1:floor(length(signal)/block_size)-1
        block = transpose(signal(1+(k-1)*block_size:2*block_size+(k-1)*block_size));
        y = MDCT * block;
        z = round(y/q);
        y_bar = z*q;
        w(:,k) = I_MDCT * y_bar;
        if (k>1)
           w2 = w(block_size+1:2*block_size , k-1);
           w3 = w(1:block_size , k);
           recon_audio = [recon_audio ; (w2+w3)/2];
        end
    end

    % calculate RMSE
    RMSE(ft) = sqrt(mean((music(:) - recon_audio(:)).^2));
end

figure('Name','RMSE');
plot(1:20,RMSE,'y','linewidth',2);
ylabel('$RMSE$','Interpreter','latex');
xlabel('$ frequence tone (\times 128) $','Interpreter','latex');






