close all
clear all

%%
load('vetor.mat');
sinal = -vetor;%INVERTIDO PQ O SINAL PPG Ã‰ NORMALMENTE ANALIZADO INVERTIDO
sinal2 = -vetor; %inverteu o sinal
L = length(sinal); %L é o tamnnho do sinal
vetor_x = [1:1:length(sinal)]'; % criou uma vetor coluna com o npumero de cada linha
vetor_x_derivada = [1:1:(length(sinal)-1)]'; % criou um vetor coluna com o nu=úmero de cada linha menos 1
Y = fft(sinal); %fez a trsnformada de fourier de todos os pontos do sinal
Yabs = abs(Y/L); %valor absoluto da FFT
Yang = angle(Y); %angulos da FFT
Fs = 128; %definiu a frequencia de amostragem em 128

P1 = Yabs(1 : L/2+1);
P1(2 : end-1) = 2*P1(2 : end-1); %fft da variável "sinal"

f = Fs*(0 : (L/2))/L;

%% FILTRAGEM
% wp_A = 0.00008*pi;  %frequÃªncia de passagem 0,01 Hz (0,01/128)*pi
% ws_A = 0.04*pi;  %frequÃªncia de rejeiÃ§Ã£o 5 Hz (5/128)*pi
% Ap_A = 0.1;  %atenuaÃ§Ã£o de passagem
% As_A = 50; %atenuaÃ§Ã£o de rejeiÃ§Ã£o
% h = filtro_passa_alta(wp_A, ws_A, Ap_A, As_A); %coeficiente do filtro passa ALTA
% Filtrado_PA = conv(sinal, h);
% Filtrado_PA = Filtrado_PA(1:length(sinal));
% 
% wp_B = 0.08*pi;  %frequÃªncia de passagem 10 Hz (10/128)*pi
% ws_B = 0.15*pi;  %frequÃªncia de rejeiÃ§Ã£o 20 Hz (20/128)*pi
% Ap_B = 0.1;  %atenuaÃ§Ã£o de passagem
% As_B = 50; %atenuaÃ§Ã£o de rejeiÃ§Ã£o
% h2 = filtro_passa_baixa(wp_B, ws_B, Ap_B, As_B); %coeficiente filtro passa BAIXA
% Filtrado_PB = conv(Filtrado_PA, h2);
% Filtrado_PB = Filtrado_PB(1:length(sinal));

h3 = filtro_PA_IIR(); %filtro IIR PA gerado com Filter Design
Filtrado_PA_2 = filter(h3, sinal2);
h4 = filtro_PB_IIR(); %filtro IIR PB gerado com Filter Design
Filtrado_PB_2 = filter(h4, Filtrado_PA_2);
%passou a variável "sinal2" pelo filtro passa alta e depois pelo filtro passa baixa



%% for para calcular a derivada
delta_t = 1/Fs;
derivada = zeros(length(Filtrado_PB_2), 1);
for i=1 : length(Filtrado_PB_2)-1
    derivada(i) = (Filtrado_PB_2(i+1)-Filtrado_PB_2(i));
end
%criou um vetor com a derivada do sinal - diferença entre o próximo ponto e
%o atual ------dddd

%% for onde identifica o inÃ­cio de cada ciclo
vetor_zero_cross = zeros(length(derivada), 1);
zero_cross_anterior = 0;
for i = 1 : length(derivada)-1
   if derivada(i+1)>0 && derivada(i)<0
      distancia_zero_cross = vetor_x(i+1) - zero_cross_anterior;
      if distancia_zero_cross > 80 %verifica se a distancia do zero-cross anterior Ã© maior do 70 amostras
        %plot(vetor_x(i+1), derivada(i+1), 'x'); 
        vetor_zero_cross(i+1) = 1;
        zero_cross_anterior = vetor_x(i+1);
      end
   end
end
%retorna um vetor com os pontos no qual teve uma distancia mair que 80
%amostras entre os picos

%% for que cria o vetor para identificar cada ciclo
vetor_ciclos = zeros(length(derivada), 1);
id_ciclo = 1;
for i = 1 : length(derivada)-1
    if vetor_zero_cross(i) == 1
        id_ciclo = id_ciclo + 1;
        vetor_ciclos(i+1) = id_ciclo;
    else
        vetor_ciclos(i+1) = id_ciclo;
    end    
end
%preenche o vetor_ciclos no inicio com 1 e toda vez que o vetor_zero_cros
%tever 1 eler soma 1 no vetor_ciclos e mantem repetidno este núermo até o
%proximo zero cross que aumenta mais 1 - 1111112222222223333333344444444


%% neste for em cada loop temos um vetor para o ciclo atual
distancia_picos = 0;
x_pico_anterior = 0;
y_pico_anterior = 0;
x_pico_atual = 0;
y_pico_atual = 0;
vetor_tempo_picos = zeros(vetor_ciclos(end), 1); %cria um vetor para armezenar os tempos de cada um dos picos encontrados, vetor com o mesmo número de linhas que o número de picos 

for i = 1 : vetor_ciclos(end)%vai repetir o mesmo núero de vezes que apareceu um pico
    vetor_id_ciclo = vetor_ciclos == i;
    y_ciclo_atual = Filtrado_PB_2(vetor_id_ciclo);
    x_ciclo_atual = vetor_x(vetor_id_ciclo);
    [x_pico_atual, y_pico_atual] = detectar_pico(x_ciclo_atual, y_ciclo_atual); %testa de o pr[óximo ponto é um queda após o pico
    %plot(x_ciclo_atual, (y_ciclo_atual+1300)*0.3); %plot de cada ciclo
    %plot(x_pico_atual, (y_pico_atual+1300)*0.3, 'o');
    
    vetor_tempo_picos(i) = (x_pico_atual - x_pico_anterior)/Fs; %salva o tempo entre o pico atual e o pico anterior
    
    x_pico_anterior = x_pico_atual;
    y_pico_anterior = y_pico_atual;
end

BPM = 60/mean(vetor_tempo_picos)

%subplot(3,1,1);
%plot(Filtrado_PA);
%xlim([450 2500]); %limita o valor do eixo X no plot
% subplot(4,1,1);
% plot(Filtrado_PA);
% xlim([2000 5000]); %limita o valor do eixo X no plot
% %hold on
% 
% subplot(4,1,2);
% %plot(vetor_x, Filtrado_PB_2);%soma e multiplicaÃ§Ã£o para deixar alinhados no plot
% %ylim([40 200]); %limita o valor do eixo Y no plot
% plot(vetor_x, Filtrado_PA_2);
% xlim([2000 5000]); %limita o valor do eixo X no plot

%hold on
subplot(2,1,1);
%plot(vetor_x, derivada);
plot(vetor_x, Filtrado_PA_2);
%ylim([40 200]); %limita o valor do eixo Y no plot
xlim([2000 3000]); %limita o valor do eixo X no plot

subplot(2,1,2);
plot(vetor_x, Filtrado_PB_2);
xlim([2000 3000]);

%yline(0);



%% FUNÃ‡Ã•ES
function [x_pico, y_pico] = detectar_pico(x_ciclo_atual, y_ciclo_atual)
    pico_detectado = 0;
    for i = 1 : length(y_ciclo_atual)-1
       if y_ciclo_atual(i+1) < y_ciclo_atual(i) && pico_detectado == 0
          x_pico = x_ciclo_atual(i);
          y_pico = y_ciclo_atual(i);
          pico_detectado = 1;
       end
    end
end

function h = filtro_passa_alta(wp, ws, Ap, As)
    %% Filtro PASSA ALTAS---------------  HZ ---------------------------------
    deltap = (10^(Ap/20)-1)/(10^(Ap/20)+1);
    deltas = (1+deltap)/(10^(As/20));
    delta = min(deltap, deltas);
    A = -20*log10(delta);
    Deltaw = ws - wp;
    omegac = (ws+wp)/2;
    numDeCoef = ceil(6.6*pi/Deltaw)+1; %ceil arredonda e somou 1 para ficar com nÃºmero Ã­mpar de coeficientes
    M = numDeCoef-1;

    alpha = (numDeCoef-1)/2; 
    n = [0:1:(numDeCoef-1)];
    m = n - alpha; 

    fc = omegac/pi;

    hd = -fc*sinc(fc*m); %passa altas
    hd(1,ceil(numDeCoef/2))=1-fc;%passa altas

    h = hd.*hamming(numDeCoef)';
end

function h2 = filtro_passa_baixa(wp, ws, Ap, As)
    deltap = (10^(Ap/20)-1)/(10^(Ap/20)+1);
    deltas = (1+deltap)/(10^(As/20));
    delta = min(deltap, deltas);
    A = -20*log10(delta);
    Deltaw = ws - wp;
    omegac = (ws+wp)/2;
    numDeCoef = ceil(6.6*pi/Deltaw)+1; %ceil arredonda e somou 1 para ficar com nÃºmero Ã­mpar de coeficientes
    M = numDeCoef-1;

    alpha = (numDeCoef-1)/2; 
    n = [0:1:(numDeCoef-1)];
    m = n - alpha; 

    fc = omegac/pi;

    hd = fc*sinc(fc*m); %passa baixas

    h2 = hd.*hamming(numDeCoef)';
end

function h3 = filtro_PA_IIR()
    Fs = 128;  % Sampling Frequency

    Fstop = 0.1;         % Stopband Frequency
    Fpass = 0.5;           % Passband Frequency
    Astop = 80;          % Stopband Attenuation (dB)
    Apass = 1;           % Passband Ripple (dB)
    match = 'stopband';  % Band to match exactly

    % Construct an FDESIGN object and call its BUTTER method.
    h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
    h3 = design(h, 'cheby2', 'MatchExactly', match);
end

function h4 = filtro_PB_IIR()
    Fs = 128;  % Sampling Frequency

    Fpass = 5;          % Passband Frequency
    Fstop = 10;          % Stopband Frequency
    Apass = 1;           % Passband Ripple (dB)
    Astop = 80;          % Stopband Attenuation (dB)
    match = 'stopband';  % Band to match exactly

    % Construct an FDESIGN object and call its BUTTER method.
    h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
    h4 = design(h, 'cheby2', 'MatchExactly', match);
end
