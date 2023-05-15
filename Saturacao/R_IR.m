%PRIMEIRA COLUNA DOS DADOS É O NÚMERO DA AMOSTRA
%SEGUNDA COLUNA DOS DADOS SÃO AS AMOSTRAS COM A COR VERMELHA
%TERCEIRA COLUNA DOS DADOS SÃO AS AMOSTRAS COM A COR INFRAVERMELHA
clear all
load('dados_R_IR_5.mat');


vetor_vermelho = dados(:,2);
vetor_infra = dados(:,3);
vetor_x = dados(:,1);

Fs = 128;
defasagem = 2500; % ESTE VALOR É ONDE IRÁ INICIAR TODOS OS 'FOR', PARA NÃO PEGAR A DEFASAGEM CAUSADA PELA FILTRAGEM;

%% filtragem dos sinais
vermelho_filtrado_PA = filter(filtro_PA_IIR(), vetor_vermelho); %filtragem sinal vermelho
vermelho_filtrado_PB = filter(filtro_PB_IIR(), vermelho_filtrado_PA);

infra_filtrado_PA = filter(filtro_PA_IIR(), vetor_infra); %filtragem sinal infravermelho
infra_filtrado_PB = filter(filtro_PB_IIR(), infra_filtrado_PA);

%% plots
% subplot(6, 1, 1);
% plot(dados(:,1), vetor_vermelho);
% xlim([2000 8000]);
% title('Vermelho - Original');
% subplot(6, 1, 2);
% plot(dados(:,1), vermelho_filtrado_PA);
% xlim([2000 8000]);
% title('Vermelho - PA');
% subplot(6, 1, 3);
% plot(dados(:,1), vermelho_filtrado_PB);
% xlim([2000 8000]);
% title('Vermelho - PB');
% 
% subplot(6, 1, 4);
% plot(dados(:,1), vetor_infra);
% xlim([2000 8000]);
% title('Infra - Original');
% subplot(6, 1, 5);
% plot(dados(:,1), infra_filtrado_PA);
% xlim([2000 8000]);
% title('Infra - PA');
% subplot(6, 1, 6);
% plot(dados(:,1), infra_filtrado_PB);
% xlim([2000 8000]);
% title('Infra - PB');

%% ---------------------------------SINAL VERMELHO--------------------------------------

%% NORMALIZAÇÃO DO SINAL VERMELHO
valor_max_vermelho = 0;
valor_min_vermelho = 0;
valor_max_infra = 0;
valor_min_infra = 0;
%for para encontrar o maior e menor valor do sinal VERMELHO
for i=defasagem : length(vermelho_filtrado_PB)
    if(vermelho_filtrado_PB(i) > valor_max_vermelho)
       valor_max_vermelho = vermelho_filtrado_PB(i);
    end
    if(vermelho_filtrado_PB(i) < valor_min_vermelho)
       valor_min_vermelho = vermelho_filtrado_PB(i); 
    end
end
%for para encontrar o maior e menor valor do sinal INFRAVERMELHO
for i=defasagem : length(infra_filtrado_PB)
    if(infra_filtrado_PB(i) > valor_max_infra)
       valor_max_infra = infra_filtrado_PB(i);
    end
    if(infra_filtrado_PB(i) < valor_min_infra)
       valor_min_infra = infra_filtrado_PB(i); 
    end
end

if(valor_max_vermelho > valor_max_infra) %este if verifica qual dos sinais tem o maior valor para ser utilizado na normalização dos dois
   valor_max = valor_max_vermelho;
else
   valor_max = valor_max_infra;
end

if(valor_min_vermelho < valor_min_infra) %este if verifica qual dos sinais tem o menor valor para ser utilizado na normalização dos sinais
   valor_min = valor_min_vermelho;
else
   valor_min = valor_min_infra;
end

vermelho_normalizado = zeros(length(vermelho_filtrado_PB), 1);
infra_normalizado = zeros(length(infra_filtrado_PB), 1);
for i=defasagem : length(vermelho_filtrado_PB)
    vermelho_normalizado(i) = (vermelho_filtrado_PB(i) - valor_min)/(valor_max - valor_min);
    infra_normalizado(i) = (infra_filtrado_PB(i) - valor_min)/(valor_max - valor_min);
end

%% for para calcular a derivada vermelho
delta_t = 1/Fs;
derivada_vermelho = zeros(length(vermelho_normalizado), 1);
for i=defasagem : length(vermelho_normalizado)-1
    derivada_vermelho(i) = (vermelho_normalizado(i+1)-vermelho_normalizado(i));
end

plot(vetor_x, vermelho_normalizado);
%% for onde identifica o início de cada ciclo
vetor_zero_cross = zeros(length(derivada_vermelho), 1);
zero_cross_anterior = 0;
for i = defasagem : length(derivada_vermelho)-1
   if derivada_vermelho(i+1)>0 && derivada_vermelho(i)<0
      distancia_zero_cross = vetor_x(i+1) - zero_cross_anterior;
      if distancia_zero_cross > 80 %verifica se a distancia do zero-cross anterior é maior do 80 amostras
        %plot(vetor_x(i+1), derivada_vermelho(i+1), 'x'); 
        vetor_zero_cross(i+1) = 1;
        zero_cross_anterior = vetor_x(i+1);
      end
   end
end

%% for que cria o vetor para identificar cada ciclo
vetor_ciclos_vermelho = zeros(length(derivada_vermelho), 1);
id_ciclo = 1;
for i = defasagem : length(derivada_vermelho)-1
    if vetor_zero_cross(i) == 1
        id_ciclo = id_ciclo + 1;
        vetor_ciclos_vermelho(i+1) = id_ciclo;
    else
        vetor_ciclos_vermelho(i+1) = id_ciclo;
    end    
end

%% for para analisar cada ciclo individualmente
media_int_diastole = 0;
media_int_sistole = 0;
for i = 1 : vetor_ciclos_vermelho(end)
    vetor_id_ciclo_vermelho = vetor_ciclos_vermelho == i;
    
    y_ciclo_atual = vermelho_normalizado(vetor_id_ciclo_vermelho);
    x_ciclo_atual = vetor_x(vetor_id_ciclo_vermelho);
    intensidade_diastole = 0;
    intesidade_sistole = 0;
    for k=1 : length(x_ciclo_atual)
        if(y_ciclo_atual(k) > intensidade_diastole)
            intensidade_diastole = y_ciclo_atual(k);
            x_int_diastole = x_ciclo_atual(k);
        end
        intensidade_sistole = y_ciclo_atual(1); %pega o primeiro valor do vetor, onde será o valor da intensidade da sistole
        x_int_sistole = x_ciclo_atual(1);
    end
    media_int_diastole = media_int_diastole + intensidade_diastole;
    media_int_sistole = media_int_sistole + intensidade_sistole;
    %[x_pico_atual, y_pico_atual] = detectar_pico(x_ciclo_atual, y_ciclo_atual);
    plot(x_ciclo_atual, y_ciclo_atual); %plot de cada ciclo
    hold on
    plot(x_int_diastole, intensidade_diastole, 'x');
    hold on
    plot(x_int_sistole, intensidade_sistole, 'o');
    hold on
    
end
media_int_diastole = media_int_diastole/id_ciclo;
media_int_sistole = media_int_sistole/id_ciclo;
valor_AC_vermelho = media_int_diastole - media_int_sistole;

valor_DC_vermelho = 0;
for i=2000 : length(vermelho_normalizado)
    valor_DC_vermelho = valor_DC_vermelho + vermelho_normalizado(i);
end
valor_DC_vermelho = valor_DC_vermelho/(length(vermelho_normalizado)-defasagem);
%yline(valor_DC_vermelho);


%% ---------------------------------SINAL INFRAVERMELHO---------------------------------

%% NORMALIZAÇÃO DO SINAL INFRAVERMELHO
% valor_max = 0;
% valor_min = 0;
% for i=defasagem : length(infra_filtrado_PB)
%     if(infra_filtrado_PB(i) > valor_max)
%        valor_max = infra_filtrado_PB(i);
%     end
%     if(infra_filtrado_PB(i) < valor_min)
%        valor_min = infra_filtrado_PB(i); 
%     end
% end
% infra_normalizado = zeros(length(infra_filtrado_PB), 1);
% for i=defasagem : length(infra_normalizado)
%     infra_normalizado(i) = (infra_filtrado_PB(i) - valor_min)/(valor_max - valor_min);
% end

%% for para calcular a derivada infra
delta_t = 1/Fs;
derivada_infra = zeros(length(infra_normalizado), 1);
for i=defasagem : length(infra_normalizado)-1
    derivada_infra(i) = (infra_normalizado(i+1)-infra_normalizado(i));
end

plot(vetor_x, infra_normalizado);
%% for onde identifica o início de cada ciclo
vetor_zero_cross = zeros(length(derivada_infra), 1);
zero_cross_anterior = 0;
for i = defasagem : length(derivada_infra)-1
   if derivada_infra(i+1)>0 && derivada_infra(i)<0
      distancia_zero_cross = vetor_x(i+1) - zero_cross_anterior;
      if distancia_zero_cross > 80 %verifica se a distancia do zero-cross anterior é maior do 80 amostras
        %plot(vetor_x(i+1), derivada_infra(i+1), 'x'); 
        vetor_zero_cross(i+1) = 1;
        zero_cross_anterior = vetor_x(i+1);
      end
   end
end

%% for que cria o vetor para identificar cada ciclo
vetor_ciclos_infra = zeros(length(derivada_infra), 1);
id_ciclo = 1;
for i = defasagem : length(derivada_infra)-1
    if vetor_zero_cross(i) == 1
        id_ciclo = id_ciclo + 1;
        vetor_ciclos_infra(i+1) = id_ciclo;
    else
        vetor_ciclos_infra(i+1) = id_ciclo;
    end    
end

%% for para analisar cada ciclo individualmente
media_int_diastole = 0;
media_int_sistole = 0;
for i = 1 : vetor_ciclos_infra(end)
    vetor_id_ciclo_infra = vetor_ciclos_infra == i;
    
    y_ciclo_atual = infra_normalizado(vetor_id_ciclo_infra);
    x_ciclo_atual = vetor_x(vetor_id_ciclo_infra);
    intensidade_diastole = 0;
    intesidade_sistole = 0;
    for k=1 : length(x_ciclo_atual)
        if(y_ciclo_atual(k) > intensidade_diastole)
            intensidade_diastole = y_ciclo_atual(k);
            x_int_diastole = x_ciclo_atual(k);
        end
        intensidade_sistole = y_ciclo_atual(1); %pega o primeiro valor do vetor, onde será o valor da intensidade da sistole
        x_int_sistole = x_ciclo_atual(1);
    end
    media_int_diastole = media_int_diastole + intensidade_diastole;
    media_int_sistole = media_int_sistole + intensidade_sistole;
    plot(x_ciclo_atual, y_ciclo_atual); %plot de cada ciclo
    hold on
    plot(x_int_diastole, intensidade_diastole, 'x');
    hold on
    plot(x_int_sistole, intensidade_sistole, 'o');
    hold on
    
end
media_int_diastole = media_int_diastole/id_ciclo;
media_int_sistole = media_int_sistole/id_ciclo;
valor_AC_infra = media_int_diastole - media_int_sistole;

valor_DC_infra = 0;
for i=2000 : length(infra_normalizado)
    valor_DC_infra = valor_DC_infra + infra_normalizado(i);
end
valor_DC_infra = valor_DC_infra/(length(infra_normalizado)-defasagem);
%yline(valor_DC_infra);

%% CÁLCULO DO SpO2
%
% R =  AC(R)  parte AC do sinal vermelho
%     ------
%      DC(R)  parte DC do sinal vermelho
% -------------
%     AC(IR)  parte AC do sinal infravermelho
%     ------
%     DC(IR)  parte DC do sinal infravermelho

R = (valor_AC_vermelho/valor_DC_vermelho)/(valor_AC_infra/valor_DC_infra);

%SpO2 = 110 - 25*R
SpO2 = -45.060*R*R + 30.354 *R + 94.845 %formula do calculo encontrada na biblioteca do MAX30102



%% FUNÇÕES
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

function h1 = filtro_PA_IIR()
    Fs = 128;  % Sampling Frequency

    Fstop = 0.1;         % Stopband Frequency
    Fpass = 0.5;         % Passband Frequency
    Astop = 80;          % Stopband Attenuation (dB)
    Apass = 1;           % Passband Ripple (dB)
    match = 'stopband';  % Band to match exactly

    % Construct an FDESIGN object and call its CHEBY2 method.
    h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
    h1 = design(h, 'cheby2', 'MatchExactly', match);
end

function h2 = filtro_PB_IIR()
    Fs = 128;  % Sampling Frequency

    Fpass = 5;          % Passband Frequency
    Fstop = 10;          % Stopband Frequency
    Apass = 1;           % Passband Ripple (dB)
    Astop = 80;          % Stopband Attenuation (dB)
    match = 'stopband';  % Band to match exactly

    % Construct an FDESIGN object and call its CHEBY2 method.
    h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
    h2 = design(h, 'cheby2', 'MatchExactly', match);
end