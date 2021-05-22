%%% ///////////////////////////////////////////////////////////////////////////////////
%%%% EC1 - PEA3306 - 21/05/2021 - TRANSFORMADORES
%%%%
%%%% >>>>> 11260811 Gabriel Macias de Oliveira
%%%%       11261656 Leonardo Isao Komura
%%%%       10297265 Rodrigo Ryuji Ikegami
%%%%       11259715 Vanderson da Silva dos Santos
%%%%
%%%% FATOR DE POTENCIA DA CARGA: 0,7 CAPACITIVO
%%%%
%%%% VARIAVEIS
%%%%
%%%% SN   Potência nominal                    2.5   [MVA]
%%%% VNP  Tensão Nominal no Primário          19.1  [kV]
%%%% VNS  Tensão Nominal no Secundário        3.81  [kV]
%%%% f    Frequência de Operação              50    [Hz]
%%%% Vv   Tensão no Ensaio em Vazio           3810  [V]
%%%% Iv   Corrente no Ensaio em Vazio         9.86  [A]
%%%% Pv   Potência Ativa no Ensaio em Vazio   8140  [W]
%%%% Vcc  Tensão no Ensaio em Curto           854   [V]
%%%% Icc  Corrente no Ensaio em Curto         131   [A]
%%%% Pcc  Potência no Ensaio em Curto         8890  [W]
%%%% a    Relação de Transformação
%%%% fp   Fator de Potência da Carga          0.7 cap
%%%%
%%%% ///////////////////////////////////////////////////////////////////////////////////
clear; clc;
nusp = 11260811;
grupo = 14; % numero do grupo(1 a 60)

fp = 0.7;
SN = 2.5e6;
VNP = 19.1e3;
VNS = 3.81e3;
f = 50;
Vv = 3810;
Iv = 9.86;
Pv = 8140;
Vcc = 854;
Icc = 131;
Pcc = 8890;
a = VNP / VNS;

%%% ////////////////////////////////////////////////////
%%% CALCULO DOS PARAMETROS
%%% ////////////////////////////////////////////////

%%% Parâmetros calculados a partir do ensaio em vazio (baixa tensão)

fpv = Pv / (Vv * Iv);                  % Fator de potência no ensaio em vazio (fp = P/|S| = P/(|V|*|I|))
Ipv = Iv * fpv;                        % Corrente de perdas no ferro no ensaio em vazio [A]
Imv = Iv * sqrt(1 - fpv^2);            % Corrente de magnetização no ensaio em vazio [A]

rp = Vv/Ipv                            % Resistência de perdas por histerese e correntes de Foucault no núcleo [Ohm]
xm = Vv/Imv                            % Reatância de magnetização no núcleo. [Ohm]

%%% Parâmetros calculados a partir do ensaio em curto (alta tensão)

Zcc = Vcc/Icc;                         % Impedância de curto circuito [Ohm]
rcc = Pcc/(Icc^2);                     % Resistência de curto circuito [Ohm]
xcc = sqrt(Zcc^2 - rcc^2);             % Reatância de curto circuito [Ohm]

r1 = rcc / 2                           % Resistência do enrolamento primário [Ohm]
x1 = xcc / 2                           % Reatância de dispersão de fluxo no primário [Ohm]

%%% ////////////////////////////////////////////////////
%%% PARAMETROS NA AT
%%% ////////////////////////////////////////////////

rp_at = rp * a^2                       % rp refletido ao lado de alta tensão [Ohm]
xm_at = xm * a^2                       % xm refletido ao lado de alta tensão [Ohm]

%%% ////////////////////////////////////////////////////
%%% PARAMETROS NA BT
%%% ////////////////////////////////////////////////

r2 = rcc / (2 * a^2)                   % Resistência do enrolamento secundário [Ohm]
x2 = xcc / (2 * a^2)                   % Reatância de dispersão de fluxo no secundário [Ohm]

# Refletir Parâmetros

Zeq_nucleo = 1/(1/rp_at + 1/(j*xm_at));    % Impedância equivalente do núcleo do transformador [Ohm]
VNS_r = VNS * a;                           % Tensão no secundário refletida ao primário [V]
r2_r = r2 * a^2;                           % r2 refletido ao primário [Ohm]
x2_r = x2 * a^2;                           % x2 refletido ao primário [Ohm]

%%% ////////////////////////////////////////////////////
%%% CALCULO DE potencias
%%% ////////////////////////////////////////////////

# Funções para cálculo com matrizes
mat_abs = @(x) abs(x);                    % Calcula o módulo de cada elemento de uma dada matriz
mat_conj = @(x) conj(x);                  % Retorna o conjugado de cada elemento de uma dada matriz
mat_real = @(x) real(x);                  % Retorna a parte real de cada elemento de uma dada matriz

# Intervalos de interesse
x = 0.3:0.01:1.4;             % Intervalo [0.3; 1.4] para o gráfico da regulação
x_eta = x(1:101);             % Intervalo [0.3; 1.3] para o gráfico do rendimento

# Cálculo dos valores para os gráficos de rendimento e regulação para fp = 1
Scarga_fp1 = x * SN;                % Potência da carga, variando de 0.3 a 1.4 da potência nominal
Icarga_fp1 = Scarga_fp1 / VNS_r;    % Corrente na carga, refletida ao lado de alta tensão, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga

VT_fp1 = VNS_r .+ (r2_r + j * x2_r) * Icarga_fp1;        % Tensão no lado de baixa tensão, refletida ao lado de alta tensão, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga
I0_fp1 = VT_fp1 / Zeq_nucleo;                            % Corrente que passa pelo ramo em paralelo, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga
I1_fp1 = Icarga_fp1 .+ I0_fp1;                           % Corrente que sai do gerador, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga
V1_fp1 = VT_fp1 + (r1 + j * x1) * I1_fp1;                % Tensão no gerador, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga

VTab_fp1 = (Zeq_nucleo / (Zeq_nucleo + r1 + j * x1)) .* V1_fp1;          % Tensão em aberto no lado de baixa tensão, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga

Pcu_fp1 = mat_real((VT_fp1 - VNS_r) .* mat_conj(Icarga_fp1)) + mat_real((V1_fp1 .- VT_fp1) .* mat_conj(I1_fp1));      % Perdas Joule no transformador, associadas ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga
Pfe_fp1 = mat_real(VT_fp1 .* mat_conj(I0_fp1));                          % Perdas por histerese e correntes Foucault no núcleo, associadas ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga

eta_fp1 = Scarga_fp1 ./ (Scarga_fp1 .+ Pcu_fp1 .+ Pfe_fp1);              % Rendimento do transformador, associado ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga
eta_fp1 = eta_fp1(1:101);                                                % Rendimento do transformador, associado ao intervalo de 0.3 a 1.3 da potência nominal sobre a carga
reg_fp1 = (mat_abs(VTab_fp1) - VNS_r) ./ VNS_r * 100;                    % Regulação do transformador, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga

# Cálculo dos valores para os gráficos de rendimento e regulação para fp = 0.7 capacitivo
Scarga_fp = x * (SN * fp - j * SN * sin(acos(fp)));           % Potência da carga, variando de 0.3 a 1.4 da potência nominal
Icarga_fp = mat_conj(Scarga_fp / VNS_r);                      % Corrente na carga, refletida ao lado de alta tensão, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga

VT_fp = VNS_r .+ (r2_r + j * x2_r) * Icarga_fp;               % Tensão no lado de baixa tensão, refletida ao lado de alta tensão, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga
I0_fp = VT_fp / Zeq_nucleo;                                   % Corrente que passa pelo ramo em paralelo, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga
I1_fp = Icarga_fp .+ I0_fp;                                   % Corrente que sai do gerador, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga
V1_fp = VT_fp + (r1 + j * x1) * I1_fp;                        % Tensão no gerador, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga

VTab_fp = (Zeq_nucleo / (Zeq_nucleo + r1 + j * x1)) .* V1_fp; % Tensão em aberto no lado de baixa tensão, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga

Pcu_fp = mat_real((VT_fp - VNS_r) .* mat_conj(Icarga_fp)) + mat_real((V1_fp .- VT_fp) .* mat_conj(I1_fp));           % Perdas Joule no transformador, associadas ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga
Pfe_fp = mat_real(VT_fp .* mat_conj(I0_fp));                               % Perdas por histerese e correntes Foucault no núcleo, associadas ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga

eta_fp = mat_real(Scarga_fp) ./ (mat_real(Scarga_fp) .+ Pcu_fp .+ Pfe_fp); % Rendimento do transformador, associado ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga
eta_fp = eta_fp(1:101);                                                    % Rendimento do transformador, associado ao intervalo de 0.3 a 1.3 da potência nominal sobre a carga
reg_fp = (mat_abs(VTab_fp) - VNS_r) ./ VNS_r * 100;                        % Regulação do transformador, associada ao intervalo de 0.3 a 1.4 da potência nominal sobre a carga

# Gráfico do rendimento do transformador em função da potência sobre a carga
fig=figure(1, "position", [100,-50,1920,1080]);
titulo=['EC1 - Grupo z',num2str(grupo),' - ', num2str(nusp),' - ', date()];
subplot(2,1,1) % first subplot RENDIMENTO
plot(x_eta, eta_fp1, 'r-', x_eta, eta_fp, 'b-', x_eta(71), eta_fp1(71), 'mo', 'markersize', 10, 
x_eta(71), eta_fp(71), 'mo', 'markersize', 10, [1,1],[0.982,eta_fp1(71)], 'k--');
title(titulo);
text (x(71), eta_fp1(71) + 0.0005, "Condições Nominais", 'fontsize', 15, "horizontalalignment", "center");
text (x(71) + 0.005, eta_fp1(71) - 0.0005, ["(", num2str(eta_fp1(71)), ")"], 'fontsize', 15);
text (x(71) + 0.005, eta_fp(71) - 0.0005, ["(", num2str(eta_fp(71)), ")"], 'fontsize', 15);
set(gca, 'FontSize', 10)
xlabel('% Snom');
ylabel('Rendimento (%)');
ylim([0.982,0.995]);
grid on;
legend("Rendimento para fp = 1", 'Rendimento para fp = 0.7 cap', "box", "on", "location", "southeastoutside", "fontsize", 10);

# Gráfico da regulação do transformador em função da potência sobre a carga
subplot(2,1,2) % second subplot REGULACAO
plot(x, reg_fp1, 'r-', x, reg_fp, 'b--', x(71), reg_fp1(71), 'mo', 'markersize', 10, 
x(71), reg_fp(71), 'mo', 'markersize', 10, [1,1],[-4.5,reg_fp1(71)], 'k--');
text (x(71), reg_fp1(71) + 0.3, "Condições Nominais", 'fontsize', 15, "horizontalalignment", "center");
text (x(71) + 0.005, reg_fp1(71) - 0.3, ["(", num2str(reg_fp1(71)), ")"], 'fontsize', 15);
text (x(71) + 0.005, reg_fp(71) + 0.3, ["(", num2str(reg_fp(71)), ")"], 'fontsize', 15);
set(gca, 'FontSize', 10)
xlabel('% Snom');
ylabel('Regulação (%)');
ylim([-4.5,2]);
grid on;
legend("Regulação para fp = 1", 'Regulação para fp = 0.7 cap', "box", "on", "location", "southeastoutside", "fontsize", 10);

arq=['EC1_2021_PEA3306_z',num2str(grupo),'_',num2str(nusp),'.png']; % GERACAO DO ARQ IMAGEM
print(fig, arq); % Print dos gráficos
