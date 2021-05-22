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
%%%% SN   Pot�ncia nominal                    2.5   [MVA]
%%%% VNP  Tens�o Nominal no Prim�rio          19.1  [kV]
%%%% VNS  Tens�o Nominal no Secund�rio        3.81  [kV]
%%%% f    Frequ�ncia de Opera��o              50    [Hz]
%%%% Vv   Tens�o no Ensaio em Vazio           3810  [V]
%%%% Iv   Corrente no Ensaio em Vazio         9.86  [A]
%%%% Pv   Pot�ncia Ativa no Ensaio em Vazio   8140  [W]
%%%% Vcc  Tens�o no Ensaio em Curto           854   [V]
%%%% Icc  Corrente no Ensaio em Curto         131   [A]
%%%% Pcc  Pot�ncia no Ensaio em Curto         8890  [W]
%%%% a    Rela��o de Transforma��o
%%%% fp   Fator de Pot�ncia da Carga          0.7 cap
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

%%% Par�metros calculados a partir do ensaio em vazio (baixa tens�o)

fpv = Pv / (Vv * Iv);                  % Fator de pot�ncia no ensaio em vazio (fp = P/|S| = P/(|V|*|I|))
Ipv = Iv * fpv;                        % Corrente de perdas no ferro no ensaio em vazio [A]
Imv = Iv * sqrt(1 - fpv^2);            % Corrente de magnetiza��o no ensaio em vazio [A]

rp = Vv/Ipv                            % Resist�ncia de perdas por histerese e correntes de Foucault no n�cleo [Ohm]
xm = Vv/Imv                            % Reat�ncia de magnetiza��o no n�cleo. [Ohm]

%%% Par�metros calculados a partir do ensaio em curto (alta tens�o)

Zcc = Vcc/Icc;                         % Imped�ncia de curto circuito [Ohm]
rcc = Pcc/(Icc^2);                     % Resist�ncia de curto circuito [Ohm]
xcc = sqrt(Zcc^2 - rcc^2);             % Reat�ncia de curto circuito [Ohm]

r1 = rcc / 2                           % Resist�ncia do enrolamento prim�rio [Ohm]
x1 = xcc / 2                           % Reat�ncia de dispers�o de fluxo no prim�rio [Ohm]

%%% ////////////////////////////////////////////////////
%%% PARAMETROS NA AT
%%% ////////////////////////////////////////////////

rp_at = rp * a^2                       % rp refletido ao lado de alta tens�o [Ohm]
xm_at = xm * a^2                       % xm refletido ao lado de alta tens�o [Ohm]

%%% ////////////////////////////////////////////////////
%%% PARAMETROS NA BT
%%% ////////////////////////////////////////////////

r2 = rcc / (2 * a^2)                   % Resist�ncia do enrolamento secund�rio [Ohm]
x2 = xcc / (2 * a^2)                   % Reat�ncia de dispers�o de fluxo no secund�rio [Ohm]

# Refletir Par�metros

Zeq_nucleo = 1/(1/rp_at + 1/(j*xm_at));    % Imped�ncia equivalente do n�cleo do transformador [Ohm]
VNS_r = VNS * a;                           % Tens�o no secund�rio refletida ao prim�rio [V]
r2_r = r2 * a^2;                           % r2 refletido ao prim�rio [Ohm]
x2_r = x2 * a^2;                           % x2 refletido ao prim�rio [Ohm]

%%% ////////////////////////////////////////////////////
%%% CALCULO DE potencias
%%% ////////////////////////////////////////////////

# Fun��es para c�lculo com matrizes
mat_abs = @(x) abs(x);                    % Calcula o m�dulo de cada elemento de uma dada matriz
mat_conj = @(x) conj(x);                  % Retorna o conjugado de cada elemento de uma dada matriz
mat_real = @(x) real(x);                  % Retorna a parte real de cada elemento de uma dada matriz

# Intervalos de interesse
x = 0.3:0.01:1.4;             % Intervalo [0.3; 1.4] para o gr�fico da regula��o
x_eta = x(1:101);             % Intervalo [0.3; 1.3] para o gr�fico do rendimento

# C�lculo dos valores para os gr�ficos de rendimento e regula��o para fp = 1
Scarga_fp1 = x * SN;                % Pot�ncia da carga, variando de 0.3 a 1.4 da pot�ncia nominal
Icarga_fp1 = Scarga_fp1 / VNS_r;    % Corrente na carga, refletida ao lado de alta tens�o, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga

VT_fp1 = VNS_r .+ (r2_r + j * x2_r) * Icarga_fp1;        % Tens�o no lado de baixa tens�o, refletida ao lado de alta tens�o, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga
I0_fp1 = VT_fp1 / Zeq_nucleo;                            % Corrente que passa pelo ramo em paralelo, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga
I1_fp1 = Icarga_fp1 .+ I0_fp1;                           % Corrente que sai do gerador, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga
V1_fp1 = VT_fp1 + (r1 + j * x1) * I1_fp1;                % Tens�o no gerador, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga

VTab_fp1 = (Zeq_nucleo / (Zeq_nucleo + r1 + j * x1)) .* V1_fp1;          % Tens�o em aberto no lado de baixa tens�o, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga

Pcu_fp1 = mat_real((VT_fp1 - VNS_r) .* mat_conj(Icarga_fp1)) + mat_real((V1_fp1 .- VT_fp1) .* mat_conj(I1_fp1));      % Perdas Joule no transformador, associadas ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga
Pfe_fp1 = mat_real(VT_fp1 .* mat_conj(I0_fp1));                          % Perdas por histerese e correntes Foucault no n�cleo, associadas ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga

eta_fp1 = Scarga_fp1 ./ (Scarga_fp1 .+ Pcu_fp1 .+ Pfe_fp1);              % Rendimento do transformador, associado ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga
eta_fp1 = eta_fp1(1:101);                                                % Rendimento do transformador, associado ao intervalo de 0.3 a 1.3 da pot�ncia nominal sobre a carga
reg_fp1 = (mat_abs(VTab_fp1) - VNS_r) ./ VNS_r * 100;                    % Regula��o do transformador, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga

# C�lculo dos valores para os gr�ficos de rendimento e regula��o para fp = 0.7 capacitivo
Scarga_fp = x * (SN * fp - j * SN * sin(acos(fp)));           % Pot�ncia da carga, variando de 0.3 a 1.4 da pot�ncia nominal
Icarga_fp = mat_conj(Scarga_fp / VNS_r);                      % Corrente na carga, refletida ao lado de alta tens�o, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga

VT_fp = VNS_r .+ (r2_r + j * x2_r) * Icarga_fp;               % Tens�o no lado de baixa tens�o, refletida ao lado de alta tens�o, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga
I0_fp = VT_fp / Zeq_nucleo;                                   % Corrente que passa pelo ramo em paralelo, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga
I1_fp = Icarga_fp .+ I0_fp;                                   % Corrente que sai do gerador, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga
V1_fp = VT_fp + (r1 + j * x1) * I1_fp;                        % Tens�o no gerador, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga

VTab_fp = (Zeq_nucleo / (Zeq_nucleo + r1 + j * x1)) .* V1_fp; % Tens�o em aberto no lado de baixa tens�o, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga

Pcu_fp = mat_real((VT_fp - VNS_r) .* mat_conj(Icarga_fp)) + mat_real((V1_fp .- VT_fp) .* mat_conj(I1_fp));           % Perdas Joule no transformador, associadas ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga
Pfe_fp = mat_real(VT_fp .* mat_conj(I0_fp));                               % Perdas por histerese e correntes Foucault no n�cleo, associadas ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga

eta_fp = mat_real(Scarga_fp) ./ (mat_real(Scarga_fp) .+ Pcu_fp .+ Pfe_fp); % Rendimento do transformador, associado ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga
eta_fp = eta_fp(1:101);                                                    % Rendimento do transformador, associado ao intervalo de 0.3 a 1.3 da pot�ncia nominal sobre a carga
reg_fp = (mat_abs(VTab_fp) - VNS_r) ./ VNS_r * 100;                        % Regula��o do transformador, associada ao intervalo de 0.3 a 1.4 da pot�ncia nominal sobre a carga

# Gr�fico do rendimento do transformador em fun��o da pot�ncia sobre a carga
fig=figure(1, "position", [100,-50,1920,1080]);
titulo=['EC1 - Grupo z',num2str(grupo),' - ', num2str(nusp),' - ', date()];
subplot(2,1,1) % first subplot RENDIMENTO
plot(x_eta, eta_fp1, 'r-', x_eta, eta_fp, 'b-', x_eta(71), eta_fp1(71), 'mo', 'markersize', 10, 
x_eta(71), eta_fp(71), 'mo', 'markersize', 10, [1,1],[0.982,eta_fp1(71)], 'k--');
title(titulo);
text (x(71), eta_fp1(71) + 0.0005, "Condi��es Nominais", 'fontsize', 15, "horizontalalignment", "center");
text (x(71) + 0.005, eta_fp1(71) - 0.0005, ["(", num2str(eta_fp1(71)), ")"], 'fontsize', 15);
text (x(71) + 0.005, eta_fp(71) - 0.0005, ["(", num2str(eta_fp(71)), ")"], 'fontsize', 15);
set(gca, 'FontSize', 10)
xlabel('% Snom');
ylabel('Rendimento (%)');
ylim([0.982,0.995]);
grid on;
legend("Rendimento para fp = 1", 'Rendimento para fp = 0.7 cap', "box", "on", "location", "southeastoutside", "fontsize", 10);

# Gr�fico da regula��o do transformador em fun��o da pot�ncia sobre a carga
subplot(2,1,2) % second subplot REGULACAO
plot(x, reg_fp1, 'r-', x, reg_fp, 'b--', x(71), reg_fp1(71), 'mo', 'markersize', 10, 
x(71), reg_fp(71), 'mo', 'markersize', 10, [1,1],[-4.5,reg_fp1(71)], 'k--');
text (x(71), reg_fp1(71) + 0.3, "Condi��es Nominais", 'fontsize', 15, "horizontalalignment", "center");
text (x(71) + 0.005, reg_fp1(71) - 0.3, ["(", num2str(reg_fp1(71)), ")"], 'fontsize', 15);
text (x(71) + 0.005, reg_fp(71) + 0.3, ["(", num2str(reg_fp(71)), ")"], 'fontsize', 15);
set(gca, 'FontSize', 10)
xlabel('% Snom');
ylabel('Regula��o (%)');
ylim([-4.5,2]);
grid on;
legend("Regula��o para fp = 1", 'Regula��o para fp = 0.7 cap', "box", "on", "location", "southeastoutside", "fontsize", 10);

arq=['EC1_2021_PEA3306_z',num2str(grupo),'_',num2str(nusp),'.png']; % GERACAO DO ARQ IMAGEM
print(fig, arq); % Print dos gr�ficos
