% Parametri del segnale
n = 1000; % Lunghezza del segnale
t = (0:n-1)'; % Tempo

% Genera tre segnali sinusoidali con frequenze casuali
freq1 = rand()*0.1; % Frequenza 1
freq2 = rand()*0.1; % Frequenza 2
freq3 = rand()*0.1; % Frequenza 3
sinusoid1 = sin(2*pi*freq1*t);
sinusoid2 = sin(2*pi*freq2*t);
sinusoid3 = sin(2*pi*freq3*t);

% Combina i segnali sinusoidali
signal = sinusoid1 + sinusoid2 + sinusoid3;

% Aggiungi del rumore gaussiano
mean_val = 0; % Valore medio del rumore
std_dev = 0.1; % Deviazione standard del rumore
noise = mean_val + std_dev * randn(n, 1);
signal = signal + noise;

% Calcola l'ACF
[acf,lag] = autocorr(signal,100); % Calcola l'ACF normalizzata

% Calcola la PDF
nbins = 50; % Numero di bin per l'istogramma
pdf_edges = linspace(min(signal), max(signal), nbins+1); % Bins edges
pdf_centers = (pdf_edges(1:end-1) + pdf_edges(2:end)) / 2; % Centri dei bin
pdf = histcounts(signal, pdf_edges, 'Normalization', 'pdf'); % Calcola l'istogramma normalizzato

% Calcola la CDF
cdf = cumsum(pdf); % Calcola la CDF come la somma cumulativa della PDF

% Plot
figure;
plot(signal, LineWidth=2);
xlabel('Time [s]','FontSize',16);
ylabel('Velocity [m/s]','FontSize',16);

figure
plot(lag,acf, LineWidth=2);
xlabel('Lag','FontSize',16);
ylabel('ACF','FontSize',16);

figure
bar(pdf_centers, pdf, 'hist');
xlabel('Velocity [m/s]','FontSize',16);
ylabel('PDF','FontSize',16);

figure;
plot(pdf_centers, cdf, LineWidth=2);
xlabel('Velocity [m/s]','FontSize',16);
ylabel('CDF','FontSize',16);