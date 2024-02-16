 clc

FoldIn = '/Users/camillabandinelli/Desktop'; 

% Specifica il nome dell'immagine
RootIn = 'angolo_IMU_pitot.jpeg'; 

% Componi il percorso completo dell'immagine
sIn = fullfile(FoldIn, RootIn);

% Carica l'immagine utilizzando imread
I = imread(sIn); 

% Mostra l'immagine a video in scala di grigi
imshow(rgb2gray(I));


%dierction vector of x-axes
initial_point1=[490,506];
final_point1=[202,506];



%direction vector pitot tube