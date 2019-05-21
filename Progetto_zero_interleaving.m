%PROGETTO ZERO INTERLEAVING

%
% Lo script deve ricostruire una sequenza x_n a partire da una sua versione "zero interleaved" indicata con y_n
% in cui ogni M campioni sono inseriti M-1 zeri (al posto dei campioni originari)

% Esempio:
% sequenza originaria                       x_n = 1,2,3,4,5,6,7,8,9,10...
% sequenza zero_interleaved (fattore M = 3)  y_n= 1,0,0,4,0,0,7,0,0,10...

% La sequenza di partenza x_n ha 400 campioni ed è fornita nel file zerointerleaving.mat

% Si richiede che lo script Matlab esegua le seguenti operazioni:

% - generi la sequeza zero-interleaved y_n con M variabile (M = 2,3,...),
% cioé si deve poter scegliere il numero di campioni (M-1) da azzerare.

% - permetta di scegliere la posizione dei campioni non nulli, es. si chiede la
% possibilità di generare M possibili sequenze con campioni non nulli in posizione
% Mk, Mk+1, Mk+2, ... Mk+(M-1) con k = 0,1,2,...

% - rappresenti graficamente la sequenza x_n e le sequenze y_n (zero_interleaved)
% nel dominio del tempo e delle frequenze (trasformata discreta di Fourier).
%
% - rappresenti graficamente nel dominio del tempo e delle frequenze il filtro impiegato per la ricostruzione
% della sequenza x_n

% - rappresenti graficamente la sequenza originaria x_n e la sua versione
% ricostruita, mostrando che la sequenza ricostruita è la stessa qualunque
% sia la scelta della posizione dei campioni non nulli (purché non vi siano fenomeni di alias).

% - permetta di valutare il massimo valore di M che non produce distorsione del segnale ricostrutito.

% Strutturare lo script Matlab in modo tale che possa essere rapidamente adattato per gestire sequenze in ingresso
% diverse da quella assegnata (es. diverso numero di campioni).
% Definire in modo chiaro le variabili utilizzate e commentare sinteticamente i vari passi dello script.

%START
close all;
%clear;
clc;
prompt = ('Inserisci nome file :  ');
nome_file = input (prompt, 's');
load (nome_file);
%si può trovare un modo per caricare direttamente in una variabile?
%load zerointerleaving.mat
y = x; %crea un secondo vettore per non modificare sequenza originale
dim = length(x);
M = input('Inserisci M: ');
f_s = 1/M; %frequenza di campionamento, ma potremmo toglierla perché non serve
%f  = -f_s : 1 : f_s;

%ZERO-INTERLEAVING
y_n = zeros(M,dim); %crea una matrice di zeri contenente le M sequenze lungo le righe
Yf_n = y_n; %copia la matrice creata in quella che sarà la trasformata
for j = 1:M %ciclo per creare sequenze campionate
    i = j;
    while i<dim+1
        y_n(j,i) = y(i);
        i = i+M;
    end
end

%TRASFORMAZIONI 
n = (0:dim-1);
Xf=fft(y);
figure
subplot(2,1,1)
stem (n,y);
xlabel ('Campioni')
title('Sequenza di Partenza')
pause
subplot(2,1,2)
stem (n,real(Xf));
xlabel ('Frequenza (Hz)')
title('Trasf seq X')
pause
for k = 1:M
    titolo = 'Sequenza con primo elemento diverso da 0 in posizione %d';
    pos = k-1;
    Yf_n(k,:) = fft(y_n(k,:)); 
    figure(2)
    subplot(2,1,1)
    stem(n,y_n(k,:));
    xlabel ('Campioni')
    title(sprintf(titolo,pos))
    subplot(2,1,2)
    stem(n,real(Yf_n(k,:)));
    xlabel ('Frequenza (Hz)')
    title('Sequenza trasformata')
    pause
end

%FILTRO
figure
t=-dim/2:dim/2-1;
filtro_t=sinc(t/M);
filtro=abs(fft(filtro_t));


subplot (2,1,1)
stem(filtro);
title('Filtro in frequenza')
subplot (2,1,2)
plot(t,filtro_t);
grid on;
title('Filtro nei tempi')
pause
Zf_n=zeros(M,dim);
z=zeros(M,(2*dim)-1);
u=zeros(M,dim);

for i=1:M
    figure(3)
    Zf_n(i,:)=Yf_n(i,:).*filtro;
    z(i,:) = conv(y_n(i,:),filtro_t);
    u(i,:)=z(i,(dim/2)+1:(3/2)*dim);
    
    subplot(2,1,1)
    stem(u(i,:));
    subplot(2,1,2)
    stem(real(Zf_n(i,:)));
%     titolo='Errore quadratico medio sequenza %d = %0.6f';
%     disp(sprintf(titolo,i,errorequadratico));
        
%         subplot(4,1,4)
%         stem(errorequadraticofrequenza(i,:));
        %pause
    %la convoluzione ancora è da controllare,
    %ma almeno la forma c'è; 
    %il picco è un po' più in alto
    %Non torna il numero di campioni del plot, 
    %dovrebbe essere 400 anche dopo la convoluzione
    
    %eliminare i campioni in eccesso
end


r_n=zeros(dim);
%calcolo la matrice contenente tutte le sequenze con primo campione non
%nullo, al variare di M (che va da 1 a dim)
%calcolo la matrice contentente le sequenze ricostruite per ogni riga della
%matrice r_n

for i=1:dim

l=1;
    while l<dim+1
    r_n(i,l)=y(l);
    l=l+i;
    
    end
   rec_sign(i,:) = conv(r_n(i,:),filtro_t);
     sign(i,:)=rec_sign(i,(dim/2)+1:(3/2)*dim);     
         
    

end
%calcolo errore quadratico medio per tutti gli M
errorequadratico=zeros(1,dim);
for i=1:dim
somma=0;
for j=1:dim
somma=somma+(abs(y(j)-sign(i,j))^2);
end
errorequadratico(i)=(errorequadratico(i)+somma)/dim;

end


% 
