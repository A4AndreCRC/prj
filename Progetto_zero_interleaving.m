%PROGETTO ZERO INTERLEAVING

%
% Lo script deve ricostruire una sequenza x_n a partire da una sua versione "zero interleaved" indicata con y_n
% in cui ogni M campioni sono inseriti M-1 zeri (al posto dei campioni originari)

% Esempio:
% sequenza originaria                       x_n = 1,2,3,4,5,6,7,8,9,10...
% sequenza zero_interleaved (fattore M = 3)  y_n= 1,0,0,4,0,0,7,0,0,10...

% La sequenza di partenza x_n ha 400 campioni ed Ã¨ fornita nel file zerointerleaving.mat

% Si richiede che lo script Matlab esegua le seguenti operazioni:

% - generi la sequeza zero-interleaved y_n con M variabile (M = 2,3,...),
% cioÃ© si deve poter scegliere il numero di campioni (M-1) da azzerare.

% - permetta di scegliere la posizione dei campioni non nulli, es. si chiede la
% possibilitÃ  di generare M possibili sequenze con campioni non nulli in posizione
% Mk, Mk+1, Mk+2, ... Mk+(M-1) con k = 0,1,2,...

% - rappresenti graficamente la sequenza x_n e le sequenze y_n (zero_interleaved)
% nel dominio del tempo e delle frequenze (trasformata discreta di Fourier).
%
% - rappresenti graficamente nel dominio del tempo e delle frequenze il filtro impiegato per la ricostruzione
% della sequenza x_n

% - rappresenti graficamente la sequenza originaria x_n e la sua versione
% ricostruita, mostrando che la sequenza ricostruita Ã¨ la stessa qualunque
% sia la scelta della posizione dei campioni non nulli (purchÃ© non vi siano fenomeni di alias).

% - permetta di valutare il massimo valore di M che non produce distorsione del segnale ricostrutito.

% Strutturare lo script Matlab in modo tale che possa essere rapidamente adattato per gestire sequenze in ingresso
% diverse da quella assegnata (es. diverso numero di campioni).
% Definire in modo chiaro le variabili utilizzate e commentare sinteticamente i vari passi dello script.

close all
clear all
clc

prompt = ('Inserisci nome file :  ');
nome_file = input (prompt, 's');
load (nome_file);

y = x; %crea un secondo vettore per non modificare sequenza originale
dim = length(x);
M = input('Inserisci M: ');
y_n = zeros(M,dim); %crea una matrice di zeri contenente le M sequenze lungo le righe
Yf_n = y_n; %copia la matrice di zeri creata in quella che sarà  la trasformata
n = (0:dim-1); %intervallo di rappresentazione
Xf=fft(y);%trasformata della sequenza originaria
Zf_n=zeros(M,dim);%matrice delle trasformate delle singole sequenze
z=zeros(M,(2*dim)-1);
u=zeros(M,dim);

%definizione asse tempi per filtro
if mod(dim,2)==0
    
    t=-dim/2:dim/2-1;
    
else
    
    t=-floor(dim/2):floor(dim/2);
    
end

%creazione matrice zerointerleaved

for j = 1:M %scorre le M righe corrispondenti alle M sequenze zerointerleaved
    
    i = j;
    while i<dim+1 %ciclo per creare sequenze campionate

        y_n(j,i) = y(i);
        i = i+M;
        
    end
    
end

%rappresentazione delle sequenze zerointerleaved nel dominio dei tempi e
%frequenze

figure(1)
subplot(2,1,1)
stem (n,y);
xlabel ('Campioni')
title('Sequenza di partenza')
subplot(2,1,2)
stem (n,real(Xf));
xlabel ('Campioni')
title('Trasformata della sequenza di partenza')
pause
for k = 1:M
    titolo = 'Sequenza con primo elemento diverso da 0 in posizione %d';
    pos = k;
    Yf_n(k,:) = fft(y_n(k,:));
    figure (2)
    subplot(2,1,1)
    stem(n,y_n(k,:));
    xlabel ('Campioni')
    title(sprintf(titolo,pos))
    subplot(2,1,2)
    stem(n,real(Yf_n(k,:)));
    xlabel ('Campioni')
    title('Sequenza trasformata')
    pause
end


%creazione del filtro di ricostruzione
filtro_t=sinc(t/M);%filtro nei tempi
filtro=abs(fft(filtro_t));%filtro in frequenza

%rappresentazione filtro 
figure (3)
subplot (2,1,2)
stem(filtro);
title('Filtro in frequenza')
subplot (2,1,1)
plot(t,filtro_t);
grid on;
title('Filtro nei tempi')
pause


%ricostruzione della sequenza originaria a partire dalle sequenze
%zerointerleaved
for i=1:M
    
    figure (4)
    Zf_n(i,:)=Yf_n(i,:).*filtro;%ricostruzione in frequenza
    z(i,:) = conv(y_n(i,:),filtro_t);%ricostruzione nei tempi
    u(i,:)=z(i,floor(dim/2)+1:fix((3/2)*dim));%eliminazione campioni aggiunti
    % dalla convoluzione

    %rappresentazione della sequenza ricostruita al variare del primo elemento 
    %non nullo con confronto con la sequenza originaria
    titolo1='Confronto sequenze nei tempi con primo elemento diverso da 0 in posizione %d';
    titolo2='Confronto sequenze nelle frequenze con primo elemento diverso da 0 in posizione %d';
    subplot(2,1,1)
    stem(n,u(i,:),'filled','DisplayName','Sequenza ricostruita');
    xlabel ('Campioni')
    legend;
    title(sprintf(titolo1,i))
    hold on
    stem (n,x,'r','DisplayName','Sequenza di partenza');
    xlabel ('Campioni')
    legend;
    hold off
    
    subplot(2,1,2)
    stem(n,real(Zf_n(i,:)),'filled','DisplayName','Sequenza ricostruita');
    xlabel('Campioni')
    legend;
   title(sprintf(titolo2,i))
    hold on
    stem (n,real(Xf),'r','DisplayName','Sequenza di partenza');
    xlabel('Campioni')
    legend;
    hold off
    pause
    
    
end




%inizializzazione variabili per il calcolo dell'errore
r_n=zeros(dim);
filter_t=zeros(dim);
rec_sign=zeros(M,(2*dim)-1);
sign=zeros(dim);
errorequadratico=zeros(1,dim);


for i=1:dim
    
    filter_t(i,:)=sinc(t/i);
    
end

for i=1:dim
%calcolo la matrice contenente tutte le sequenze con primo campione non
%nullo, al variare di M (che va da 1 a dim)
 
    l=1;
    while l<dim+1
        r_n(i,l)=y(l);
        l=l+i;
        
    end
    %calcolo la matrice rec_sign contentente le sequenze ricostruite per
    % ogni riga della matrice r_n
    rec_sign(i,:) = conv(r_n(i,:),filter_t(i,:));
    %eliminazione campioni aggiunti dalla convoluzione
    sign(i,:)=rec_sign(i,floor(dim/2)+1:fix((3/2)*dim));
end

%calcolo errore quadratico medio per tutti gli M
for i=1:dim
    
    errorequadratico(i,:)=immse(y,sign(i,:));
    
end

flag=0;%variabile flag per ripetere il menù

while flag==0%apertura menù di scelta
    
    err=input('inserisci errore quadratico medio massimo accettabile: ');
    i=1;
    while(errorequadratico(i)<err && i<dim)
        i=i+1;
        
    end
    %gestione caso in cui si inserisce un valore di errore troppo alto
    if i==dim
        fprintf("errore inserito troppo alto!\n");
        scelta=input('ripetere operazione per errore diverso? s--> si, altri tasti--> no ','s');
        if scelta~='s'
            flag=1;
        end
        
    else
        
        erroreq=zeros(1,i);%crea vettore di zeri di lunghezza pari alla lunghezza i
       for j=1:i
            
         erroreq(j) = errorequadratico(j);%copia i primi i valori del 
         %vettore errorequadratico in erroreq
            
            end
        threshold (1:j+1) = err;
        
        %rappresentazione grafica dell'errore
        figure (5)
        stem(1:i,erroreq,'DisplayName','Errore Quadratico Medio')
        title('Correlazione tra M e MSE')
        xlabel ('Valori di M')
        ylabel ('MSE')
        xticks(1:1:i+1)
        xlim([1 i+1])
        grid on
        hold on
        plot (threshold,'r','LineStyle','--','DisplayName','Threshold')
        legend('Location','North');
        hold off
        output = '\nIl massimo valore di M tale per cui la distorsione è inferiore al valore richiesto risulta essere %d ';
        fprintf(output,i-1)
        scelta=input('\nripetere operazione per errore diverso? s--> si, altri tasti--> no ','s');
        if scelta~='s'
            flag=1;
        end
        
    end
end
