function [ vector_caracteristicas, vector_intervalos, columnas ] = obtieneCaracteristicas2( x1 )
% fs = 200;              % tasa de muestreo
% fs = 360;              % tasa de muestreo de MIT-BIH
fs = 250;              % tasa de muestreo de Euro ST-T

N = length (x1);       % longitud de la señal
t = [0:N-1]/fs;        % indice temporal

%  figure(1)
% % subplot(2,1,1)
% %  xlabel('Segundos');ylabel('mVolts');title('señal ECG de entrada')
% % 
% %subplot(2,1,2)
%  plot(t,x1)
%  grid on
% xlabel('segundos');ylabel('mVolts');title('segundos 1-3 en la señal de entrada')
%  xlim([0 3])
%xlim([1 5])
 
%NORMALIZACION Y CANCELACION DEL COMPONENTE DE CD
x1 = x1 - mean (x1);        % cancelando componente CD
x1 = x1 / max( abs(x1));    % normalizando

% eliminacion  de corrimiento de linea
% n=72;
% Ventana de 200ms (?)
n = (200 * fs)/1000;
ki = 1;
for k=1:length(x1)
    linea0 = zeros(1,n);
    lim1 = k - (n/2);
    lim2 = k + (n/2) - 1;
    if lim2 > length(x1)
        lim2 = length(x1);
    end
    if lim1 <= 0
        for k2 = (n/2) + lim1: lim2
            linea0(k2) = x1(k);
            k=k+1;
        end
    else
        kl = 1;
        for k2 = lim1: lim2
            linea0(kl) = x1(k2);
            kl=kl+1;
        end 
    end
    linea0 = sort(linea0);
    if lim1 <= 0 
        prom = (linea0(n+lim1-1)+linea0(n+lim1))/2;
    else
        prom = (linea0(n/2)+linea0(n/2 +1))/2;
    end
    linea1(ki) = prom;
    ki=ki+1;
end

% ventana 600mseg
% n=216;
% Cálculo de número de muestras para ventana de 600ms
n = (600 * fs)/1000;
ki = 1;
for k=1:length(x1)
    linea0 = zeros(1,n);
    lim1 = k - (n/2);
    lim2 = k + (n/2) - 1;
    if lim2 > length(x1)
        lim2 = length(x1);
    end
    if lim1 <= 0
        for k2 = (n/2) + lim1: lim2
            linea0(k2) = linea1(k);
            k=k+1;
        end
    else
        kl = 1;
        for k2 = lim1: lim2
            linea0(kl) = linea1(k2);
            kl=kl+1;
        end 
    end
    linea0 = sort(linea0);
    if lim1 <= 0 
        prom = (linea0(n+lim1-1)+linea0(n+lim1))/2;
    else
        prom = (linea0(n/2)+linea0(n/2+1))/2;
    end
    linea2(ki) = prom;
    ki=ki+1;
end

x1=x1-linea2;

% lineaBase = medfilt1(x1,72); % filtro mediana de 200 mseg de ventana (72 muestras
% lineaBase = medfilt1(lineaBase,216);% filtro de mediana de 600 mseg de ventana (216 muestras)
% x1= x1 - lineaBase;

%  figure(2)
% % % %subplot(2,1,1)
%  plot(t,x1)
%  xlabel('Segundos');ylabel('mVolts');title(' Señal ECG despues de cancelacion y normalizacion del componente de CD')
% % % %subplot(2,1,2)
% % % %plot(t(360:1800),x1(360:1800))
%  grid on
% % % %xlabel('segundos');ylabel('mVolts');title('Segundos 1-5 de señal ECG')
%   xlim([0 3])
% %xlim([1 3])

%FILTRO PASO BAJO (LPF)

% LPF (1-z^-6)^2/(1-z^-1)^2
b=[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a=[1 -2 1];

h_LP=filter(b,a,[1 zeros(1,12)]);   % Función de transferencia del LPF

x2 = conv (x1 ,h_LP);
x2 = x2/ max( abs(x2 ));            % Normalizar.

x22=x2;

%  figure(3)
% %subplot(2,1,1)
%  plot([0:length(x2)-1]/fs,x2)
% xlabel('Segundos');ylabel('mVolts');title(' ECG despues de LPF')
%  grid on
% xlim([ 0 3])
% %xlim([0 max(t)])
%subplot(2,1,2)
%plot(t(200:600),x2(200:600))
%plot(t(360:1800),x2(360:1800))
%xlabel('segundos');ylabel('Volts');title(' ECG de 1-5 segundos')
%xlim([1 3])
%xlim([1 5])

%FILTRADO PASO ALTO (HPF)

% HPF = Allpass-(Lowpass) = z^-16-[(1-z^-32)/(1-z^-1)]
b = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a = [1 -1];

h_HP=filter(b,a,[1 zeros(1,32)]);   % Respuesta a señal impulso del HPF

x3 = conv (x2 ,h_HP);
x3 = x3/ max( abs(x3 ));

%local = 'C:\ProyectoCardio\BDArchivosMIT-BIH\';
%ruta = strcat(local,'senalFiltrada.dat');
%ex=fopen(ruta,'w');
%for fila=1:length(x2)
%    fwrite(ex,x2(fila),'float');
%end
%fclose(ex);

%  figure(4)
% % % subplot(2,1,1)
%   plot([0:length(x3)-1]/fs,x3)
%   xlabel('Segundos');ylabel('mVolts');title(' Señal ECG despues de HPF')
%   grid on
%  xlim([ 0 3])
% %xlim([0 max(t)])
% %subplot(2,1,2)
%plot(t(360:1800),x3(360:1800))
%xlabel('segundos');ylabel('Volts');title(' ECG de 1-5 segundos')
%xlim([1 3])
%xlim([1 5])
% local = 'C:\ProyectoCardio\BDArchivosMIT-BIH\';
% ruta = strcat(local,'senal100.dat');
% ex=fopen(ruta,'w');
% for fila=1:length(x2)
%     
%     fwrite(ex,x2(fila),'float');
%   
% end
% fclose(ex);


%FILTRO DERIVATIVO

% Respuesta al Impulso
h = [-1 -2 0 2 1]/8;

% Aplicando Filtro
x4 = conv (x3 ,h);
x4 = x4 (2+[1: N]);
x4 = x4/ max( abs(x4 ));

%  figure(5)
% % %subplot(2,1,1)
%  plot([0:length(x4)-1]/fs,x4)
%  xlabel('Segundos');ylabel('mVolts');title(' Señal ECG despues de la Derivada')
%  grid on
%  xlim([ 0 3])
%subplot(2,1,2)
%plot(t(360:1800),x4(360:1800))
%xlabel('segundos');ylabel('Volts');title(' ECG de 1-5 segundos')
%xlim([1 3])
%xlim([1 5])

%SQUARING

x5 = x4 .^2;
x5 = x5/ max( abs(x5) );

%  figure(6)
% % %subplot(2,1,1)
%  plot([0:length(x5)-1]/fs,x5)
%  xlabel('Segundos');ylabel('mVolts');title(' Señal ECG al Cuadrado')
%  grid on
%  xlim([0 3])
% subplot(2,1,2)
% plot(t(360:1800),x5(360:1800))
% xlabel('segundos');ylabel('Volts');title(' ECG de 1-5 segundos')
% xlim([1 3])
% xlim([1 5])

%INTEGRACION DE VENTANA TEMPORAL

%  Respuesta al Impulso
h = ones (1 ,31)/31;
Delay = 15; % Retraso en Muestras

% Aplicando Filtro
x6 = conv (x5 ,h);
x6 = x6 (15+[1: N]);
x6 = x6/ max( abs(x6 ));

%  figure(7)
% % subplot(2,1,1)
%  plot([0:length(x6)-1]/fs,x6)
%  xlabel('Segundos');ylabel('mVolts');title(' Señal ECG Integrada')
%  grid on
% xlim([0 3])
%subplot(2,1,2)
%plot(t(360:1800),x6(360:1800))
%xlabel('segundos');ylabel('Volts');title(' ECG de 1-5 segundos')
%xlim([1 3])
%xlim([1 5])


% ondoleta = calculaOndoleta(x2,4,0);
% ondoleta = ondoleta(1:length(x2),1);
% iCruce=1;
% puntoCero=0;
% bandera = 0; % no ha cruzado
% for ceroO = 1:length(ondoleta)
%     if bandera == 0  && ondoleta(ceroO) > 0
%             vectorCruces(iCruce) = ceroO;
%             iCruce = iCruce+1;
%             bandera = 1; % ya cruzo
%     else
%         if bandera == 1 && ondoleta(ceroO) > 0
%             continue;
%         else if bandera == 1 && ondoleta(ceroO) < 0
%                 bandera = 0;
%             end
%         end
%     end
% end
   
%    figure(77)
% % subplot(2,1,1)
% plot([0:length(ondoleta)-1]/fs,ondoleta)
% xlabel('Segundos');ylabel('mVolts');title(' ondoleta')
% grid on
% xlim([3 6])

%ENCONTRANDO PUNTOS QRS
max_h = max(x6);
thresh = mean (x6);
poss_reg = (x6>thresh*max_h)';

left = find(diff([0;poss_reg])==1);
right = find(diff([poss_reg;0])==-1);

left = left-(6+20);     % Retraso debido al LP y HP
right = right-(6+12);   % Retraso debido al LP y HP

%elimina los negativos y reajusta los vectores left y right
leftPos = find(left > 0);
left = left(leftPos);
right = right(leftPos);

% iko = 1; 
% desde=1;
% for puntosC=1:length(vectorCruces)
%    for puntosR = desde:length(left)
%          if left(puntosR) < vectorCruces(puntosC) && vectorCruces(puntosC) < right(puntosR)
%             correctoR(iko) = vectorCruces(puntosC);
%             iko = iko+1;
%             desde = puntosR;
%             break;
%          else
%             break;
%         end
%     end
% end

k=1;
for i=1:length(left)
    if right(i) - left(i) > 30
        [R_value(k), R_loc(k)] = max( x2(left(i):right(i)) );
        %[R_value(k), R_loc(k)] = max( x2(left(i):right(i)) );
        % R_value(k) = x2(R_loc(k));
        if (R_value(k) < 0)                 % picos negativos de R
            [R_value(k), R_loc(k)] = min( x2(left(i):right(i)) );
            R_loc(k) = R_loc(k)-1+left(i);  % Agregando offset
            
            [Q_value(k), Q_loc(k)] = max( x2(left(i):R_loc(k)) );
            Q_loc(k) = Q_loc(k)-1+left(i);  % Agregando offset

            [S_value(k), S_loc(k)] = max( x2(R_loc(k):right(i)) );
            S_loc(k) = S_loc(k)-1+R_loc(k); % Agregando offset
            k = k+1;
            continue;
        else                                % picos positivos de R
            % [R_value(k), R_loc(k)] = max( x2(left(i):right(i)) );
            R_loc(k) = R_loc(k)-1+left(i);  % Agregando offset

            [Q_value(k), Q_loc(k)] = min( x2(left(i):R_loc(k)) );
            Q_loc(k) = Q_loc(k)-1+left(i);  % Agregando offset

            [S_value(k), S_loc(k)] = min( x2(R_loc(k):right(i)) );
            S_loc(k) = S_loc(k)-1+R_loc(k); % Agregando offset
        end
        k = k+1;
    end
end

% No hay onda selectiva

Q_loc=Q_loc(find(Q_loc~=0));
R_loc=R_loc(find(R_loc~=0));
S_loc=S_loc(find(S_loc~=0));


% #### si hay mas de 4 latidos entonces procesa, de otra manera se regresa
if (length(R_loc) > 4)

% figure(9)
% %subplot(2,1,1)
% plot ([0:length(x22)-1]/fs,x22  , t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o'); 
% title('Señal ECG con puntos Q,R,S');
% grid on
% xlabel('Segundos');ylabel('mVolts')
% legend('ECG','R','S','Q');
% xlim([0 3])
%subplot(2,1,2)
%plot (t,x1/max(x1) , t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o'); title('Señal ECG con puntos Q,R,S');
%xlabel('Segundos');ylabel('mVolts')
%xlim([0 3])
%xlim([0 5])
%figure(10)
%subplot(2,1,1)
%plot (t,x1/max(x1) , t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o'); title('Señal ECG con puntos Q,R,S');
%xlabel('Segundos');ylabel('mVolts')
%xlim([0 1])


%Puntos Q,R y S
%puntos Q
%disp ('puntos Q')
%disp (t(Q_loc))
%puntos R
%disp ('puntos R')
%disp (t(R_loc))
%puntos S
%disp ('puntos S')
%disp (t(S_loc))

%calculo duracion intervalos R-R
RR = 1:length(R_loc)-1;
for k = 1:length(R_loc)-1
    RR(k) = t(R_loc(k+1)) - t(R_loc(k));
end
%calculo duracion intervalos Q-S
QS = 1:length(Q_loc)-1;
for k = 1:length(Q_loc)-1
    QS(k) = t(S_loc(k)) - t(Q_loc(k));
end

%%%% deteccion anchura complejo QRS

longECG = length(x1);
tsamp = 0.0028; % frecuencia muestreo de 360 hz
% Ventana cálculo QRS
msVentanaQRS = 100;
offsetVentanaQRS = round( (msVentanaQRS * fs )/1000 );
% ajuste de vectores R para la ventana
mayoresVentana = find(R_loc > offsetVentanaQRS);
R_loc = R_loc(mayoresVentana);
R_value = R_value(mayoresVentana);
long_R  = length(R_loc);

%indices
iOndaP = 0;
iOndaT = 0;

for i=1:long_R
    %%%% deteccion anchura complejo QRS
    if (R_loc(i)+offsetVentanaQRS+7) < longECG  % 7 desplazamiento
        % vector=x1(t(R_loc(i))-36:t(R_loc(i))+36); %72*0.0028=0.20s 
        vectorV = x2(R_loc(i)-offsetVentanaQRS:R_loc(i)+offsetVentanaQRS);
        vectorT = t(R_loc(i)-offsetVentanaQRS:R_loc(i)+offsetVentanaQRS);
        % wavelet
        % ondoleta2=cwt(vectorV,2,'gaus1');
        ondoleta = calculaOndoleta(vectorV,2,0);
        ondoleta = ondoleta(1:length(vectorV),1);
        % gaussiana = 1:length(vectorV);
        % a=2;
        % b=0;
        % indT=1;
        % calcula primera derivada de una gaussiana
        % for tf=-length(vectorV)/2:(length(vectorV)/2)-1
            % exponente = -0.5*(((tf-b)/a)*((tf-b)/a));
            % display('exp');display(tf);display(exponente);
            % gaussiana(indT) = ((-tf-b)/a)*exp(exponente);
            % display(exp(exponente));
            % exponente = -1*(((tf-b)/a)*((tf-b)/a));
            % gaussiana(indT) = ((-2*tf-b)/a)*((2/3.141516)^0.25)*exp(exponente);
            % indT=indT+1;
        % end
 
        % convolucion
        % ondoletaPirata=zeros(2*length(vectorV));
        % matlab
        % ondoletaPirata = -1*conv(vectorV,gaussiana,'same');
        % implementada
        % fprintf('The inner product of u and v is %1.2f\n', dot(vectorV,gaussiana));
        % fprintf('The inner product of u and v is %1.2f\n', dot(vectorV,ondoleta));
  
        % my_length = length(vectorV) + length(gaussiana) - 1;
        % ondoletaPirata = zeros(my_length, 1 );
        % for iC = 1:my_length
        %   for jC = 1:length(vectorV) 
        %       if( (iC-jC+1) > 0 && (iC-jC+1) < length(gaussiana) ) 
        %           ondoletaPirata(iC) =  ondoletaPirata(iC)+ vectorV(jC) * gaussiana(iC-jC+1); 
        %       end
        %   end
        % end
  
        % jC=1; 

        % for iC = ceil(my_length/4):ceil(3*my_length/4)
        %   ondoleta(jC) = -1*ondoletaPirata(iC);
        %   jC=jC+1;
        % end

        figure(16);
        subplot(2,1,1);
        plot (vectorV); title('Senal ECG');
        grid on;
        subplot(2,1,2);
        plot (ondoleta); title('Wavelet');
        grid on;
 

        q_on = min(ondoleta);
        q_on25 = (q_on*0.25);

        s_off = max(ondoleta);
        % s_off25 = s_off+s_off *0.25;
        s_off25 = s_off*0.25;

        % calcula inicio de onda Q
        for loc=1:length(ondoleta)
            if q_on == ondoleta(loc)
                break;
            end
        end
 
        for loc_min=loc:-1:1
            if ondoleta(loc_min)< q_on25
                continue;
            else
                break;
            end
        end 
 
        q_on_loc(i)=R_loc(i)-loc_min+8; % el 7 es el ajuste por desplazamiento
        if q_on_loc(i) < 0 
            q_on_loc(i) = 1;
        end

        % calcula fin de onda S
        for loc=loc_min:length(ondoleta)
            if s_off == ondoleta(loc)
                break;
            end
        end
 
        for loc_max=loc:length(ondoleta)
            % if ondoleta(loc_max) > 0
            if ondoleta(loc_max) > s_off25
                continue;
            else
                break;
            end
        end
 
        s_off_loc(i) = R_loc(i)-offsetVentanaQRS+loc_max+7; % 7 es el ajuste del desplazamiento

        anchoQRS(i) = t(s_off_loc(i))-t(q_on_loc(i));
 
        % a cero QRS
        % for kk=Q_loc(i): S_loc(i)
        for kk=R_loc(i)-30: R_loc(i)+57
            x2(kk)=0;
        end
 
        % figure(100)
        % plot ([0:length(x2)-1]/fs,x2); 
        % title('Complejo QRS a 0')
        % grid on;
        % xlabel('Segundos');ylabel('mVolts')
        % xlim([0 2])

        inicioP = 0;
        
        %% DETECCIÓN ONDA P
        
        % Ventana búsqueda P
        msBusqP = 200;
        offsetBusqP = round( (msBusqP * fs )/1000 );
        msBusqP2 = 50;
        offsetBusqP2 = round( (msBusqP2 * fs )/1000 );
        
        if (R_loc(i) - offsetBusqP) > offsetBusqP2  % se usa 90
            % vectorP = x2(R_loc(i)-72:q_on_loc(i));
            vectorP = x2(R_loc(i)-offsetBusqP:q_on_loc(i));
            % ondoletaP=cwt(vectorP,6,'gaus1');
 
            ondoletaP = calculaOndoleta(vectorP,6,0);
            ondoletaP = ondoletaP(1:length(vectorP),1);

            max_P = max(ondoletaP);
            max_P05 = max_P*0.95;
  
            % encuentra valor minimo en la primer mitad de la ondoleta
            min_P = ondoletaP(1);
            for indice=2:length(ondoletaP)/2
                if ondoletaP(indice) < min_P
                    min_P = ondoletaP(indice);
                end
            end
            min_P05 = min_P*0.95;
  
            % figure(12);
            % subplot(2,1,1);
            % plot (vectorP); title('Ventana P');
            % grid on;
            % subplot(2,1,2);
            % plot (ondoletaP); title('Ondoleta P');
            % grid on;
  
            % encuentra punto minimo
            for iP = 1:length(ondoletaP)/2
                if min_P == ondoletaP(iP)
                    break;
                end
            end
            for iPCero = iP:length(ondoletaP)
                if ondoletaP(iPCero) > 0
                    break;
                end
            end
 
            % calcula inicio onda P
            for loc_minP=iP:-1:1
                if min_P05 > ondoletaP(loc_minP)
                    continue;
                else
                    break;
                end
            end
            iOndaP = iOndaP+1;

            % picoP(iOndaP) = R_loc(i)-72+iPCero-6;
            % p_inicio(iOndaP) = R_loc(i)-72+loc_minP-6;
            picoP(iOndaP) = R_loc(i)-offsetBusqP+iPCero;
            p_inicio(iOndaP) = R_loc(i)-offsetBusqP+loc_minP-6;
            inicioP =  p_inicio(iOndaP);
            %encuentra punto maximo
            for iP = iPCero:length(ondoletaP)
                if max_P == ondoletaP(iP)
                    break;
                end
            end

            for loc_maxP=iP:length(ondoletaP)
                if max_P05 < ondoletaP(loc_maxP) 
                    continue;
                else
                    break;
                end
            end
            % p_final(iOndaP) = R_loc(i)-90+loc_maxP+3;
            % p_final(iOndaP) = R_loc(i)-90+iP-3;
            % p_final(iOndaP) = R_loc(i)-72+iP-6;
            p_final(iOndaP) = R_loc(i)-offsetBusqP+iP;
            if p_final(iOndaP) > q_on_loc(i)
                p_final(iOndaP) =  q_on_loc(i)-6;
            end
            duracionP(iOndaP) = t(p_final(iOndaP)) - t(p_inicio(iOndaP));
            intervaloPR(iOndaP) = t(q_on_loc(i))-t(p_inicio(iOndaP));

            % a cero onda P
            for kk=p_inicio(iOndaP): p_final(iOndaP)
                x2(kk)=0;
            end
            
            % figure(101)
            % plot ([0:length(x2)-1]/fs,x2); 
            % title('P a 0')
            % grid on;
            % xlabel('Segundos');ylabel('mVolts')
            % xlim([0 2])
        end
        
        % TERMINA ONDA P

        %% DETECCIÓN ONDA T
        
        msVentanaDetecT = round((0.4 * mean(RR))*1000);
        % mean(RR) está en segs, se pasa a ms
%         msVentanaDetecT = 300;
        % 480ms según reporte
        % ~420ms según valor de 150 muestras con fs 360
        % 300 ms
        % 0.4 * RR, según literatura
        ventanaDetecT = round( (msVentanaDetecT*fs)/1000 );
  
        if ((s_off_loc(i) + ventanaDetecT) < longECG )
            vectorOT = x2(s_off_loc(i):s_off_loc(i)+ventanaDetecT);
            % ondoletaT=cwt(vectorOT,10,'gaus1');
    
            ondoletaT = calculaOndoleta(vectorOT,6,0);
            ondoletaT = ondoletaT(1:length(vectorOT),1);
    
            % min_T = min(ondoletaT);
            % min_T05 = min_T*0.95;

            % encuentra valor minimo en la primer mitad de la ondoleta
            min_T = min(ondoletaT);
            for indice=2:length(ondoletaT)/2
                if ondoletaT(indice) < min_T
                    min_T = ondoletaT(indice);
                end
            end
            min_T05 = min_T*0.05;
  
            % figure(13);
            % subplot(2,1,1);
            % plot (vectorOT); title('Ventana T');
            % grid on;
            % subplot(2,1,2);
            % plot (ondoletaT); title('Ondoleta T');
            % grid on;
            %   
            % max_T = max(ondoletaT);
            % max_T05 = max_T*0.95;
            % encuentra punto minimo
            for iT = 1:length(ondoletaT)
                if min_T == ondoletaT(iT)
                    break;
                end
            end
            for iTCero = iT:length(ondoletaT)
                if ondoletaT(iTCero) > 0
                    break;
                end
            end
 
            % calcula inicio onda T
            for loc_minT=iT:-1:1
                if min_T05 > ondoletaT(loc_minT)
                    continue;
                else
                    break;
                end
            end
            iOndaT = iOndaT+1;
            picoT(iOndaT) = s_off_loc(i)+iTCero-5;
            % t_inicio(iOndaT) = s_off_loc(i)+loc_minT-6;
            t_inicio(iOndaT) = s_off_loc(i)+loc_minT-5;
    
            % encuentra punto maximo
            max_T = 0;
            for iTT = iTCero:length(ondoletaT)
                if max_T >= ondoletaT(iTT)
                    break;
                else
                    max_T = ondoletaT(iTT);
                end
            end
            max_T05 = max_T *0.95;
            for loc_maxT=iTT:length(ondoletaT)
                if max_T05 < ondoletaT(loc_maxT)
                    continue;
                else
                    break;
                end
            end
 
            % t_final(iOndaT) = s_off_loc(i)+loc_maxT-6;
            t_final(iOndaT) = s_off_loc(i)+loc_maxT+8;
            duracionT(iOndaT) = t(t_final(iOndaT)) - t(t_inicio(iOndaT));
    
            % if i == 2
            %     ex=fopen('vectorOT.txt','w');
            %     for kk=1:length(vectorOT)
            %            fprintf(ex,'%d,',vectorOT(kk));
            %     end
            %     fclose(ex);
            % end
  
            segmentoST(i) = t(t_inicio(iOndaT))-t(s_off_loc(i));
            intervaloQT(i) = t(t_final(iOndaT))-t(q_on_loc(i));

        end
        
        % FIN DETECCIÓN T
    end 
end

%% Graficación P, QRS y T

%display(long_R);

% hace la ventana
%grafX1=[valorX1 valorX1];
%grafY1 = [0 0.5];

%valorX2=t(p_inicio(1));
%grafX2=[valorX2 valorX2];
%grafY2 = [0 0.5];

%grafX3 = [valorX1 valorX2];
%grafY3 = [0.5 0.5];

figure(10)
%subplot(3,1,1)
%complejo qrs
plot ([0:length(x22)-1]/fs,x22 , t(R_loc) ,R_value , 'r^', t(q_on_loc) ,x22(q_on_loc),'+',t(s_off_loc),x22(s_off_loc),'x'); 
title('Complejo QRS')
legend('ECG','R','Q-on','S-off');
grid on;
xlabel('Segundos');ylabel('mVolts')
%xlim([274 276.5])
xlim([0 3])

figure(11)
%subplot(3,1,2)
%onda p
plot ([0:length(x22)-1]/fs,x22,  t(picoP) ,x22(picoP) , '+',t(p_inicio),x22(p_inicio),'o',t(p_final),x22(p_final),'*');
%plot (t,x1 ,  t(picoP) ,x1(picoP) , 'o', t(picoT) ,x1(picoT) , '*');%onda p
title('Onda P')
%title('Onda P y Onda T')
legend('ECG','P','Pi','Pf');
%legend('ECG','P','T');
grid on;
%xlabel('Segundos');ylabel('mVolts')
%xlim([274 276.5])
xlim([0 3])

%subplot(3,1,3)
figure(122)
%onda t
plot ([0:length(x22)-1]/fs,x22 ,  t(picoT) ,x22(picoT) , '+',t(t_inicio),x22(t_inicio),'o',t(t_final),x22(t_final),'*');
title('Onda T')
legend('ECG','T','Ti','Tf');
grid on;
%xlabel('Segundos');ylabel('mVolts')
%xlim([274 276.5])
xlim([0 3])

%% Cálculos

%calculo duracion intervalos P-P
intervaloPP = 1:length(picoP)-1;
for k=1:length(picoP)-1
    intervaloPP(k) = t(picoP(k+1))-t(picoP(k));
end

% calculo promedio de amplitudes, intervalos y segmentos
picoPromOndaP = mean(x22(picoP));
picoOndaP = x22(picoP);
picoOndaR = x22(R_loc);
picoOndaS = x22(S_loc);
picoOndaQ = x22(Q_loc);
picoOndaT = x22(picoT);

duracionPromOndaP = mean(duracionP);
amplitudPromOndaT = mean(x22(picoT));
duracionPromOndaT = mean(duracionT);
duracionPromIntervaloRR = mean(RR);
duracionPromIntervaloPR = mean(intervaloPR);
duracionPromComplejoQRS = mean(anchoQRS);
duracionPromSegmentoST = mean(segmentoST);
duracionPromIntervaloQT = mean(intervaloQT);

%verifica frecuencias cardiacas
frecuenciaCardiacaVentricular = 60/mean(RR);
frecuenciaCardiacaAuricular = 60/mean(intervaloPP);
relacionAuricularVentricular = round(frecuenciaCardiacaAuricular)/round(frecuenciaCardiacaVentricular);
if relacionAuricularVentricular > 1
    relacionFrecuencia = 1; %frecuencia auricular es mayor que la ventricular
else
    if relacionAuricularVentricular < 1
        relacionFrecuencia = -1; %frecuencia auricular es menor que la ventricular
    else
        relacionFrecuencia = 0; %relacion frecuencia es 1:1
    end
end
    
%display(relacionFrecuencia);

% verifica distancia constante de intervalos PP
% el error debe ser menor al 10% del promedio de los intervalos PP para
% que sea regular = 1 de otra manera es irregular = 0
promPP = mean(intervaloPP);
for ipp = 1:length(intervaloPP)
    diferencia = promPP - intervaloPP(ipp);
    % if diferencia > promPP*0.10;
    if diferencia < 0.04
        ritmoAuricular(ipp) = 1;
    else
        ritmoAuricular(ipp) = 0;
    end
end
tasaAuricular = promPP;

% verifica distancia constante de intervalos RR
% el error debe ser menor al 10% del promedio de los intervalos RR para
% que sea regular = 1 de otra manera es irregular = 0
for ipp = 1:length(RR)
    diferencia=duracionPromIntervaloRR - RR(ipp);
    %if diferencia > duracionPromIntervaloRR*0.10;
    if diferencia < 0.04
        ritmoVentricular(ipp) = 1;
    else
        ritmoVentricular(ipp) = 0;
    end
end

tasaVentricular = duracionPromIntervaloRR;

% Complejo QRS, ancho normal de QRS entre 60 y 120 mseg.
for iqrs = 1:length(anchoQRS)
    if anchoQRS(iqrs) >= 0.060 && anchoQRS(iqrs) <= 0.120
        normalQRS(iqrs) = 1;
    else
        if anchoQRS(iqrs) < 0.060
            normalQRS(iqrs) = 0;
        else
            normalQRS(iqrs) = 2;
        end
    end
end

% Intervalo PR, ancho normal de PR entre 120 y 200 mseg.
for ipr = 1:length(intervaloPR)
    if intervaloPR(ipr) >= 0.120 && intervaloPR(ipr) <= 0.200
        normalPR(ipr) = 1;
    else if intervaloPR(ipr) < 0.120
        normalPR(ipr) = 0;
        else
             normalPR(ipr) = 2;
        end
    end
end

% Intervalo QT, ancho normal de QT entre (360 y 440 mseg) o 1/2 RR.
if length(RR) < length(intervaloQT)
    limSQT = length(RR);
else
    limSQT = length(intervaloQT);
end

normalQT = zeros(1,limSQT);
for ipr = 1:limSQT
    if (intervaloQT(ipr) >= 0.360 && intervaloQT(ipr) <= 0.440) || (intervaloQT(ipr) < 0.5*RR(ipr))
        normalQT(ipr) = 1;
    else
        if intervaloQT(ipr) < 0.360
            normalQT(ipr) = 0;
        else
            normalQT(ipr) = 2;
        end
    end
end

% Onda P, ancho normal de P menor de 100 mseg.(limite normal en personas
% adultas es 110-20 ms a 110+20 ms)
normalP = zeros(1,length(duracionP));
for ipr = 1:length(duracionP)
    if  duracionP(ipr) > 0.0 && duracionP(ipr) <= 0.130
        normalP(ipr) = 1;
    else
        if  duracionP(ipr) < 0.090
            normalP(ipr) = 0;
        else
            normalP(ipr) = 2;
        end
    end
end

% obtiene vector de caracteristicas
iC=1;
%for latido = 3:length(picoOndaR)-1
for latido = 3:length(picoOndaT)
    if latido <= length(RR)
        % ################### Intervalos RR ###################
        rr_interval_1 = RR(latido-1); % x1=lj -lj-1 rr
        rr_interval_2 = RR(latido-2); % x2=lj-1 - lj-2 pre
        rr_interval_3 = RR(latido); % x3=lj+1 -lj pos
    
        px4 = rr_interval_1 - rr_interval_2;  %  x1-x2
        px5 = rr_interval_3 - rr_interval_1; % x3-x1
    
        if length(picoOndaP) == length(picoOndaR)
            pico_p = picoOndaP(latido-1);
        else
            pico_p = picoOndaP(latido-2);
        end
    
        pico_q = picoOndaQ(latido-1);
        pico_r = picoOndaR(latido-1);
        pico_s = picoOndaS(latido-1);
        pico_t = picoOndaT(latido-1);
  
        rr_category = 1;
%        rr1 = rr_interval_1/fs;
%        rr2 = rr_interval_2/fs;
%        rr3 = rr_interval_3/fs;
        rr1 = rr_interval_1;
        rr2 = rr_interval_2;
        rr3 = rr_interval_3;

        if(rr2 < 0.6 && 1.8*rr2 < rr1)
            rr_category = 3;
            % TODO restante del algoritmo para la prueba con grabaciones reales
            % para correcta deteccion de episiodos de aleteo ventricular 
        end
    
        if( (rr1 > rr2*1.15 && rr3 > rr2*1.15) || ...
            (abs(rr1-rr2) > 0.3 && rr2 < 0.8 && rr1 < 0.8 && 1.2*((rr1+rr2)/2) < rr3) || ...
            (abs(rr2-rr3) > 0.3 && rr2 < 0.8 && rr3 < 0.8 && 1.2*((rr3+rr2)/2) < rr1) )
        
            rr_category = 2;
        end
    
        if((2.2 < rr2 && rr2 < 3) && (abs(rr1-rr2) < 0.2 || abs(rr2-rr3) < 0.2))
            rr_category = 4;
        end
        
%       prematurity =   (rr_interval_3/rr_interval_1)^2 + ...  
%                       (rr_interval_2/rr_interval_1)^2 - ...
%                       wentropy([rr_interval_1, rr_interval_2, rr_interval_3], 'shannon');
                
        prematurity =   (rr_interval_3/rr_interval_2)^2 + ...  % x6
                        (rr_interval_1/rr_interval_2)^2 - ...
                        0.333*(((rr_interval_1^2)*log(rr_interval_1^2))+((rr_interval_2^2)*log(rr_interval_2^2))+((rr_interval_3^2)*log(rr_interval_3^2)));
    
        % energia qrs
        energia_qrs = 0;
        puntosQRS = 1;
        existeQRS = 0;
        for e = q_on_loc(latido-1):s_off_loc(latido-1)
            energia_qrs = energia_qrs+abs(x22(e).^2);
            amplitudesQRS(puntosQRS) = x22(e);
            puntosQRS = puntosQRS +1;
            existeQRS = 1;
        end
    
        % polaridad complejo QRS
        if existeQRS == 1
            polaridad_qrs = abs(max(amplitudesQRS)/min(amplitudesQRS));
        else
            polaridad_qrs=0;
        end
    
        if polaridad_qrs > 0
            polaridad_qrs = 1;
        else
            polaridad_qrs = -1;
        end
   
        % muestras amplitud pqr y rst
%       pqr_amp = [];
%       rst_amp = [];
%       for i=1:14
%           if ((R_loc(latido-1)-i*7) < 1 )
%               pqr_amp(i) = 0.0;
%               rst_amp(i) = 0.0;
%               continue;
%           else
%               pqr_amp(i) = x22(R_loc(latido-1)-i*7);
%               rst_amp(i) = x22(R_loc(latido-1)+i*7);
%           end
%       end
    
        %morfologia en QRS (12 muestras)
        qrs_amp = [];
        iQRS=1;
        indQRS = [];
        for i= 1:12
            qrs_amp(iQRS) = x22((R_loc(latido-1)-20)+(i*5));
            indQRS(iQRS) = (R_loc(latido-1)-20)+(i*5);
            iQRS=iQRS+1;
        end
   
        %morfologia en T (9 muestras)
        t_amp = [];
        iQRS=1;
        rangoT = t_final(latido-1) - t_inicio(latido-1);
        pasoT = round(rangoT/9);
        indTT = [];
        for i=1:9
            t_amp(iQRS) = x22((t_inicio(latido-1))+(i*pasoT)); 
            indTT(iQRS) = (t_inicio(latido-1))+(i*pasoT);
            iQRS=iQRS+1;
        end
    
%       figure(12);
%       plot ([0:length(x22)-1]/fs,x22, t(indQRS) ,x22(indQRS) , '*',t(indTT) ,x22(indTT) , '+'); %complejo qrs
%       title('Ventana QRS');
%       grid on;
%       xlim([0 2])

        %onda P
        amplitudP =0;
        if p_inicio(1) < R_loc(1)
            duracionP = t(p_final(latido-1)) - t(p_inicio(latido-1));
            duracionPQ = t(q_on_loc(latido-1)) - t(p_inicio(latido-1));
        else
            duracionP = t(p_final(latido-1)) - t(p_inicio(latido-1));
            duracionPQ = t(q_on_loc(latido))- t(p_inicio(latido-1));
        end
    
        if pico_p > x22(p_final(latido-1)) && pico_p > x22(p_inicio(latido-1))
            amplitudP = 1;
        else
            amplitudP = 0;
        end
    
        if t(p_final(latido-1)) > t(Q_loc(latido-1))
            verificaExisteP = 0;
        else
            verificaExisteP = 1;
        end
    
        if duracionP > 0.0 && duracionP <= 0.120 && amplitudP == 1 && verificaExisteP == 1
            existeP = 1; % presencia P positiva
        else
            if duracionP > 0.100 && amplitudP == 1 && verificaExisteP == 1
                existeP = 1; % presencia
            else
                if amplitudP == 0 || verificaExisteP == 0
                    existeP = 0;  % ausencia
                end
            end
        end
    
        %complejo QRS
        duracionQRS = t(s_off_loc(latido-1)) - t(q_on_loc(latido-1)); ... % duracion complejo QRS
        if duracionQRS >= 0.060 && duracionQRS <= 0.120
            tipoQRS = 1;
            existeQRS = 1;
        else
            if duracionQRS < 0.060
                tipoQRS = 0;
                existeQRS=0;
            else
                tipoQRS = 2;
                existeQRS = 1;
            end
        end
    
        % onda T
        conT = t_inicio(latido-1);
        for ti = 1:t_final(latido-1)-t_inicio(latido-1)
            ondaTT(ti) = x22(conT);
            conT = conT+1;
        end
     
        duracionT = t(t_final(latido-1)) - t(t_inicio(latido-1));
   
%       figure(133)
%       plot (ondaTT);%onda t
%       title('Onda T')
%       grid on;
%       xlabel('Segundos');ylabel('mVolts')
%       xlim([274 276.5])

        if duracionT >= 0.120 && duracionT <= 0.160
            tipoT = 1;
        else
            if duracionT < 0.120
                tipoT = 0;
            else
                tipoT = 2;
            end
        end
     
        % intervalo QTc
        duracionQT = t(t_final(latido-1)) - t(q_on_loc(latido-1));
        duracionQTc = duracionQT / sqrt(rr3);
        if duracionQTc >= 0.340  && duracionQTc <= 0.450;
            tipoQTc = 1;
        else 
            tipoQTc = 0;
        end
     
        % PQ %%%%%%%%%%%
        if existeP == 1
            if duracionPQ >= 0.120 && duracionPQ <= 0.2
                tipoPQ = 1;
            else
                if duracionPQ < 0.120
                    tipoPQ = 0;
                else
                    tipoPQ = 2;
                end
            end
        else
            duracionPQ = 0;
            tipoPQ = 2;
        end
         
        % frecuencia
        if frecuenciaCardiacaVentricular < 60
            cf = 0; % bradicardia
        else
            if frecuenciaCardiacaVentricular >= 60 && frecuenciaCardiacaVentricular <= 100
                cf = 1; % normal
            else
                cf = 2; % taquicardia;
            end
        end
  
        % prueba si P es positivo y existe antes de QRS
        if pico_p > 0.0 && existeP == 1 && existeQRS == 1
            PantesQRS = 1;
        else
            PantesQRS = 0;
        end
    
        % intervalo RP %%%%%%%
        if existeP == 1
            intervaloRP = t(p_inicio(latido-1))-t(q_on_loc(latido-2));
        else
            intervaloRP = 0;  
            intervaloPR(latido-1)=0;
        end
    
    %################### Vector de Caracteristicas #################
    vector_caracteristicas(1:39,iC:iC) = [ 
%                                               pico_p; ... % amplitud P
%                                                pico_q; ... % amplitud Q
%                                                pico_r; ... % amplitud R
%                                                pico_s; ... % amplitud S
%                                                pico_t; ... % amplitud T
                                               %duracionP; ... % duracion onda P
                                             %  cf; %frecuencia
                                               existeP; ... % presencia/ausencia onda P
                                               logsig(duracionQRS); ... % duracion complejo QRS
                                               tipoQRS; ... %tipo QRS
                                               logsig(duracionPQ); ... % duracion intervalo PQ
                                               tipoPQ; ... %tipo PQ
                                               logsig(duracionQTc); ... %duracion QTc
                                               tipoQTc; ... %tipo QTc
                                               logsig(duracionT); ... % duracion onda T
                                               tipoT; ... %tipo T
                                               rr_interval_1; ... % intervalo RR
                                               rr_interval_2; ... % intervalo pre RR
                                               rr_interval_3; ... % intervalo post RR
                                               duracionPromIntervaloRR; % promedio local
                                               rr_category; ... % categoria derivada de los intervalos RR
                                               px4; ... % incremento o decremento tasa corazon
                                               px5; ... % incremento o decremento tasa corazon
                                               prematurity; ... % prematuridad
                                               logsig(energia_qrs); ... % energia QRS
                                           %    polaridad_qrs; ...  % polaridad complejo QRS
                                               tansig(qrs_amp(:)); ... % muestras en QRS
                                               tansig(t_amp(:));     % muestras en T
        ];
    
%     duracionQRS =duracionQRS * 100;
%     duracionQRS = floor(duracionQRS);
%     duracionQRS = duracionQRS * 0.01;
%     duracionPQ =duracionPQ * 100;
%     duracionPQ = floor(duracionPQ);
%     duracionPQ = duracionPQ * 0.01;
%     duracionQTc =duracionQTc * 100;
%     duracionQTc = floor(duracionQTc);
%     duracionQTc = duracionQTc * 0.01;
%     duracionQT =duracionQT * 100;
%     duracionQT = floor(duracionQT);
%     duracionQT = duracionQT * 0.01;
    vector_intervalos(1:12,iC:iC) = [ 
                            duracionQRS; ...
                            duracionPQ; ...
                            duracionQTc; ...
                            duracionQT; ...
                            cf; ...
                            frecuenciaCardiacaVentricular; ...
                            PantesQRS; ...
                            ritmoVentricular(latido-1); ...
                            intervaloPR(latido-1);
                            tasaAuricular;
                            tasaVentricular;
                            intervaloRP;
    ];
    
    iC = iC+1;
    end
end

columnas = iC-1;
else
   columnas = 0; 
   vector_caracteristicas = zeros(39);
   vector_intervalos = zeros(12);
end

%% PRUEBA VARIABLES GLOBALES
global pruebaTiempo;
pruebaTiempo= t;
global pruebaValores;
pruebaValores= x22;
global pruebaOndaSOff;
pruebaOndaSOff = s_off_loc;
global pruebaOndaTIni;
pruebaOndaTIni = t_inicio;
global pruebaOndaTFin;
pruebaOndaTFin = t_final;
global pPicoT;
pPicoT = picoT;

global pFinalP;
pFinalP = p_final;
global pInicioQ;
pInicioQ = q_on_loc;

global R_pos;
R_pos = R_loc;
global R_val;
R_val = R_value;

%% Prueba Características ST-T
% global caracteristicasSTT;
% caracteristicasSTT = zeros(5,5);

end

