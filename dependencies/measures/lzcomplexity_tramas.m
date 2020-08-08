function [LZC,meanLZC]=lzcomplexity_tramas(registro,umbral_db,num_simbolos,tam,solapamiento)
%   function LZcomplexity=lzcomplexity_tramas(registro,umbral_db,tam,solapamiento)
%
%   La serie de entrada es dividida en segmentos, de forma que se obtiene 
%   un valor de LZC para cada segmento. Finalmente se realiza un 
%   promediado para obtener un único valor de LZC para la serie.
%%
%   Argumentos:
%
%   registro:       Serie de datos de entrada de la que estimaremos su LZC.
%   umbral_db:      Umbral empleado en la conversión binaria de la serie.
%                   Tomará los valores 'media' o 'mediana'.
%   num_simbolos:   Número de símbolos empleados en la conversión binaria
%                   de la serie de datos original.
%   tam:            Número de muestras de los segmentos en que se dividirá la 
%                   serie de datos original. Si tam = 0, se aplicará el
%                   algoritmo sobre la serie original sin segmentar.
%   solapamiento:   valor entre 0 y 1 que indica el tanto por uno de solapamiento
%
%   Variables de salida:
%   LZC:        Valor de LZC de la serie de datos de entrada.
%
%
%   Última actualización: 17/01/2013
%   Autor: Carlos Gómez Peña
%

% Calculamos el número de veces a aplicar el algoritmo según el número de 
% segmentos de longitud tam muestras en los que se puede dividir la serie 
% de datos de entrada.

warning off

if tam==0   % Si tam = 0 el algoritmo se calcula sobre el registro completo
    n_tramas=1;
    tam=length(registro);
    solapamiento=0; % Independientemente de lo que haya introducido el usuario
else
    n_tramas=fix((length(registro)-tam)/(tam*(1-solapamiento))+1);
end

% Calculamos LZC para cada trama
for l=1:n_tramas
    % Extraemos la trama correspondiente del registro original.
    if l==n_tramas
        trama=registro((1+(l-1)*tam*(1-solapamiento)):end);
    else
        trama=registro((1+(l-1)*tam*(1-solapamiento)):((l-1)*tam*(1-solapamiento)+tam));
    end
    
    % Número de muestras del segmento.
    n=length(trama);
    
    % Transformación de la trama en una secuencia binaria
    trama=transforma_binaria(trama,num_simbolos,umbral_db);
    
    if num_simbolos == 2
        b=n/log2(n);
    else
        b=n/(log(n)/log(3));
    end

    % Inicializamos las variables empleadas en el cálculo la complejidad LZ.
    c = 1;          % Valor inicial del contador de complejidad.
    S = trama(1);   % Inicialización de la subsecuencia S.
    Q = trama(2);   % Inicialización de la subsecuencia Q.
    
    for i = 2:tam
        % Concatenamos ambas subsecuencias.
        SQ = [S,Q];
        % Eliminamos el último caracter de la subsecuencia resultado de la
        % concatenación.
        SQ_pi = [SQ(1:(length(SQ)-1))];
        
        % Comprobamos si la subsecuencia Q se encuentra contenida en SQ_pi
        indice = findstr(Q,SQ_pi); % Nos da los índices en los que Q empieza dentro de SQ_pi
        
        if length(indice)==0
            % Q no se encuentra al inspeccionar SQ_pi: Q es una secuencia nueva.
            c = c+1;                % Incrementamos el contador de complejidad.
            if (i+1)>tam            % Si llegamos al final de la serie hemos terminado.
                break;
            else                    % Si no hemos llegado al final concatenamos las subsecuencias S y Q.
                S = [S,Q];          % Formamos una nueva subsecuencia S.
                Q = trama(i+1);     % Actualizamos la subsecuencia Q.
            end
        else
            % Q forma parte de SQ_pi.
            if (i+1)>tam            % Si llegamos al final de la serie hemos terminado.
                break;
            else
                Q = [Q,trama(i+1)];	% Extendemos la subsecuencia Q.
            end
        end
    end
    % Normalizamos el contador de complejidad c de tal forma que 0 <= c/b <= 1
    LZC(l) = c/b;
    meanLZC = mean(LZC);
end

