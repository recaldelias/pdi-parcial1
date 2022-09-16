function image = QHELC_IR(f,gamma)  %                               PASO 1
%                                                                   PASO 2
    [count,binLocation]=imhist(f);
    histo=[binLocation, count ];
    pK = zeros(256,2);
    [x,y]=size(f);
    
    for i=1:256
        pK(i,1)=i-1;
        pK(i,2)=histo(i,2)/(x*y);
    end
    
%    SE HACE LA PRIMERA DIVISION                                    PASO 3

    
    [m,~]=size(histo);
    promedio=average(histo,1,m);
    centroide=inicializarCentroide(promedio);
    %[test,c,sum,D]=kmeans(histo,2,'start',centroide);
    [test,~,~,~]=kmeans(pK,2,'start',centroide);
    SP=spFinder(test);
    
%   SE APLICA KMEANS SOBRE LA PRIMERA MITAD                         PASO 4
    
    promedio1=average(histo,1,SP);
    centroide=inicializarCentroide(promedio1);
    [test1,~,~,~]=kmeans(pK(1:SP,:),2,'start',centroide);
    SPL=spFinder(test1);

%   SE APLICA KMEANS A LA SEGUNDA MITAD                             PASO 5 
    promedio2=average(histo,SP+1,m);
    centroide=inicializarCentroide(promedio2);
    [test2,~,~,~]=kmeans(pK(SP+1:m,:),2,'start',centroide);
    SPU=spFinder(test2)+SP;    

%   SE HALLAN LOS CUTOFF LIMITS                                     PASO 5.1

    CL1 = cutOffLimit(histo, 1, SPL, gamma);
    CL2 = cutOffLimit(histo, SPL+1, SP, gamma);
    CU1 = cutOffLimit(histo, SP+1, SPU, gamma);
    CU2 = cutOffLimit(histo, SPU+1, 256, gamma);

% SE CALCULAN LOS PIXELES QUE ESTAN ENCIMA DE LOS CUTOFF POINTS    PASO 5.2

    TI = aboveCutOff(histo, SP, SPL, SPU, CL1, CL2, CU1, CU2);

    
    
%   SE CALCULA EL AVERAGE INCREMENT                                 PASO 5.3

    AI = averageIncrement( TI, SP, SPL, SPU);
    
    %[AI,TI]
    
%   SE SACAN LOS PIXELES QUE ESTAN ENCIMA DE LOS CUTOFF POINTS      PASO 6
    
    hmod = trimmedH(histo, 1, SPL, CL1, AI);
    hmod = trimmedH(hmod, SPL+1, SP, CL2, AI);
    hmod = trimmedH(hmod, SP+1, SPU, CU1, AI);
    hmod = trimmedH(hmod, SPU+1, 256, CU2, AI);
    

    
%   EQUALIZACION
    
    equalizado= equalizacion(hmod, SP, SPL, SPU);
    

%   MODIFICACION DE IMAGEN                                          PASO 7

%   se modifican los niveles de gris por los hallados en el vector
%   equalizacion
    image=f;
    for i=1:x
        for j=1:y
            if equalizado(image(i,j)+1,2)~=0
                image(i,j)=equalizado(image(i,j)+1,2);
            end
        end
    end
    
    figure; graficarHistograma(histo,SP, SPL, SPU, CL1, CL2, CU1, CU2, 'original');
    figure; graficarHistograma(hmod,SP, SPL, SPU, CL1, CL2, CU1, CU2, 'modificado');
    figure; imshow(f); title('original');
    figure; imshow(image); title('modificado');
end
    
function promedio= average(histograma,inicio,fin)
%   funcion que halla el promedio de un histograma
    x=0;
    y=0;
    for i=inicio:fin
        x=x+histograma(i,1)*histograma(i,2);
        y=y+histograma(i,2);
    end
    promedio=x/y;
    
end

function centroide= inicializarCentroide(promedio)
%   funcion que inicializa la matriz de centroides para utilizar en
%   en algoritmo de kmeans
    matriz=zeros(2,2);

    for i=1:2
        matriz(i)=ceil(promedio);
        matriz(i,2)=ceil(promedio);
    end
    centroide=matriz;
end

function SP= spFinder(test)
%   funcion para hallar las divisiones entre los clusters
    i=1;
    while test(i)==test(1)
        SP=i;
        i=i+1;
    end
end

function CI = cutOffLimit(histo,LL,UP,gamma)
%   Establece los limites del histograma
    N=0;
    for i=LL:UP
        N= histo(i,2 ) +N;
    end
    
    I1 = UP - LL;
    CI= ceil(N/I1) + round(gamma * (N - (N/I1)));

end

function T = aboveCutOff(histo, SP, SPL, SPU, CL1, CL2, CU1, CU2)
%   Encuentra la cantidad de pixeles sobre el cutoff limit
    T=zeros(256,1);    
    for i=1:256
        if i>=0 && i<=SPL

            T(i,1) = max(histo(i,2)-CL1,0);
        elseif i>=(SPL+1) && i<=SP

            T(i,1) = max(histo(i,2)-CL2,0);
        elseif i>=(SP+1) && i<=SPU
            
            T(i,1) = max(histo(i,2)-CU1,0);
        elseif i>=(SPU+1) && i<=256
            
            T(i,1) = max(histo(i,2)-CU2,0);
        end
    end
end




function A = averageIncrement( TI, SP, SPL, SPU)
%   calcula el incremento promedio de la intensidad de los distintos
%   niveles de gris
    A=zeros(256,1);

    for i=1:256
        if i>=0 && i<=SPL
            
            longitud = SPL;
            A(i,1) = floor(TI(i,1)/longitud);
        elseif i>=(SPL+1) && i<=SP

            longitud = SP -SPL -1;
            A(i,1) = floor(TI(i,1)/longitud);
        elseif i>=(SP+1) && i<=SPU
            
            longitud = SPU -SP -1;
            A(i,1) = floor(TI(i,1)/longitud);
        elseif i>=(SPU+1) && i<=256
            
            longitud = 255 -SPU -1;
            A(i,1) = floor(TI(i,1)/longitud);
        end
    end
end


function h = trimmedH (histo, LL, UP, CL, AI)
%   funcion que modifica los distintos histogramas retirando los pixeles
%   que esten encima de la barrera del cutoff limit
    for i = LL:UP
        if histo(i,2) > (CL-AI(i,1))
            histo(i,2) = CL;
            
        else
            histo(i,2)= histo(i,2)+AI(i,1);

        end
    end
    h=histo;
end


function e = equalizacion(hmod, SP, SPL, SPU)
%   funcion que adapta el histograma para modificar la imagen
    e=zeros(256,2);
    for i=1:256
        e(i,1) = i-1;
    end
    
    for i=1:256
        
        if e(i,1)>=0 && e(i,1)<=SPL

            e(i,2)= 0 + (SPL - 0)*acumulativo(e(i,1),hmod, 1, SPL); 
        elseif e(i,1)>=(SPL+1) && e(i,1)<=SP

            e(i,2)= SPL+1 + (SP - SPL -1)*acumulativo(e(i,1),hmod, SPL+1, SP); 
        elseif e(i,1)>=(SP+1) && e(i,1)<=SPU
            
            e(i,2)= SP+1 + (SPU - SP - 1)*acumulativo(e(i,1),hmod, SP+1, SPU);
        elseif e(i,1)>=(SPU+1) && e(i,1)<=256
            
            e(i,2)= SPU+1 + (255 - SPU - 1)*acumulativo(e(i,1),hmod, SPU+1, 256); 

        end
    end
end


function a= acumulativo(q, histo, LL, UP)
%   funcion C(q) en el paper
    a=0;
    for i=LL:q
        a= a + pdf(i, histo, LL, UP);
    end
end

function p=pdf(q, histo, LL, UP)
% funcion que halla el probability density function de un histograma
    N=0;
    
    for i=LL:UP
        
        N= histo(i,2 ) +N;
    end
    
    p=histo(q,2)/N;
end

function  graficarHistograma(histo, SP, SPL, SPU, CL1, CL2, CU1, CU2, titulo)
%   funcion para graficar el histograma
    
    plot(histo(:,1),histo(:,2))
    %line([SP SP],[0 3000],'color','red')
    xline(SP,'color','red')
    text(SP,0,'SP')
    xline(SPL,'color','green')
    text(SPL,0,'SPL')
    xline(SPU,'color','blue')
    text(SPU,0,'SPU')
    text(256,0,'L')
    line([0 SPL],[CL1 CL1],'color','magenta')
    text(0,CL1,'CL1')
    line([SPL+1 SP],[CL2 CL2],'color','magenta')
    text(SPL,CL2,'CL2')
    line([SP+1 SPU],[CU1 CU1],'color','magenta')
    text(SP+1,CU1,'CU1')
    line([SPU+1 256],[CU2 CU2],'color','magenta')
    text(SPU+1,CU2,'CU2')
    title(titulo)
end