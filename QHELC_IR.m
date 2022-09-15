function image = QHELC_IR(f,gamma)
    
    [count,binLocation]=imhist(f);
    histo=[binLocation,count];
    
    [cantidad, clasificacion] = histcounts(f(:), 'Normalization', 'probability','NumBins',256 );
    clasificacion(end)=[];
    figure; plot(count, binLocation);
    figure; plot( binLocation, count);
%    SE HALLA EL PRIMER THRESHOLD

    transpuesto=transpose([clasificacion;cantidad]);
    [m,n]=size(transpuesto);
    promedio=average(histo,1,m);
    centroide=inicializarCentroide(promedio);
    [test,c,sum,D]=kmeans(transpuesto,2,'start',centroide);
    SP=spFinder(test);
        

%   SE APLICA KMEANS A LA PRIMERA MITAD
    promedio1=average(histo,1,SP);
    centroide=inicializarCentroide(promedio1);
    [test1,c1,sum1,D1]=kmeans(transpuesto(1:SP,:),2,'start',centroide);
    SPL=spFinder(test1);
    
    
% SE APLICA KMEANS A LA SEGUNDA MITAD
    promedio2=average(histo,SP+1,m);
    centroide=inicializarCentroide(promedio2);
    [test2,c2,sum2,D2]=kmeans(transpuesto(SP+1:m,:),2,'start',centroide);
    SPU=spFinder(test2)+SP;
    
    
% SE HALLAN LOS CUTOFF LIMITS

    CL1 = cutOffLimit(histo, 1, SPL, gamma);
    CL2 = cutOffLimit(histo, SPL+1, SP, gamma);
    CU1 = cutOffLimit(histo, SP+1, SPU, gamma);
    CU2 = cutOffLimit(histo, SPU+1, 256, gamma);

%  SE CALCULAN LOS PIXELES QUE ESTAN ENCIMA DE LOS CUTOFF POINTS

    TL1 = aboveCutOff(histo, 1, SPL, CL1);
    TL2 = aboveCutOff(histo, SPL+1 , SP, CL2);
    TU1 = aboveCutOff(histo, SP+1, SPU, CU1);
    TU2 = aboveCutOff(histo, SPU+1, 256, CU2);
    
    
%   SE CALCULA EL AVERAGE INCREMENT

    AL1 = averageIncrement(TL1, 1, SPL);
    AL2 = averageIncrement(TL2, SPL+1, SP);
    AU1 = averageIncrement(TU1, SP+1, SPU);
    AU2 = averageIncrement(TU2, SPU+1, 256);
    
%   SE SACAN LOS PIXELES QUE ESTAN ENCIMA DE LOS CUTOFF POINTS
    
    hmod = trimmedH(histo, 1, SPL, CL1, AL1);
    hmod = trimmedH(hmod, SPL+1, SP, CL2, AL2);
    hmod = trimmedH(hmod, SP+1, SPU, CU1, AU1);
    hmod = trimmedH(hmod, SPU+1, 256, CU2, AU2);
    

    
%   EQUALIZACION
    
    equalizado= equalizacion(hmod, SP, SPL, SPU);
    

%   MODIFICACION DE IMAGEN

    [m,n]=size(f);
    image=f;
    for i=1:m
        for j=1:n
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
    x=0;
    y=0;
    for i=inicio:fin
        x=x+histograma(i,1)*histograma(i,2);
        y=y+histograma(i,2);
    end
    promedio=x/y;
    
end

function centroide= inicializarCentroide(promedio)
    matriz=zeros(2,2);

    for i=1:2
        matriz(i)=ceil(promedio);
        matriz(i,2)=ceil(promedio);
    end
    centroide=matriz;
end

function SP= spFinder(test)
    i=1;
    while test(i)==test(1)
        SP=i;
        i=i+1;
    end
end

function CI = cutOffLimit(histo,LL,UP,gamma)
    N=0;
    for i=LL:UP
        N= histo(i,2 ) +N;
    end
    
    I1 = UP - LL;
    CI= ceil(N/I1) + round(gamma * (N - (N/I1)));

end

function T = aboveCutOff(histo,LL,UP,CL)
T=0;    
    for i=LL:UP
        T = T + max(histo(i,2)-CL,0);
    end
end

function A = averageIncrement(T,   LL , UP)
    I = UP - LL;
    A = floor(T/I);
end


function h = trimmedH (histo, LL, UP, CL, AI)
    for i = LL:UP
        if histo(i,2) > (CL-AI)
            histo(i,2) = CL;
            
        else
            histo(i,2)= histo(i,2)+AI;
        end
    end
    h=histo;
end


function e = equalizacion(hmod, SP, SPL, SPU)
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
    a=0;
    for i=LL:q
        a= a + pdf(i, histo, LL, UP);
    end
end

function p=pdf(q, histo, LL, UP)
    N=0;
    
    for i=LL:UP
        
        N= histo(i,2 ) +N;
    end
    
    p=histo(q,2)/N;
end

function  graficarHistograma(histo, SP, SPL, SPU, CL1, CL2, CU1, CU2, titulo)
    
    plot(histo)
    line([SP SP],[0 3000],'color','red')
    text(SP,0,'SP')
    line([SPL SPL],[0 3000],'color','green')
    text(SPL,0,'SPL')
    line([SPU SPU],[0 3000],'color','blue')
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
