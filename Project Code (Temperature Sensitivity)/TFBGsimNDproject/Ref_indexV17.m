function n = Ref_indexV17(material,x,lambda,thermoOpticCoeff,temperature,count)
%refindex  Retourne l'indice de réfraction du matériau spécifié à la longueur d'onde lambda.
% Le matériau est une "string" choisie parmi les choix suivants:
%  1) SiO2-GeO2
%	   Dans ce cas, l'appel de la fonction a la forme refindex('SiO2-GeO2',x,lambda),
%		où x est la concentration molaire du dopant ( i.e. x=0.05 for 5% GeO2)
% La longueur d'onde lambda doit être entrée en microns.


switch material
    
    case 'MIP'
        n = input('What is the refractive index you want to use for the MIP layer?   ');
    
    case 'SRI'
        n = input('What is the refractive index you want to use for the outside medium?   ');
        
    case 'Fixed'
        n = input('What is the refractive index you want to use for the "Fixed" layer?   ');

	case 'SiO2-GeO2'
        AS1 = 0.69616630;
        AS2 = 0.40794260;
        AS3 = 0.89747940;
        lS1 = 0.068404300;
        lS2 = 0.11624140;
        lS3 = 9.8961610;
   
        AG1 = 0.80686642;
        AG2 = 0.71815848;
        AG3 = 0.85416831;
        lG1 = 0.068972606;
        lG2 = 0.15396605;
        lG3 = 11.841931;
   
        n = sqrt(1+(AS1+x*(AG1-AS1))./(1-(lS1+x*(lG1-lS1))^2./lambda.^2)+(AS2+x*(AG2-AS2))./(1-(lS2+x*(lG2-lS2))^2./lambda.^2)+(AS3+x*(AG3-AS3))./(1-(lS3+x*(lG3-lS3))^2./lambda.^2));
        
        %modify refractive index from baseline temp (chosen to be 0 - only care about changes in temperature)
        n = vpa(n + thermoOpticCoeff*temperature(count));
        n = double(n);
        
    case 'Air'
        n = 1.00027;
        
    case 'Water1_315'
        n = 1.315;
        
    case 'Gold'                         %Linear approximation valid for 1.25 to 1.65 band from FIMMWAVE values
        rn = 0.14206*lambda^2+0.22151*lambda-0.12745;
        in = -1*(5.97171*lambda^2-11.515*lambda+13.25);
        n = complex(rn,in);%NOTE: loss corresponds to negative imaginary part of the index
        
    case 'AuNP'
        fill=0.20;
        rn = -0.40794*fill+0.57044;
        in = 7.48422*fill-9.8966;
        n = complex(rn,in);%NOTE: loss corresponds to negative imaginary part of the index
        
    case 'CVD_Au'
        t=25;  %film thickness in nm (from P data of Wenjun JCP paper)
        rn = 1.53+0.0222*t;
        in = -(0.07+0.002);
        n = complex(rn,in);%NOTE: loss corresponds to negative imaginary part of the index
    
    case 'Pd'
        rn = 0.33071*lambda^2-0.04617*lambda+2.22709;
        in = -1*(-0.03929*lambda^2+4.38603*lambda+1.62122);
        n = complex(rn,in);%NOTE: loss corresponds to negative imaginary part of the index
        
    case 'Lossy'
        n = 1.00 -1i;
        
    case 'Water'                        %Linear approximation valid for C+L band only
        rn = -0.021*lambda+1.3505;
        in = 0.0002*lambda-0.0004;
        n = complex(rn,in);
        
    case 'Buffer'
        n = 1.35;
        
    case 'Poly'
        n=1.5;
    
    otherwise
        n = 5;
end
   
end

