function retVal = beta ( u , v , kernel, a, b)

switch kernel
    
case 1
    
    retVal	= (u .* v);                                                     % Produktkern
    
case 2
    
    retVal  = (u + v).^a./(u.*v).^b;                                        % Kapurkern
    
case 3
    
    retVal  = ( u.^(1/3) + v.^(1/3) ) .* ( u.^(-1/3) + v.^(-1/3) );         % Brownsche Bewegung
    
case 4
    
    retVal  =  ( u.^(1/3) + v.^(1/3) ).^3;                                  % Shear-Kern
    
case 5
    
    retVal  = ( u.^(1/3) + v.^(1/3) ).^2 .* abs(u.^(-1) - v.^(-1));         % EKE-Kern (hnlich)
    
case 6
    
    retVal  = ( u + v );                                                    % Summenkern
    
case 7
    
    retVal = 1;                                                             % grenunabhngiger Kern
    
case 8
        
	retVal = ( u.^(1/3) + v.^(1/3) ).^2 .* sqrt(u.^(-1) + v.^(-1));         % korrekter EKE
        
case 9
  
	retVal   = ( ( u.^(1/3) + v.^(1/3) ).^2 ) .* abs(u.^(1/6) - v.^(1/6));  % Gravitationskern
    
case 10
  
	retVal   = ( ( u.^(1/3) + v.^(1/3) ).^2 ) .* abs(u.^(2/3) - v.^(2/3));  % unbenannter Kern noname

case 11
    
    retVal   = ( ( u.^(1/3) + v.^(1/3) ).^2 ) .* abs(u.^(-1/6) - v.^(-1/6)); % Gravi like

case 12  % stokes

    retVal = 0;


end
return