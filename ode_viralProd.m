function y = ode_viralProd(t,x,k)
kUptake         = k(1);
kEscape         = k(2);
kNuclearEntry   = k(3);
kPlasmidDeg     = k(4);

kRepSyn         = k(5);
kCapSyn         = k(6);
kProteinDegRep  = k(7);
kProteinDegCap  = k(8);
kDNArep         = k(9);
kBindRCplasmid  = k(10);
kAssembly       = k(11);
kBindCapsid     = k(12);
kPack           = k(13);
kSecrete        = k(14);

mu = calc_growth_rate(t); 
%% Delivery

dx1 = [ -k(1)*x(1);...                                                      % 1. pRC extracellular
        k(1)*x(1)  -   (k(2) + k(4) + mu)*x(2);...                          % 2. pRC endosomal
        k(2)*x(2)  -   (k(3) + k(4) + mu)*x(3);...                          % 3. pRC cytosol
        k(3)*x(3)  -   mu*x(4) - k(4)*x(4) - kBindRCplasmid*x(4)*x(13)];    % 4. pRC nucleus

dx2 = [ -k(1)*x(5);...                                                      % 5. pVector extracellular
        k(1)*x(5)  -   (k(2) + k(4) + mu)*x(6);...                          % 6. pVector endosomal
        k(2)*x(6)  -   (k(3) + k(4) + mu)*x(7);...                          % 7. pVector cytosol
        k(3)*x(7)  -   mu*x(8) - k(4)*x(8)];                                % 8. pVector nucleus                        
    
dx3 = [ -k(1)*x(9);...                                                      % 9. pHelper extracellular
        k(1)*x(9)  -   (k(2) + k(4) + mu)*x(10);...                         % 10. pHelper endosomal
        k(2)*x(10)  -  (k(3) + k(4) + mu)*x(11);...                         % 11. pHelper cytosol
        k(3)*x(11)  -   mu*x(12)  - k(4)*x(12)];                            % 12. pHelper nucleus        
    
%% Viral Replication
rRepSyn     = kRepSyn*x(4)*x(12);
rCapSyn     = kCapSyn*x(4)*x(12);
rDNArep     = kDNArep*x(13)*x(8)*x(12);
rBindRC     = kBindRCplasmid*x(4)*x(13);
rAssembly   = kAssembly*x(14);
rBindCapsid = kBindCapsid*x(13)*x(16);
rPack       = kPack*x(21)*x(15);


dVP = [ rRepSyn - rBindRC  - rBindCapsid + rPack - (kProteinDegRep + mu)*x(13)  ;...            % 13. repProtein

        rCapSyn - 60*rAssembly - (kProteinDegCap + mu)*x(14);...            % 14. cap Protein

        rDNArep - rPack - mu*x(15)                         ;...             % 15. vDNA

        rAssembly -   (kSecrete + mu)*x(16)    - rBindCapsid       ;...  	% 16. empty cap nuc

        rPack - (kSecrete + mu)*x(17)                           ;...        % 17. full cap nuc

        kSecrete*x(16) - mu*x(18)                               ;...        % 18. empty capsid cyto

        kSecrete*x(17) - mu*x(19)                               ;...        % 19. full capsid cyto
        
        rBindRC                                                ;...         % 20. Rep bound RC plasmid complex
        
        rBindCapsid - rPack                                     ];          % 21. Rep-bound empty capsid nucleus                      
        

y = [dx1; dx2; dx3; dVP];
end
