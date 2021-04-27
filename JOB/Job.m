function [F] = Job(AD, AM, Y,alpha,beta,gama)
% AD£ºthe disease-disease similarity matrix
% AM: the miRNA-miRNA similarity matrix
% Y : the ground truth (the known disease-miRNA associations)
% F : the result predicted by our method

    nd = size(AD,1); % the size of disease
    nm = size(AM,1); % the size of miRNA

    F = Y ;          % F is the label matrix need to learn
    F_old = zeros(size(Y));

    NITER = 40;

   
    
    thresh = 10^-9;      % Iterative terminating condition

    for iter = 1:NITER
       SD = zeros(size(AD)); 
       SM = zeros(size(AM)); 
            
       %fixed SM and F update SD
       distd = L2_distance_1(F',F');
       AD0 = AD-diag(diag(AD));
        for i=1:nd
            ad0 = AD0(i,:);
            %idxd0 = 1:nd;
            idxd0 = find(ad0>0);
            adi = ad0(idxd0);
            dai = distd(i,idxd0);
            ad = adi-0.5*alpha*dai; 
            SD(i,idxd0) = EProjSimplex_new(ad);
        end;


        %fixed SD and F update SM
        AM0 = AM-diag(diag(AM));
        distm = L2_distance_1(F,F);
        for i=1:nm
            am0 = AM0(i,:);
            %idxm0 = 1:nm;
            idxm0 = find(am0>0);
            ami = am0(idxm0);
            mai = distm(i,idxm0);
            am = ami-0.5*beta*mai; 
            SM(i,idxm0) = EProjSimplex_new(am);
        end;


        %fixed SD and SM update F
        AD1 = SD-diag(diag(SD));
        AD20 = (AD1+AD1')/2;
        DD20 = diag(sum(AD20));             
        %LD1 = DD20 - AD20;
        idD = eye(size(AD20, 1));
        NDD20 = diag(1 ./ diag(sqrtm(DD20)));
        LD1 = idD - NDD20 * AD20 * NDD20;
        

        AM1 = SM-diag(diag(SM));
        AM20 = (AM1+AM1')/2;
        DM20 = diag(sum(AM20));             
        %LM1 = DM20 - AM20;
        idM = eye(size(AM20, 1));
        NMM20 = diag(1 ./ diag(sqrtm(DM20)));
        LM1 = idM - NMM20 * AM20 * NMM20;

        F = sylvester(2*beta*LD1+gama*eye(nd),2*alpha*LM1,gama*Y);
       
        diff = abs(sum(sum((abs(F)-abs(F_old)))));
        if diff < thresh
            break ; 
        end
        F_old = F ;
        
        %disp(iter)
        %disp(diff)
    end

end


   