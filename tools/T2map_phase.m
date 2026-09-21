function T2map = T2map_phase(img_phase, LUTs, mode)

if strcmp(mode,'phase')
    
    %disp('Phase Mode')

    T2map = 0*img_phase;
    for ii=1:size(img_phase,1)
        for jj=1:size(img_phase,2)

            [t,ind] = min(abs(LUTs.theta-img_phase(ii,jj)));
            T2map(ii,jj)=LUTs.index(1,ind,2);
        end
    end

else
    
    %disp('Magnitude Mode')
 
    T2map = 0*img_phase;
    for ii=1:size(img_phase,1)
        for jj=1:size(img_phase,2)

            [t,ind] = min(abs(LUTs.eta-img_phase(ii,jj)));
            T2map(ii,jj)=LUTs.index(1,ind,2);
        end
    end
    
end

