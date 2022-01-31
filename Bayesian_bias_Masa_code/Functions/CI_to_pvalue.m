function [pvalue_versus_original, pvalue_versus_shuffle] = CI_to_pvalue(log_odd,AUC,AUC_UCI,AUC_LCI)
samplesize = 1000;% Because bootstrapped 1000 times
for epoch = 1:3
    k = 1;
    for test = 2:5
        SD1= ((abs(AUC_UCI(epoch,1)-AUC(epoch,1)) + abs(AUC_LCI(epoch,1)-AUC(epoch,1)))/2)/tinv(0.975,samplesize-1)*sqrt(samplesize);
        SD2= ((abs(AUC_UCI(epoch,test)-AUC(epoch,test)) + abs(AUC_LCI(epoch,test)-AUC(epoch,test)))/2)/tinv(0.975,samplesize-1)*sqrt(samplesize);
        t = (AUC(epoch,1)-AUC(epoch,test))/sqrt( (SD1^2/samplesize) + (SD1^1/samplesize));
        pvalue_versus_original(epoch,k) = tcdf(abs(t),samplesize-1,'upper');
        k = k+1;
    end
    
    k = 1;
    for test = [1,2,3,5]
        SD1= ((abs(AUC_UCI(epoch,test)-AUC(epoch,test)) + abs(AUC_LCI(epoch,test)-AUC(epoch,test)))/2)/tinv(0.975,samplesize-1)*sqrt(samplesize);
        SD2= ((abs(AUC_UCI(epoch,4)-AUC(epoch,4)) + abs(AUC_LCI(epoch,4)-AUC(epoch,4)))/2)/tinv(0.975,samplesize-1)*sqrt(samplesize);
        t = (AUC(epoch,test)-AUC(epoch,4))/sqrt( (SD1^2/samplesize) + (SD1^1/samplesize));
        pvalue_versus_shuffle(epoch,k) = tcdf(abs(t),samplesize-1,'upper');
        k = k+1;
    end
end

end

