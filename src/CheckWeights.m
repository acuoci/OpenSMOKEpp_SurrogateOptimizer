function checkWeights(optimizedProperties,doIDTLBV)
    % Put to 0 the weights of disabled properties
    global optimizer
    if optimizedProperties(1)==false, optimizer.weight_MW=0; end
    if optimizedProperties(2)==false, optimizer.weight_HC=0; end
    if optimizedProperties(3)==false, optimizer.weight_CN=0; end
    if optimizedProperties(4)==false, optimizer.weight_TSI=0; end
    if optimizedProperties(5)==false, optimizer.weight_mu=0; end
    if optimizedProperties(6)==false, optimizer.weight_YSI=0; end
    if optimizedProperties(7)==false, optimizer.weight_DC=0; end
    if optimizedProperties(8)==false, optimizer.weight_rho=0; end
    if optimizedProperties(9)==false, optimizer.weight_idt=0; end
    if optimizedProperties(10)==false,optimizer.weight_lbv=0; end
    if doIDTLBV == false
        optimizer.weight_idt = 0; optimizer.weight_lbv = 0; 
    end
end