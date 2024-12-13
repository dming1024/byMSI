
#byMSI 

p_msi=0.05
posterior_probability<-function(len_indels,lambda_MSI,lambda_MSS,p_msi){
  p_mss=1-p_msi
  (p_msi * dpois(len_indels,lambda_MSI)) / ( p_msi * dpois(len_indels,lambda_MSI)  + p_mss * dpois(len_indels,lambda_MSS)) 
}

MSIp<-function(sequence,allele_sequence,prevalence=0.05){
  #prevalence, 0.05,by default
  lambda_msi=3.58 #locus specific
  lambda_mss=1.22 #locus specific
  
  ref_alternation=abs(sequence)
  p_ref=posterior_probability(ref_alternation,lambda_msi,lambda_mss,prevalence)
  if(is.null(allele_sequence)){
    p_ = p_ref^2
  }else{
    allele_alternation=abs(allele_sequence)
    p_allel=posterior_probability(allele_alternation,lambda_msi,lambda_mss,prevalence)
    p_=sqrt(p_ref*p_allel) 
  }
  
  return(p_)
}

#input the length of indels in chr3:124951171
p_ = MSIp(10,5)