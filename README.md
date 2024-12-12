
# byMSI 

Bayesian Inference of Microsatellite Instability Status Based on Indel Variations in Tumor Genomes

## optional: update parameters

+ $\lambda_{msi-h}$
+ $\lambda_{mss}$
+ $prevalece_msi-h$,  ($p_{msi-h}$)
+ $prevalence_mss$, ($p_{mss}$)

## Usage

```R

p_msi=0.05
posterior_probability<-function(len_indels,lambda_MSI,lambda_MSS,p_msi){
  p_mss=1-p_msi
  (p_msi * dpois(len_indels,lambda_MSI)) / ( p_msi * dpois(len_indels,lambda_MSI)  + p_mss * dpois(len_indels,lambda_MSS)) 
}

MSIp<-function(sequence,allele_sequence,prevalence=0.05){
  #prevalence, 0.05
  lambda_msi=3.58
  lambda_mss=1.22
  
  ref_alternation=abs(sequence)
  p_ref=posterior_probability(ref_alternation,lambda_msi,lambda_mss,prevalence)
  if(is.null(allele_sequence)){
    p_ = p_ref^2
  }else{
    allele_alternation=abs(allele_sequence)
    p_allel=posterior_probability(allele_alternation,lambda_msi,lambda_mss,prevalence)
    p_=sqrt(p_ref*p_allel) 
    #p_=max(c(p_ref,p_allel))
  }
  
  return(p_)
}

#input the length of indels in chr3:124951171
p_ = MSIp(10,5)
```

+ The null hypothesis ($H_0$) is rejected if the `p_` < 0.05, in which case the sample is considered to not be MSI-H
+ Alternative hypothesis ($H_1$) assumes the sample is not MSI-H