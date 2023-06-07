.sg <- function(probs)
{
  #input: probs = extremile order in (0,1)
  #output: power $s(probs)$ in the definition of the measure $K_{probs}$ [see daouia & Gijbels (Extremiles)]
  
  return(log(1/2)/log(1-probs));
}

.Kg <- function(probs,t)
{
  #input: probs = extremile order in (0,1)
  #           t = real number in the support [0,1] of the $K_{probs}$
  #output: measure $K_{probs}(t)$ [see daouia & Gijbels (Extremiles)]
  res=((1-(1-t)^(.sg(probs)))*(0<probs)*(probs<=0.5))+((t^(.sg(1-probs)))*(0.5<probs)*(probs<1));     
  return(res);
}

extremile <- function(ytab,probs)
{
  #(Extremiles: by Daouia.A, May 4, 2015)
  #inputs:
  #ytab = sample of nx1 obs
  # probs = extremile order in (0,1)
  #output: $probs$th sample extremile          
  
  sytab=sort(ytab);
  res=0;
  n=length(ytab);
  for (i in 1:n)
  {res = res + ((.Kg(probs,i/n)-.Kg(probs,(i-1)/n))*sytab[i]);}
  return(res);
}