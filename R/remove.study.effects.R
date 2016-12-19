remove.study.effects <- function(effectsize, se, meta.effectsize, meta.se) {
  #' Function to remove an individual study's effect contribution to a 
  #' fixed-effect meta-analysis summary statistics
  #' @keywords internal
  
  if(cor(effectsize, meta.effectsize) < 0) warning("Correlation of effectsize and meta.effectsize is negative")
  
  # prec stands for 'precision'
  prec.meta <- 1/meta.se^2
  prec <- 1/se^2
  prec.new <- prec.meta - prec
  problems <- prec.new <= 0
  prec.new[problems] <- 0
  
  eff.new <- (prec.meta * meta.effectsize - prec * effectsize)/prec.new
  eff.new[is.infinite(eff.new)] <- 0
  
  return(list(effectsize=eff.new, se=1/sqrt(prec.new), problems=problems))
  
}