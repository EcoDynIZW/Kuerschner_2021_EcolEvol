
# title: "source_lib"
# author: "Tobias Kuerschner"
# date: "2 December 2019"


### Packages

for (
      pckg in 
        c
          (
             'viridis','scico','tidyverse','ggforce','svglite','egg'
          )
    ) 
  { 
  if (!require(pckg, character.only=T)) install.packages(pckg, dependencies = TRUE)
  require(pckg, character.only=T)
  }


