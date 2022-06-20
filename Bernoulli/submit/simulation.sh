#!/bin/bash
Rscript -e "rmarkdown::render('simulation.Rmd', 
                  params = list(confi_index = '$confi_index'),
                  output_file = 'simulation.html')"