library(optparse)
library(markdown)
option_list <- list(
  make_option(c("-t", "--template"), type = "character",
              default = "../Templates/BiasEvaluationTemplate.Rmd", help = "config file path",
              metavar="character"),
  make_option(c("-d", "--dir"), type = "character",
              default = "../data/InterceptBiasSimulation/", help = "config file path",
              metavar="character"),
  make_option(c("-n", "--num"), type = "integer",
              default = 1, help = "config file path",
              metavar="integer")
)


opt_parser <- OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
if(!file.exists(paste0(opt$dir, "Sim", opt$num, "/template.Rmd"))){
  file.create(paste0(opt$dir, "Sim", opt$num, "/template.Rmd"))
}

file.copy(opt$template, paste0(opt$dir, "Sim", opt$num, "/template.Rmd"), overwrite = TRUE)
rmarkdown::render(paste0(opt$dir, "Sim", opt$num, "/template.Rmd"))
