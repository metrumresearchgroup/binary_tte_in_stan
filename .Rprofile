options(repos = c(
  MPN = "https://mpn.metworx.com/snapshots/stable/2021-02-01/", 
  # a value must be set to CRAN or R will complain, so we'll point both to MPN
  CRAN = "https://mpn.metworx.com/snapshots/stable/2021-02-01/"
  ),
  # set some bbr opinionated defaults, these won't impact users who don't use bbr
  'bbr.bbi_exe_path' = file.path(getwd(), "bin", "bbi")
)

# source after setting the repos to make sure renv will see those repo versions
source("renv/activate.R")

if(interactive()){
  message("repos set to: \n\t", paste0(unique(getOption('repos')), collapse = "\n\t"))
  message("library paths set to: \n\t", paste0(.libPaths(), collapse = "\n\t"))
  
  local({
    bbi_path <- getOption('bbr.bbi_exe_path')
    
    if (!file.exists(bbi_path)) {
      warning(
        sprintf("bbi not found at path `%s` ", bbi_path),
        "either run `bbr::use_bbi()` to install bbi or ",
        "check what you entered in the `bbr.bbi_exe_path` option and make sure it is correct.",
        call. = FALSE
      )
    }
  })
}
