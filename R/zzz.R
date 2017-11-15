.onAttach <- function (lib, pkg){
  packageStartupMessage("This is biotic_cr, version ",
                        utils::packageDescription("biotic_cr", fields="Version"),
                        appendLF = TRUE)
}
