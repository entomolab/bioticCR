.onAttach <- function (lib, pkg){
  packageStartupMessage("This is bioticCR, version ",
                        utils::packageDescription("bioticCR", fields="Version"),
                        appendLF = TRUE)
}
