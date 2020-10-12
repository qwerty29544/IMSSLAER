# get_packages ------------------------------------------------------------

#' Get Packages function
#' @description For downloading and load to projects GPU torch functions
#' @return
#' @export
#'
#' @examples IMSSLAER::get_packages()
get_packages <- function() {
    if ("torch" %in% rownames(installed.packages()) == FALSE) {
        cat("\nTake it easy, torch is installing on your device now ...\n")
        install.packages("torch")
    }
    library(torch)
    cat("\nDONE!\n")
}
