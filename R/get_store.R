#' Get the targets store on various machines
#'
#' @returns
#' @export
#'
#' @examples
get_store <- function(){
  if(dir.exists("D:/Data_For_Remi/narw_ais")) {
    # on the HDD
    store = "D:/Data_For_Remi/narw_ais"
  } else if(dir.exists("//ci-WPNSBIO9039519-smb-1.mar.dfo-mpo.ca/ocean_data/SPA/NARW_AIS")) {
    # on the NAS from windows
    store = "//ci-WPNSBIO9039519-smb-1.mar.dfo-mpo.ca/ocean_data/SPA/NARW_AIS"
  } else if(dir.exists("~/ocean_data")) {
    # on the NAS from linux
    store = normalizePath("~/ocean_data/SPA/NARW_AIS")
  } else if (dir.exists("/srv/sambashare/NARW/SPA/NARW_AIS")) {
    # on a linux machine on the sambashare
    store = "/srv/sambashare/NARW"
  } else if (
    dir.exists(
      "//wpnsbio9039519.mar.dfo-mpo.ca/sambashare/NARW"
    )
  ) {
    # on a windows machine on the sambashare
    store <- "//wpnsbio9039519.mar.dfo-mpo.ca/sambashare/NARW"
  } else {
    store <-  getwd()
  }

  return(store)
}
