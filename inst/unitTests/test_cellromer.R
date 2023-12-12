get_pkgname <- function(){
  get0("pkgname",ifnotfound="CellRomeR")
}

find_extdata_dir <- function(){
  root <- system.file(package=get_pkgname())
  Filter(file.exists,
         c(
           file.path(root,"inst","extdata"),
           file.path(root,"extdata")))[[1]]
}

get_pkg_function <- function(funcname)
  do.call(":::", list(get_pkgname(),
                      substitute(funcname)))

## TrackMate XML import
test.loading_xml <- function(){
  get_pkg_function(import_XML)(tmXML =
                                         file.path(
                                           find_extdata_dir(),
                                           "ExampleTrackMateData.xml"),
                           dataset = "Z0_T00_C1", 
                           experiment = "ExpT", sample = "S-test",
                           condition = "Test", MigrDatObj = NULL)
  RUnit::checkTrue(TRUE)
}
