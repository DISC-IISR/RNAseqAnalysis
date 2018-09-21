#'@title loading data
#'
#'@description copies the data needed to run the example form the system folder to the current working directory
#'@export
#'@seealso \code{\link{file.copy}}
#'@examples
#'
#'\dontrun{
#'import_dataset()
#'}
#'
import_dataset <- function()
{
  currentfiles <- system.file("dataset", package="newtest")
  dir.create("dataset")
  new_folder <- "./dataset"
  list_of_files <- list.files(currentfiles)
  file.copy(file.path(currentfiles,list_of_files), new_folder)
}
