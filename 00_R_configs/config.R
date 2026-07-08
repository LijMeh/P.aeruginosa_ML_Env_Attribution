# Auto-detects computer and sets appropriate data path
# This is a generic template - replace computer names and paths with your own

computer_name <- Sys.info()["nodename"]

if (computer_name == "personal-desktop") { # Home desktop
  DATA_PATH <- "C:/Users/username/Cloud Storage/project_folder/data"

} else if (computer_name == "personal-laptop") { # Personal laptop
  DATA_PATH <- "C:/Users/username/Cloud Storage/project_folder/data"
  
} else if (computer_name == "work-computer") { # Work computer
  DATA_PATH <- "C:/Users/username/Cloud Storage/project_folder/data"

} else if (computer_name == "remote-server.domain.com") { # Remote server
  DATA_PATH <- "/home/username/data/project_folder/"

} else {
  stop(paste("Unknown computer:", computer_name, 
             "\nPlease add this computer's data path to config.R"))
}

# Helper function to access data
data_here <- function(...) {
  file.path(DATA_PATH, ...)
}