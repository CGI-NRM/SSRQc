load_data <- function(file_path, sheet = 1, na_strings = c("NA", "-99", "0", "000", "No peaks in locus
", "No peaks"
)) {
  # Helper function read facilitate data import of different file formats
  # Returns a dataframe of the imported data
  #
  # file_path: Name of the file to imported, including path
  # sheet: Which sheet to read data from in excel and ods file formats
  # na_strings: Which entries should be converted to NA values on import
  #
  # NB! The function uses file endings to guess file formats. Will not work if file ending
  # is not consistent with actual file type
    if (endsWith(file_path, ".xls") | endsWith(file_path, ".xlsx")) {
        raw_data <- readxl::read_excel(path = file_path,
                                       col_names = TRUE,
                                       na = na_strings,
                                       sheet = sheet)
    } else if (endsWith(file_path, ".ods")) {
        raw_data <- readODS::read_ods(path = file_path,
                                      col_names = TRUE,
                                      na = na_strings,
                                      sheet = sheet)
    } else if (endsWith(file_path, ".tsv") |
               endsWith(file_path, ".tdf") |
               endsWith(file_path, ".txt")) {
        raw_data <- read.table(file = file_path,
                               header = TRUE,
                               na.strings = na_strings,
                               sep = "\t",
                               stringsAsFactors = FALSE)
    } else {
        raw_data <- read.table(file = file_path,
                               header = TRUE,
                               na.strings = na_strings,
                               sep = ",",
                               stringsAsFactors = FALSE)
    }
     data.frame(raw_data)
}
