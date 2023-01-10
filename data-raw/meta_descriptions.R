## code to prepare `meta_descriptions` dataset goes here
readme <- as.data.frame(read_xlsx("T:/Research/Aqualog_Data/0_Codes and Files/metatable_manual.xlsx", sheet = "description"))

usethis::use_data(meta_descriptions, overwrite = TRUE, internal=TRUE)
