

## ---- build_autocross.R (run from ~/across/autocross) ----

## 1. Unload previous version
if ("package:autocross" %in% search())
  detach("package:autocross", unload = TRUE, character.only = TRUE)

if ("autocross" %in% loadedNamespaces())
  unloadNamespace("autocross")

## 2. Clean possible stale DLL (macOS/Linux use .so, Windows uses .dll)
so_file <- file.path("src", ifelse(.Platform$OS.type == "windows",
                                   "autocross.dll", "autocross.so"))
if (file.exists(so_file))
  file.remove(so_file)

## 3. Build and document
library(devtools)
document()
tarball <- build()      # returns full path to .tar.gz (Mac/Linux) or .zip (Win)

## 4. Make sure the local library exists
dir.create(".r-lib", showWarnings = FALSE)

## 5. Install into project-local library
install.packages(tarball, repos = NULL, type = "source", lib = ".r-lib")

## 6. Load from that library
library(autocross, lib.loc = ".r-lib")
