library(reticulate)

# Install uci_datasets if needed
# reticulate::py_install("uci_datasets")

uci <- import("uci_datasets")

datasets <- c("wine", "energy", "parkinsons", "stock", "pumadyn32nm", 
              "housing", "forest", "bike", "solar", "gas")

dir.create("./data/uci", recursive = TRUE, showWarnings = FALSE)

for (name in datasets) {
  tryCatch({
    ds    <- uci$Dataset(name)
    xy    <- ds$get_split(split = 0L)  # returns (x, y) as numpy arrays
    x     <- as.data.frame(xy[[1]])
    y     <- as.data.frame(xy[[2]])
    df    <- cbind(x, y)
    saveRDS(df, file = paste0("./data/uci/", name, ".RDS"))
    cat("Saved:", name, "\n")
  }, error = function(e) {
    cat("Failed:", name, "-", conditionMessage(e), "\n")
  })
}