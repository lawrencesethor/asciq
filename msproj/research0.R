

## Reading the simulation result

sim_result <- read.table("sixtables.csv", header = T, sep = ",")

## round the Coverage probability to two decimal places
NAME <- names(sim_result)
sim_result1 <- as.data.frame(matrix(NA,48,17))
## give column names
colnames(sim_result1) <- NAME
for(d in 1:17){
  if(
    NAME[d] == "Norm" ||
    NAME[d] == "Exp"  ||
    NAME[d] == "Cau"  ||
    NAME[d] == "Lap"  ||
    NAME[d] == "GEV"  ||
    NAME[d] == "NMix"
    ){
    sim_result1[,d] <- round(sim_result[,d],3)
  } else {
    sim_result1[,d] <- sim_result[,d]
  }
}

sim_result1 

write.csv(sim_result1, "sim_result1.csv")


## Copied from Hadley Wickham's R Packages and modified

tabularFun <- function(df, ...) {
  stopifnot(is.data.frame(df))
  
  align <- function(x) if (is.numeric(x)) "r" else "l"
  col_align <- vapply(df, align, character(1))
  
  cols <- lapply(df, format, ...)
  contents <- do.call("paste",
                      c(cols, list(sep = " & ", collapse = " & \n  ")))
  
  paste("\\tabular{", paste(col_align, collapse = ""), "}{\n  ",
        contents, " &" , "\n}\n", sep = "")
}

## testing the function 
cat(tabularFun(mtcars[1:5, 1:5]))

## Applying the function on my data -- sim_result1
## this helps a lot with copying the data into latex format
cat(tabularFun(sim_result1))










