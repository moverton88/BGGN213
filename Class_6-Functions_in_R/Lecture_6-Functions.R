
### Working with bio3d package
# First, we installed the bio3d package using install.packages(), then loaded it using library()
# The below code is copied from the BGGN213 worksheet. It downloads .pdb protein structure files, 
# trims them to only include chain A, and plots the atomic, 3d structure of each. However it does contain
# errors, which I have corrected

library(bio3d)

s1 <- read.pdb("4AKE")  # kinase with drug
s2 <- read.pdb("1AKE")  # kinase no drug
s3 <- read.pdb("1E4Y")  # kinase with drug

s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s3, chain="A", elety="CA")

s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b

plotb3(s1.b, sse=s1.chainA, typ="l", col="red", ylab="Bfactor")
plotb3(s2.b, sse=s2.chainA, typ="l", col="blue", ylab="Bfactor")
plotb3(s3.b, sse=s3.chainA, typ="l", ylab="Bfactor")
points(s1.b, typ="l", col="red")


# The above code works fine, but lets see if we can improve it

# Inputs include a list of PDB codes, which the function uses to download protein information from the 
# PDB database, as well as the protein chain to focus on and the elety
PDB_codes <- list("4AKE", "1AKE", "1E4Y")
protein_chain = "A"
protein_elety = "CA"

# The Bfactor_plot function takes in the codes, chain, and elety inputs, downloads the protein data, and
# organizes the data into a list. 
Bfactor_plot <- function(codes, chain="A", elety="CA") {
  p_list <- lapply(codes, function(p) read.pdb(p))
  # The list is manipulated with the lapply function to parse the protein
  # data based on the chain and elety desired.
  p_list <- lapply(p_list, function(p) trim.pdb(p, chain=chain, elety=elety))
  # A color palette is created to make the plot lines distinguishable
  col_pal <- data.frame(num=c(1:5), col=c("black", "red", "blue", "green", "purple"), stringsAsFactors = F)
  # Protein Bfactor data is plotted, with each color representing a different protein, with the output
  # a plot of Bfactor and secondary structure along the length of the protein
  plotb3(p_list[[1]]$atom$b, sse=p_list[[1]], typ="l", ylab="Bfactor")
  for (i in 2:length(p_list)){
    points(p_list[[i]]$atom$b, typ="l", col=col_pal[i,2])
  }
}


Bfactor_plot(codes=PDB_codes)






