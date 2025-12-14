library(data.table)

# Read PANDA output 
edges <- fread("LUSC_panda_output.txt") 

# Check the range of force values
edges_top200 <- edges[order(-abs(force))][1:200]

# Save to file
fwrite(edges_top200, "processed_data/panda_edges_LUSC_top200.txt", sep = "\t")



#Obtain differential regulatory edges 

# Load both networks (long format: tf, gene, motif, force)
luad <- fread("LUAD_panda_output.txt")
lusc <- fread("LUSC_panda_output.txt")

# Merge by TF and gene
diffRes <- merge(luad[, .(tf, gene, force_LUAD = force)],
                 lusc[, .(tf, gene, force_LUSC = force)],
                 by = c("tf", "gene"))

# Compute differential edge weight
diffRes[, delta_force := force_LUAD - force_LUSC]

# Tag which condition has the higher force
diffRes[, LUAD_higher := delta_force > 0]

# Number of differential edges
nrow(diffRes)

# Top 200 stronger in LUAD
n <- 200
diffResLUAD <- diffRes[LUAD_higher == TRUE][order(-delta_force)][1:n]

# Top 200 stronger in LUSC
diffResLUSC <- diffRes[LUAD_higher == FALSE][order(delta_force)][1:n]

# Combine for Cytoscape
diffRes_vis <- rbind(diffResLUAD, diffResLUSC)
fwrite(diffRes_vis, "processed_data/diffRes_LUAD_vs_LUSC.txt", sep="\t")