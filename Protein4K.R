install.packages(c("devtools", "tidyverse", "esquisse", "ggplot2", "esquisse", "igraph"))
library(tidyverse)
library(esquisse)
library(ggplot2)
mpg
esquisse::esquisser()


# second step
#Use to download or update protein8k
devtools::install_github("SimonLiles/protein8k", build_vignettes = FALSE)

library(protein8k)

#Working with PDB Files ########################################################
#Load the data
fileName <- "data/1aieH"

my_protein <- p53_tetramerization

#Generates a static 3D plot of a protein
plot3D(my_protein)

#Animate the protein structure spinning
plot3D(my_protein, animated = TRUE, type = "p", groups = residue_name,
       image_width = 300, image_height = 300)

#Generate a model of each plane of the protein structure
plotModels(my_protein)

#Working with NCBI Metadata ####################################################
library(protein8k)

#Loading the data
#String representing the file path from the working directory to the report
jsonListName <- "data_report_small.jsonl"

#Convert JSONL into a single large data frame
report_df <- report_as_dataframe(fromJSONL(jsonListName))

##Visualizing the data #########################################################
#Make a partition graph and alluvial plot

#Prepare the raw data for subsequent visualizations
report_clean_df <- report_df
report_clean_df$geo_Location <- gsub(":.*", "", report_df$geo_Location)
report_clean_df$isolate_source <- as.character(report_clean_df$isolate_source)
report_clean_df$root <- "root"

#Load libraries
library(ggraph)
install.packages("igraph")
library(igraph)
library(tidyverse)

#Make list of nodes with frequencies for sizes
nodes <- report_clean_df %>%
  pivot_longer(cols = c(root, isolate_source, geo_Region, geo_Location),
               values_to = "label") %>%
  count(label, name = "size") %>%
  rowid_to_column("id")

#Adjust scale so that it is easier to read
nodes$size[nodes$label == "root"] <- 1
nodes$size <- log10(nodes$size)

#Move ids back by 1 so that they start at 0
nodes$id <- nodes$id - 1

#Making the edge list
#First create a list with the explicit hierarchy
hierarchy_list <- list()

for(record in 1:nrow(report_clean_df)) {
  branch <- c(report_clean_df$root[record],
              report_clean_df$isolate_source[record],
              report_clean_df$geo_Region[record],
              report_clean_df$geo_Location[record])
  hierarchy_list[[record]] <- branch
}

#Create edge list ("from" in first column / "to" in the second)
d <- do.call(rbind, hierarchy_list)
edges <- rbind(d[,1:2], d[,2:3], d[,3:4])

edges <- as.data.frame(edges)

edges <- edges %>%
  left_join(nodes, by = c("V1" = "label")) %>%
  rename(from = id)

edges <- edges %>%
  left_join(nodes, by = c("V2" = "label")) %>%
  rename(to = id)

#Clean up edges dataframe
edges <- select(edges, from, to)

#Get weights for edges data frame
edges <- group_by(edges, from, to) %>%
  count() %>%
  rename(weights = n) %>%
  ungroup()

#Make the network for the partition graph
loc_graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

#Circular partition graph using ggraph and ggplot
ggraph(loc_graph, 'partition', circular = TRUE) +
  geom_node_arc_bar(aes(fill = size, color = depth), size = 0.25) +
  geom_node_text(aes(label = label), size = 2.5) +
  theme_void() +
  theme(legend.position = "right") +
  scale_color_continuous(guide = "none") +
  scale_fill_viridis_c(direction = 1)

#Making an alluvial plot
library(ggalluvial)

#Prepare a new edge list for the alluvial plot
edge_alvl <- d
edge_alvl <- edge_alvl[,2:4]
edge_alvl <- as.data.frame(edge_alvl)

#Add weights column and rename columns to be more meaningful
edge_alvl <- group_by(edge_alvl, V1, V2, V3) %>%
  count() %>%
  ungroup() %>%
  rename(source = V1) %>%
  rename(region = V2) %>%
  rename(location = V3) %>%
  rename(weight = n)

#Remove NA values from source column
edge_alvl <- edge_alvl[complete.cases(edge_alvl$source),]

#Make weights log10 scale
edge_alvl$weight <- log10(edge_alvl$weight)

#Check the form
is_alluvia_form(edge_alvl)

#Plot with ggalluvial
#Color by source
ggplot(edge_alvl,
       aes(y = weight, axis1 = source, axis2 = region, axis3 = location)) +
  geom_alluvium(aes(fill = source), width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = FALSE,
            size = 2)

#Color by region
ggplot(edge_alvl,
       aes(y = weight, axis1 = source, axis2 = region, axis3 = location)) +
  geom_alluvium(aes(fill = region), width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = FALSE,
            size = 2)

