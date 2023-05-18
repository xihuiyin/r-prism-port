##GRAPHPAD PORT##

library(dplyr)
library(stringr)

#import file and create a dataframe df from csv file
df <- read.csv("4.28.23 plate freezing expt -norm.csv", stringsAsFactors = F, encoding="UTF-8")

#retrieve antigen list
ag_list <- df$Ag
ag_list <- lapply(ag_list, function(x) sort(x))


#remove duplicates in ag_list
ag_list <- ag_list[!duplicated(ag_list)]

# combine the list elements into a vector
ag_list <- unlist(ag_list)
is.list(ag_list)
is.factor(ag_list)
ag_list <- as.factor(ag_list)
ag_list <- sort(ag_list)

#generate a vector (1 to n) of the n # of antigens in ag_list
num_ag <- seq_along(ag_list)

#split df into list of dataframes for each antigen
df_list <- split(df, df$Ag)
df_list <- as.array(df_list)

#pulls out list of categories and samples (which is the same across antigens)
test <- df_list[[1]]
cat_list <- test[,2]
sample_list <- test[,5]

#remove duplicates in the lists
cat_list <- cat_list[!duplicated(cat_list)]
sample_list <- sample_list[!duplicated(sample_list)]

#convert character vector into factor
cat_list <- as.factor(cat_list) 
sample_list <- as.factor(sample_list)
test$Category <- as.factor(test$Category)

#check if these objects are factors, if not, convert to factors. They need to be factors so we can pull the levels from them
is.factor(cat_list) 
is.factor(sample_list)
is.factor(test$Category)
is.list(df_list)

# define an empty list to store the adapted arrays for graphpad
final.ag.list <- vector("list", length = length(df_list))

names(final.ag.list) <- names(df_list)

for (k in 1:length(df_list)) { #goes thru each dataframe for antigen in dm_list @iteration k
  final.ag.list[[k]] <- matrix(nrow=length(sample_list), ncol=length(cat_list)) #creates matrix for each antigen @ iteration k
  ag.holder <- df_list[[k]] #temporarily stores values for iteration k
  ag.holder$Category <- as.factor(ag.holder$Category) 
  for (i in 1:nrow(ag.holder)) { # loop through each row (sample) for specific antigen
    cat_match <- match(ag.holder$Category[[i]], cat_list)  # find the index of cat_list that matches the sample category
    serum_match <- match(ag.holder$Serum[[i]], sample_list)
    if (!is.na(cat_match)) {  # if the category was found in cat_list
      #print(cat_match)
      #print(k)
      #print(i)
      final.ag.list[[k]][serum_match,cat_match] <- ag.holder$MFI[i] # insert the MFI value of sample into corresponding column
    }
  }
  mat <- final.ag.list[[k]]
  rownames(mat) <- sample_list
  colnames(mat) <- cat_list
  final.ag.list[[k]] <- mat
  
  #replacing illegal characters (i.e. > if present in raw data)
  rownames(final.ag.list[[k]]) <- sub(">", "&gt;", rownames(final.ag.list[[k]])) # Replace ">" with "&gt;"
  colnames(final.ag.list[[k]]) <- sub(">", "&gt;", colnames(final.ag.list[[k]]))
  
}
str(df_list[[1]])


###time to make my txt file to be inserted into a .pzf file!##
##First, create a file called pzf port.txt##

#edit number of tables
table_list <- c('<TableSequence>')
print(table_list)
for (i in seq_along(ag_list)) {
  table_list <- c(table_list, paste0('<Ref ID = "Table',i,'"/>'))
}

print(table_list)

table_list <- c(table_list, '</TableSequence>') 

#now for the table itself
table_content <- character()

#for each antigen, print the list of samples within each antigen block
for (i in seq_along(ag_list)) {
  cat_names <- character()
  sample_str <- paste0("<d>", sample_list, "</d>\n", collapse = "") #pastes a list of the samples in <d> tags
  cat_names <- paste0("<d>",colnames(final.ag.list[[i]]),"</d>", collapse="\n")
  
  mfi_table_content <- character() 
  
  for (k in seq_len(ncol(final.ag.list[[i]]))) { #goes through mfi values under each category for current iteration (antigen)
    mfi_values <- final.ag.list[[i]][,k]
    mfi_str <- ""
    for (j in seq_along(mfi_values)) {
      if (is.na(mfi_values[j])) {
        mfi_str <- paste0(mfi_str, "<d/>\n") #if there is no value present (patient belongs to other cats), paste <d/>
      } else {
        mfi_str <- paste0(mfi_str, "<d>", mfi_values[j], "</d>\n") #if value is present, add <d> tags
      }
    }
    mfi_str <- paste0(mfi_str, '</Subcolumn>\n','</YColumn>\n') #add closing tags for the category
    mfi_table_content <- paste0(mfi_table_content, '<YColumn Width="125" Decimals="2" Subcolumns="1">\n', 
                                '<Title>', cat_list[k], '</Title>\n',
                                '<Subcolumn>\n',
                                mfi_str)
  }
  
  mfi_table_str <- paste(mfi_table_content, collapse = "")
  
  ag_table_content <- paste0('<Table ID="Table',i,'" YFormat="replicates" Replicates="1" TableType="TwoWay" EVFormat="AsteriskAfterNumber">\n',
                             '<Title>', ag_list[i], '</Title>\n',
                             '<RowTitlesColumn Width="151">\n',
                             '<Subcolumn>\n', sample_str, '</Subcolumn>\n',
                             '</RowTitlesColumn>\n',
                             mfi_table_str,'</Table>')
  
  table_content <- c(table_content, ag_table_content)
}

final.output <- c(table_list,table_content)
write.table(final.output, "table.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)