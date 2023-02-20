# Thank you Jeff Oliver for your code (https://github.com/jcoliver/biodiversity-sdm-lesson)

########################### 1. Run the setup code below ####################################
# This installs libraries, and downloads climate data from bioclim (https://www.worldclim.org/data/bioclim.html)

source(file = "src/setup.R")

############### 2. In the "src" directory, copy the contents of "run-sdm-single.R" ##############
# into a new file (still in 'src') called <species>-sdm-single.R (Rename <species> to your species)



# 3. In this new file, edit lines 14 & 15, changing MY_SPECIES to your species.




# 4. Below, write your spocc/gbif query, 
# and then use the "$" notation to create a variable targeting the data set. 


# Phidippus johnsoni (spiders)
# https://www.gbif.org/species/5173184
# limiting to 2020 just to keep the volume down
#johnsoni<-occ(query='Phidippus johnsoni', from="gbif", gbifopts = list(year="2020"))

#navigate object
#johnsoniData<-johnsoni$gbif$data$Phidippus_johnsoni



################ 5. Save to CSV #####################

# first, ensure all data is character data
#df <- apply(df,2,as.character)

#johnsoniData <- apply(johnsoniData,2,as.character)

# use write.csv to write the data frame to 'data' directory
# make sure the file name matches what you indicated in step 3 on line 14
#write.csv(johnsoniData, "data/johnsoni.csv")


# 6. Use the source() command to run the file you created in step 2 ############

source("src/rana-run-sdm-single.R")
# note: large data set could crash things....

#source("path/to/file.R")
