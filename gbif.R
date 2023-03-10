#Querying, cleaning, and mapping gbif data
#2/7/23



#spocc: interface to SPecies OCCurrence data sources

#Install packages and load libraries:

install.packages("spocc")
install.packages("tidyverse") #includes ggplot
library(spocc)
library(tidyverse)


#Commenting: it's good
# here is my comment


#Functions are built-in commands that often take "arguments",
# or variables/data you pass to the function.

sqrt(10)


#Occurrence Data: what is it?


# We can directly query occurrence data in R with the "occ" function, which is part of the spocc package.
# Take a look at the documentation here:
# https://www.rdocumentation.org/packages/spocc/versions/0.1.0/topics/occ

# gbifopts: run occ_options('gbif') in console to see possibilities()


#In groups: write a query for Rana draytonii from GBIF. 
#Any country. Any time period. Limit 5000.

myQuery<-occ(query="Rana draytonii", from="gbif", limit=4000)

myQuery



# Drill down to get the data using "$", and show from Env window

rana<-myQuery$gbif$data$Rana_draytonii



### Let's initially plot the data on a map.

ggplot()+
  geom_point(data=rana, mapping=aes(x=longitude, y=latitude), show.legend=FALSE)+
  labs(title="Species occurrences of R. draytonii", x="longitude", y="latitude")

wrld<-ggplot2::map_data("world")


# say we want to add country lines

ggplot()+
  geom_polygon(data=wrld, mapping=aes(x=long, y=lat, group=group), fill="grey75", colour="grey60")+
  geom_point(data=rana, mapping=aes(x=longitude, y=latitude), show.legend=FALSE)+
  labs(title="Species occurrences of R. draytonii", x="longitude", y="latitude")



# 5 minute break!


#cleaning data: what do we need to do based on the reading for today? 
### remove outliers (outside of normal sp range), remove duplicates, deal with NA values

# first remove outliers
noAfricaPoints <- rana %>% filter(longitude < 0)

# deal with NA values
noNAPoints <- noAfricaPoints %>% filter(latitude != "NA", longitude != "NA")

# remove duplicates
noDuplicates <- noNAPoints %>% mutate(location = paste(latitude, longitude, dateIdentified, sep = "/")) %>%
  distinct(location, .keep_all = TRUE)

cleanRana <- rana %>% 
  filter(longitude < 0) %>% 
  filter(latitude != "NA", longitude != "NA") %>% 
  mutate(location = paste(latitude, longitude, dateIdentified, sep = "/")) %>%
  distinct(location, .keep_all = TRUE)

### one more column to check of your data: is occurence status == PRESENT ? 

# set x and y limits
xmax <- max(cleanRana$longitude)
xmin <- min(cleanRana$longitude)
ymax <- max(cleanRana$latitude)
ymin <- min(cleanRana$latitude)

# re-do map code with cleaned data and x/y limits

ggplot()+
  geom_polygon(data=wrld, mapping=aes(x=long, y=lat, group=group), fill="grey75", colour="grey60")+
  geom_point(data=cleanRana, mapping=aes(x=longitude, y=latitude), show.legend=FALSE)+
  labs(title="Species occurrences of R. draytonii from 1938 - 2023", x="longitude", y="latitude") +
  coord_fixed(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  scale_size_area() +
  borders("state")

# make it better by adding a date range in the title
range(na.omit(as.Date(cleanRana$dateIdentified)))


# write occurence data to csv

cleanRana <- cleanRana %>% 
  select(longitude, latitude, eventDate)
write.csv(cleanRana, "ranaData.csv")
str(cleanRana)


